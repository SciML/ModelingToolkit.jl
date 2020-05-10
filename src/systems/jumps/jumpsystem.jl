JumpType = Union{VariableRateJump, ConstantRateJump, MassActionJump}

struct JumpSystem <: AbstractSystem
    eqs::Vector{JumpType}
    iv::Variable
    states::Vector{Variable}
    ps::Vector{Variable}
    name::Symbol
    systems::Vector{JumpSystem}
end

function JumpSystem(eqs, iv, states, ps; systems = JumpSystem[],
                                          name = gensym(:JumpSystem))
    JumpSystem(eqs, iv, convert.(Variable, states), convert.(Variable, ps), name, systems)
end



generate_rate_function(js, rate) = build_function(rate, states(js), parameters(js),
                                                  independent_variable(js),
                                                  expression=Val{false})

generate_affect_function(js, affect, outputidxs) = build_function(affect, states(js),
                                                      parameters(js),
                                                      independent_variable(js),
                                                      expression=Val{false},
                                                      headerfun=add_integrator_header,
                                                      outputidxs=outputidxs)[2]
function assemble_vrj(js, vrj, statetoid)
    rate   = generate_rate_function(js, vrj.rate)
    outputvars = (convert(Variable,affect.lhs) for affect in vrj.affect!)
    outputidxs = ((statetoid[var] for var in outputvars)...,)
    affect = generate_affect_function(js, vrj.affect!, outputidxs)
    VariableRateJump(rate, affect)
end

function assemble_crj(js, crj, statetoid)
    rate   = generate_rate_function(js, crj.rate)
    outputvars = (convert(Variable,affect.lhs) for affect in crj.affect!)
    outputidxs = ((statetoid[var] for var in outputvars)...,)
    affect = generate_affect_function(js, crj.affect!, outputidxs)
    ConstantRateJump(rate, affect)
end

function assemble_maj(js, maj::MassActionJump{U,Vector{Pair{V,W}},Vector{Pair{V2,W2}}},
                      statetoid, parammap) where {U,V,W,V2,W2}
    sr = maj.scaled_rates
    if sr isa Operation
        pval = substitute(sr,parammap)
    elseif sr isa Variable
        pval = Dict(parammap)[sr()]
    else
        pval = maj.scaled_rates
    end
    
    rs = Vector{Pair{Int,W}}()
    for (spec,stoich) in maj.reactant_stoch
        if iszero(spec)
            push!(rs, 0 => stoich)
        else
            push!(rs, statetoid[convert(Variable,spec)] => stoich)
        end
    end
    sort!(rs)

    ns = Vector{Pair{Int,W2}}()
    for (spec,stoich) in maj.net_stoch
        iszero(spec) && error("Net stoichiometry can not have a species labelled 0.")
        push!(ns, statetoid[convert(Variable,spec)] => stoich)
    end
    sort!(ns)

    MassActionJump(pval, rs, ns, scale_rates = false)
end

"""
```julia
function DiffEqBase.DiscreteProblem(sys::AbstractSystem, u0map, tspan,
                                    parammap=DiffEqBase.NullParameters; kwargs...)
```

Generates a DiscreteProblem from an AbstractSystem
"""
function DiffEqBase.DiscreteProblem(sys::AbstractSystem, u0map, tspan::Tuple,
                                    parammap=DiffEqBase.NullParameters(); kwargs...)
    u0 = varmap_to_vars(u0map, states(sys))
    p = varmap_to_vars(parammap, parameters(sys))
    DiscreteProblem(u0, tspan, p; kwargs...)
end

"""
```julia
function DiffEqBase.JumpProblem(js::JumpSystem, prob, aggregator; kwargs...)
```

Generates a JumpProblem from a JumpSystem.
"""
function DiffEqJump.JumpProblem(js::JumpSystem, prob, aggregator; kwargs...)
    vrjs = Vector{VariableRateJump}()
    crjs = Vector{ConstantRateJump}()
    majs = Vector{MassActionJump}()
    pvars = parameters(js)
    statetoid = Dict(convert(Variable,state) => i for (i,state) in enumerate(states(js)))
    parammap = map(Pair,pvars,prob.p)

    for j in equations(js)
        if j isa ConstantRateJump
            push!(crjs, assemble_crj(js, j, statetoid))
        elseif j isa VariableRateJump
            push!(vrjs, assemble_vrj(js, j, statetoid))
        elseif j isa MassActionJump
            println("here")
            push!(majs, assemble_maj(js, j, statetoid, parammap))
        else
            error("JumpSystems should only contain Constant, Variable or Mass Action Jumps.")
        end
    end
    ((prob isa DiscreteProblem) && !isempty(vrjs)) && error("Use continuous problems such as an ODEProblem or a SDEProblem with VariableRateJumps")
    jset = JumpSet(Tuple(vrjs), Tuple(crjs), nothing, isempty(majs) ? nothing : majs)
    JumpProblem(prob, aggregator, jset)
end
