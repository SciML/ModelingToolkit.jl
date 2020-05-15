JumpType = Union{VariableRateJump, ConstantRateJump, MassActionJump}

struct JumpSystem{U <: ArrayPartition} <: AbstractSystem
    eqs::U
    iv::Variable
    states::Vector{Variable}
    ps::Vector{Variable}
    name::Symbol
    systems::Vector{JumpSystem}
end

function JumpSystem(eqs, iv, states, ps; systems = JumpSystem[],
                                          name = gensym(:JumpSystem))

    ap = ArrayPartition(MassActionJump[], ConstantRateJump[], VariableRateJump[])
    for eq in eqs
        if eq isa MassActionJump 
            push!(ap.x[1], eq)
        elseif eq isa ConstantRateJump
            push!(ap.x[2], eq)
        elseif eq isa VariableRateJump
            push!(ap.x[3], eq)
        else
            error("JumpSystem equations must contain MassActionJumps, ConstantRateJumps, or VariableRateJumps.")
        end
    end

    JumpSystem{typeof(ap)}(ap, convert(Variable,iv), convert.(Variable, states), convert.(Variable, ps), name, systems)
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
        pval = simplify(substitute(sr,parammap)).value
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
    p  = varmap_to_vars(parammap, parameters(sys))
    f  = (du,u,p,t) -> du.=u    # identity function to make syms works
    df = DiscreteFunction(f, syms=Symbol.(states(sys)))
    DiscreteProblem(df, u0, tspan, p; kwargs...)
end


"""
```julia
function DiffEqBase.JumpProblem(js::JumpSystem, prob, aggregator; kwargs...)
```

Generates a JumpProblem from a JumpSystem.
"""
function DiffEqJump.JumpProblem(js::JumpSystem, prob, aggregator; kwargs...)

    statetoid = Dict(convert(Variable,state) => i for (i,state) in enumerate(states(js)))
    parammap  = map((x,y)->Pair(x(),y), parameters(js), prob.p)
    eqs       = equations(js)

    majs = MassActionJump[assemble_maj(js, j, statetoid, parammap) for j in eqs.x[1]]
    crjs = ConstantRateJump[assemble_crj(js, j, statetoid) for j in eqs.x[2]]
    vrjs = VariableRateJump[assemble_vrj(js, j, statetoid) for j in eqs.x[3]]
    ((prob isa DiscreteProblem) && !isempty(vrjs)) && error("Use continuous problems such as an ODEProblem or a SDEProblem with VariableRateJumps")
    
    jset = JumpSet(Tuple(vrjs), Tuple(crjs), nothing, isempty(majs) ? nothing : majs)

    if needs_vartojumps_map(aggregator) || needs_depgraph(aggregator)
        jdeps = asgraph(js)
        vdeps = variable_dependencies(js)
        vtoj = jdeps.badjlist
        jtov = vdeps.badjlist

        jtoj = needs_depgraph(aggregator) ? eqeq_dependencies(jdeps, vdeps).fadjlist : nothing
    else
        vtoj = nothing; jtov = nothing; jtoj = nothing
    end

    JumpProblem(prob, aggregator, jset; dep_graph=jtoj, vartojumps_map=vtoj, jumptovars_map=jtov)
end


### Functions to determine which states a jump depends on
function get_variables!(dep, jump::Union{ConstantRateJump,VariableRateJump}, variables)
    foreach(var -> (var in variables) && push!(dep, var), vars(jump.rate))
    dep
end

function get_variables!(dep, jump::MassActionJump, variables)
    jsr = jump.scaled_rates

    if jsr isa Variable
        (jsr in variables) && push!(dep, jsr)
    elseif jsr isa Operation
        foreach(var -> (var in variables) && push!(dep, var),  vars(jsr))            
    end

    for varasop in jump.reactant_stoch
        var = convert(Variable, varasop[1])
        (var in variables) && push!(dep, var)
    end

    dep
end

### Functions to determine which states are modified by a given jump
function modified_states!(mstates, jump::Union{ConstantRateJump,VariableRateJump}, sts)
    for eq in jump.affect!
        st = convert(Variable, eq.lhs)
        (st in sts) && push!(mstates, st)
    end
end

function modified_states!(mstates, jump::MassActionJump, sts)
    for (state,stoich) in jump.net_stoch
        st = convert(Variable, state)
        (st in sts) && push!(mstates, st)
    end
end