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

function generate_rate_function(js, rate)
    f = striplines(build_function(rate, states(js), parameters(js), independent_variable(js)))
end

function generate_affect_function(js, affect)
    # waiting on integrator build_function form
end

function assemble_vrj(js, vrj)
    rate   = generate_rate_function(js, vrj.rate)
    affect = generate_affect_function(js, vrj.affect)
    VariableRateJump(rate, affect)
end

function assemble_crj(js, crj)
    rate   = generate_rate_function(js, vrj.rate)
    affect = generate_affect_function(js, vrj.affect)
    ConstantRateJump(rate, affect)
end

function assemble_maj(maj, states_to_idxs, ps_to_idxs; scale_rate=false)

    # mass action scaled_rates need to be a Number, but 
    # operations are numbers, so can't check the type directly
    @assert !isa(maj.scaled_rates, Union{Operation,Variable})

    rstype = fieldtypes(eltype(maj.reactant_stoch))[2]
    rs     = Vector{Pair{valtype(states_to_idxs),rstype}}()
    for (spec,stoich) in maj.reactant_stoch
        push!(rs, states_to_idxs[spec] => stoich)
    end

    nstype = fieldtypes(eltype(maj.net_stoch))[2]
    ns     = Vector{Pair{valtype(states_to_idxs),nstype}}()
    for (spec,stoich) in maj.net_stoch
        push!(ns, states_to_idxs[spec] => stoich)
    end

    MassActionJump(maj.scaled_rates, rs, ns, scale_rates=scale_rate)
end