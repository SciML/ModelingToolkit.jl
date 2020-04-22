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

generate_affect_function(js, affect) = build_function(affect, states(js), 
                                                      parameters(js), 
                                                      independent_variable(js),
                                                      expression=Val{false},
                                                      headerfun=add_integrator_header)[2]
function assemble_vrj(js, vrj) 
    rate   = generate_rate_function(js, vrj.rate)
    affect = generate_affect_function(js, vrj.affect!)
    VariableRateJump(rate, affect)
end

function assemble_crj(js, crj)
    rate   = generate_rate_function(js, crj.rate)
    affect = generate_affect_function(js, crj.affect!)
    ConstantRateJump(rate, affect)
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
    for j in equations(js)
        if j isa ConstantRateJump
            push!(crjs, assemble_crj(js, j))
        elseif j isa VariableRateJump
            push!(vrjs, assemble_vrj(js, j))
        else
            (j isa MassActionJump) && error("Generation of JumpProblems with MassActionJumps is not yet supported.")
        end
    end
    ((prob isa DiscreteProblem) && !isempty(vrjs)) && error("Use continuous problems such as an ODEProblem or a SDEProblem with VariableRateJumps")    
    jset = JumpSet(Tuple(vrjs), Tuple(crjs), nothing, nothing)
    JumpProblem(prob, aggregator, jset)
end