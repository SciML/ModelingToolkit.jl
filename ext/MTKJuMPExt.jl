using JuMP, InfiniteOpt

"""
Convert an ODESystem with constraints to a JuMPProblem for optimal control solving.
"""
function JuMPProblem(sys::ODESystem, u0, tspan, p; dt = error("dt must be provided for JuMPProblem."))
    steps = tspan[1]:dt:tspan[2]
    nsteps = length(steps)
    cost = sys.cost
    coalesce = sys.coalesce

    @variables
    @variables
    @variables
end

function symbolic_to_jump_constraint(cons::Union{Num, BasicSymbolic})
    
end
