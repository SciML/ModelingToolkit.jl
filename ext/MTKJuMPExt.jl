using JuMP, InfiniteOpt

"""
Convert an ODESystem with constraints to a JuMPProblem for optimal control solving.
"""
function JuMPProblem(sys::ODESystem, u0, tspan, p; dt = error("dt must be provided for JuMPProblem."))
    ts = tspan[1]
    te = tspan[2]
    steps = ts:dt:te
    costs = get_costs(sys)
    consolidate = get_consolidate(sys)
    ctrls = get_ctrls(sys)
    states = unknowns(sys)
    constraints = get_constraints(get_constraintsystem(sys))

    model = Model()

    @infinite_parameter(model, t in [tspan[1],tspan[2]], num_supports = length(steps), derivative_method = OrthogonalCollocation(2))
    @variables(model, U[1:length(states)], Infinite(t), start = ts)
    @variables(model, V[1:length(ctrls)], Infinite(t), start = ts)
    @variables(model, K)

    jcost = generate_jump_cost_function(sys)
    @objective

    constraints = generate_jump_constraints(constraints)
    @constraints
end

function generate_jump_cost_function(costs, tsteps)
end

function generate_jump_constraints(constraints, jump_vars, jump_ps)
end

function t_to_tstep()
    
end
