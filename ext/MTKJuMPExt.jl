module MTKJuMPControlExt
using ModelingToolkit
using JuMP, InfiniteOpt
using DiffEqDevTools, DiffEqBase

struct JuMPProblem{uType, tType, isinplace, P, F, K} <:
       AbstractODEProblem{uType, tType, isinplace}
       f::F
       u0::uType
       tspan
       p
       model
       kwargs
end

"""
    JuMPProblem(sys::ODESystem, u0, tspan, p; dt)

Convert an ODESystem representing an optimal control system into a JuMP model
for solving using optimization. Must provide `dt` for determining the length 
of the interpolation arrays.

The optimization variables:
- a vector-of-vectors U representing the unknowns as an interpolation array
- a vector-of-vectors V representing the controls as an interpolation array

The constraints are:
- The set of user constraints passed to the ODESystem via `constraints`
- The solver constraints that encode the time-stepping used by the solver
"""
function JuMPProblem(sys::ODESystem, u0map, tspan, pmap; dt = error("dt must be provided for JuMPProblem."), solver = :Tsit5)
    ts = tspan[1]
    te = tspan[2]
    steps = ts:dt:te
    ctrls = get_ctrls(sys)
    states = unknowns(sys)

    if !isnothing(constraintsys)
        (length(constraints(constraintsys)) + length(u0map) > length(sts)) &&
            @warn "The BVProblem is overdetermined. The total number of conditions (# constraints + # fixed initial values given by u0map) exceeds the total number of states. The BVP solvers will default to doing a nonlinear least-squares optimization."
    end

    model = InfiniteModel()
    @infinite_parameter(model, t in [ts, te], num_supports = length(steps), derivative_method = OrthogonalCollocation(2))
    @variable(model, U[1:length(states)], Infinite(t), start = ts)
    @variable(model, V[1:length(ctrls)], Infinite(t), start = ts)
    @variable(model, K)

    f, u0, p = process_SciMLProblem(ODEFunction{iip, specialize}, sys, _u0map, parammap;
        t = tspan !== nothing ? tspan[1] : tspan, guesses,
        check_length, warn_initialize_determined, eval_expression, eval_module, kwargs...)

    add_jump_cost_function!(model, sys)
    add_user_constraints!(model, sys)
    add_solve_constraints!(model)

    JuMPProblem{iip}(f, u0, tspan, p, model; kwargs...)
end

function add_jump_cost_function!(model, sys)
    jcosts = get_costs(sys)
    consolidate = get_consolidate(sys)
    iv = get_iv(sys)

    stidxmap = Dict([v => i for (i, v) in enumerate(get_unknowns(sys))])
    cidxmap = Dict([v => i for (i, v) in enumerate(get_ctrls(sys))])

    for st in get_unknowns(sys)
        x = operation(st)
        t = only(arguments(st))
        idx = stidxmap[x(iv)]
        jcosts = Symbolics.substitute(costs, Dict(x(t) => model[:U][idx](t)))
    end
    
    for ct in get_ctrls(sys)
        p = operation(ct)
        t = only(arguments(ct))
        idx = cidxmap[p(iv)]
        jcosts = Symbolics.substitute(costs, Dict(p(t) => model[:V][idx](t)))
    end
    
    @objective(model, Min, consolidate(jcosts))
end

function add_user_constraints!(model, sys, u0map)
    jconstraints = get_constraints(get_constraintsystem(sys))
    iv = get_iv(sys)

    stidxmap = Dict([v => i for (i, v) in enumerate(get_unknowns(sys))])
    cidxmap = Dict([v => i for (i, v) in enumerate(get_ctrls(sys))])

    for st in get_unknowns(sys)
        x = operation(st)
        t = only(arguments(st))
        idx = stidxmap[x(iv)]
        subval = isequal(t, iv) ? model[:U][idx] : model[:U][idx](t)
        jconstraints = Symbolics.substitute(constraints, Dict(x(t) => subval))
    end

    for ct in get_ctrls(sys)
        p = operation(ct)
        t = only(arguments(ct))
        idx = cidxmap[p(iv)]
        subval = isequal(t, iv) ? model[:V][idx] : model[:V][idx](t)
        jconstraints = Symbolics.substitute(constraints, Dict(p(t) => subval))
    end

    for (i, cons) in enumerate(jconstraints)
        if cons isa Equation
            @constraint(model, user[i], cons.lhs - cons.rhs == 0)
        elseif cons.relational_op === Symbolics.geq 
            @constraint(model, user[i], cons.lhs - cons.rhs ≥ 0)
        else
            @constraint(model, user[i], cons.lhs - cons.rhs ≤ 0)
        end
    end

    # Add initial constraints.
end

function add_solve_constraints!(model, tsteps, solver)
    tableau = fetch_tableau(solver)

    for (i, t) in collect(enumerate(tsteps))
    end
end

"""
Solve JuMPProblem. Takes in a symbol representing the solver.
"""
function solve(prob::JuMPProblem, solver_sym::Symbol)
    model = prob.model
    tableau_getter = Symbol(:construct, solver)
    @eval tableau = $tableau_getter()
    add_solve_constraints!(model, tableau)
end
end
