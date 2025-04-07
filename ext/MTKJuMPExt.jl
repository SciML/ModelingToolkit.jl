module MTKJuMPControlExt
using ModelingToolkit
using JuMP, InfiniteOpt
using DiffEqDevTools, DiffEqBase

struct JuMPControlProblem{uType, tType, isinplace, P, F, K} <:
       AbstractODEProblem{uType, tType, isinplace}
       f::F
       u0::uType
       tspan::tType
       p::P
       model::Model
       kwargs::K
end

"""
    JuMPControlProblem(sys::ODESystem, u0, tspan, p; dt)

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
function JuMPControlProblem(sys::ODESystem, u0map, tspan, pmap; dt = error("dt must be provided for JuMPControlProblem."), guesses, eval_expression, eval_module)
    ts = tspan[1]
    te = tspan[2]
    steps = ts:dt:te
    ctrls = get_ctrls(sys)
    states = unknowns(sys)

    if !isnothing(constraintsys)
        (length(constraints(constraintsys)) + length(u0map) > length(sts)) &&
            @warn "The JuMPControlProblem is overdetermined. The total number of conditions (# constraints + # fixed initial values given by u0map) exceeds the total number of states. The solvers will default to doing a nonlinear least-squares optimization."
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

    stidxmap = Dict([v => i for (i, v) in enumerate(sts)])
    u0_idxs = has_alg_eqs(sys) ? collect(1:length(sts)) : [stidxmap[k] for (k, v) in u0map]
    add_initial_constraints!(model, u0, u0_idxs, tspan)

    JuMPControlProblem{iip}(f, u0, tspan, p, model, kwargs...)
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
end

function add_initial_constraints!(model, u0, u0_idxs, tspan)
    ts = tspan[1]
    @constraint(model, init_u0_idx[i in u0_idxs], U[i](ts) == u0[i])
end

is_explicit(tableau) = tableau isa DiffEqDevTools.ExplicitRKTableau

function add_solve_constraints!(prob, tableau, f, tsteps)
    A = tableau.A
    α = tableau.α
    c = tableau.c
    model = prob.model
    p = prob.p
    dt = step(tsteps)

    U = model[:U]
    if is_explicit(tableau)
        K = Any[]
        for τ in tsteps
            for (i, h) in enumerate(c)
                ΔU = sum([A[i, j] * K[j] for j in 1:i-1])
                Kₙ = f(U + ΔU*dt, p, τ + h*dt)
                push!(K, Kₙ)
            end
            @constraint(model, U(τ) + dot(α, K) == U(τ + dt))
            empty!(K)
        end
    else
        @variable(model, K[1:length(a)], Infinite(t), start = tsteps[1])
        for τ in tsteps
            ΔUs = A * K(τ)
            for (i, h) in enumerate(c)
                ΔU = ΔUs[i]
                @constraint(model, K[i](τ) == f(U(τ) + ΔU*dt, p, τ + h*dt))
            end
            @constraint(model, U(τ) + dot(α, K(τ)) == U(τ + dt))
        end
    end
end

"""
"""
struct JuMPControlSolution
    model
    sol::ODESolution
end

"""
Solve JuMPProblem. Takes in a symbol representing the solver. Acceptable solvers may be found at https://docs.sciml.ai/DiffEqDevDocs/stable/internals/tableaus/.
Note that the symbol may be different than the typical
name of the solver, e.g. :Tsitouras5 rather than Tsit5.
"""
function solve(prob::JuMPProblem, jump_solver, ode_solver::Symbol)
    model = prob.model
    f = prob.f
    tableau_getter = Symbol(:construct, ode_solver)
    @eval tableau = $tableau_getter()
    ts = prob.tspan[1]:dt:prob.tspan[2]
    add_solve_constraints!(model, ts, tableau, f)

    set_optimizer(model, jump_solver)
    optimize!(model)

    if is_solved_and_feasible(model)
        sol = DiffEqBase.build_solution(prob, ode_solver, ts, value.(U))
        JuMPControlSolution(model, sol)
    end
end

end
