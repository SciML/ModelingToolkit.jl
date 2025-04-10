module MTKJuMPControlExt
using ModelingToolkit
using JuMP, InfiniteOpt
using DiffEqDevTools, DiffEqBase
using LinearAlgebra
const MTK = ModelingToolkit

struct JuMPControlProblem{uType, tType, isinplace, P, F, K} <:
       AbstractOptimalControlProblem{uType, tType, isinplace}
    f::F
    u0::uType
    tspan::tType
    p::P
    model::InfiniteModel
    kwargs::K

    function JuMPControlProblem(f, u0, tspan, p, model; kwargs...)
        new{typeof(u0), typeof(tspan), SciMLBase.isinplace(f),
            typeof(p), typeof(f), typeof(kwargs)}(f, u0, tspan, p, model, kwargs)
    end
end

struct InfiniteOptControlProblem{uType, tType, isinplace, P, F, K} <:
       AbstractOptimalControlProblem{uType, tType, isinplace}
    f::F
    u0::uType
    tspan::tType
    p::P
    model::InfiniteModel
    kwargs::K

    function InfiniteOptControlProblem(f, u0, tspan, p, model; kwargs...)
        new{typeof(u0), typeof(tspan), SciMLBase.isinplace(f),
            typeof(p), typeof(f), typeof(kwargs)}(f, u0, tspan, p, model, kwargs)
    end
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
function MTK.JuMPControlProblem(sys::ODESystem, u0map, tspan, pmap;
        dt = error("dt must be provided for JuMPControlProblem."),
        guesses = Dict(), kwargs...)
    MTK.warn_overdetermined(sys, u0map)
    _u0map = has_alg_eqs(sys) ? u0map : merge(Dict(u0map), Dict(guesses))
    f, u0, p = MTK.process_SciMLProblem(ODEFunction, sys, _u0map, pmap;
        t = tspan !== nothing ? tspan[1] : tspan, kwargs...)
    model = init_model(sys, tspan[1]:dt:tspan[2], u0map, u0)

    JuMPControlProblem(f, u0, tspan, p, model, kwargs...)
end

"""
    InfiniteOptControlProblem(sys::ODESystem, u0map, tspan, pmap; dt)

Convert an ODESystem representing an optimal control system into a InfiniteOpt model
for solving using optimization. Must provide `dt` for determining the length 
of the interpolation arrays.

Related to `JuMPControlProblem`, but directly adds the differential equations
of the system as derivative constraints, rather than using a solver tableau.
"""
function MTK.InfiniteOptControlProblem(sys::ODESystem, u0map, tspan, pmap;
        dt = error("dt must be provided for InfiniteOptControlProblem."),
        guesses = Dict(), kwargs...)
    MTK.warn_overdetermined(sys, u0map)
    _u0map = has_alg_eqs(sys) ? u0map : merge(Dict(u0map), Dict(guesses))
    f, u0, p = MTK.process_SciMLProblem(ODEFunction, sys, _u0map, pmap;
        t = tspan !== nothing ? tspan[1] : tspan, kwargs...)

    model = init_model(sys, tspan[1]:dt:tspan[2], u0map, u0)
    add_infopt_solve_constraints!(model, sys, pmap)
    InfiniteOptControlProblem(f, u0, tspan, p, model, kwargs...)
end

function init_model(sys, tsteps, u0map, u0)
    ctrls = controls(sys)
    states = unknowns(sys)
    model = InfiniteModel()
    @infinite_parameter(model, t in [tsteps[1], tsteps[end]], num_supports=length(tsteps))
    @variable(model, U[i = 1:length(states)], Infinite(t))
    @variable(model, V[1:length(ctrls)], Infinite(t))

    add_jump_cost_function!(model, sys)
    add_user_constraints!(model, sys)

    stidxmap = Dict([v => i for (i, v) in enumerate(states)])
    u0_idxs = has_alg_eqs(sys) ? collect(1:length(states)) :
              [stidxmap[k] for (k, v) in u0map]
    add_initial_constraints!(model, u0, u0_idxs, tsteps[1])
    return model
end

function add_jump_cost_function!(model::InfiniteModel, sys)
    jcosts = MTK.get_costs(sys)
    consolidate = MTK.get_consolidate(sys)
    if isnothing(jcosts) || isempty(jcosts)
        @objective(model, Min, 0)
        return
    end
    iv = MTK.get_iv(sys)

    stidxmap = Dict([v => i for (i, v) in enumerate(unknowns(sys))])
    cidxmap = Dict([v => i for (i, v) in enumerate(controls(sys))])

    for st in unknowns(sys)
        x = operation(st)
        t = only(arguments(st))
        idx = stidxmap[x(iv)]
        subval = isequal(t, iv) ? model[:U][idx] : model[:U][idx](t)
        jcosts = map(c -> Symbolics.substitute(c, Dict(x(t) => subval)), jcosts)
    end

    for ct in controls(sys)
        p = operation(ct)
        t = only(arguments(ct))
        idx = cidxmap[p(iv)]
        subval = isequal(t, iv) ? model[:V][idx] : model[:V][idx](t)
        jcosts = map(c -> Symbolics.substitute(c, Dict(p(t) => subval)), jcosts)
    end

    @objective(model, Min, consolidate(jcosts))
end

function add_user_constraints!(model::InfiniteModel, sys)
    conssys = MTK.get_constraintsystem(sys)
    jconstraints = isnothing(conssys) ? nothing : MTK.get_constraints(conssys)
    (isnothing(jconstraints) || isempty(jconstraints)) && return nothing

    iv = MTK.get_iv(sys)
    stidxmap = Dict([v => i for (i, v) in enumerate(unknowns(sys))])
    cidxmap = Dict([v => i for (i, v) in enumerate(controls(sys))])

    for st in unknowns(conssys)
        x = operation(st)
        t = only(arguments(st))
        idx = stidxmap[x(iv)]
        subval = isequal(t, iv) ? model[:U][idx] : model[:U][idx](t)
        jconstraints = map(c -> Symbolics.substitute(c, Dict(x(t) => subval)), jconstraints)
    end

    for ct in controls(sys)
        p = operation(ct)
        t = only(arguments(ct))
        idx = cidxmap[p(iv)]
        subval = isequal(t, iv) ? model[:V][idx] : model[:V][idx](t)
        jconstraints = map(
            c -> Symbolics.substitute(jconstraints, Dict(p(t) => subval)), jconstriants)
    end

    for (i, cons) in enumerate(jconstraints)
        if cons isa Equation
            @constraint(model, cons.lhs - cons.rhs==0, base_name="user[$i]")
        elseif cons.relational_op === Symbolics.geq
            @constraint(model, cons.lhs - cons.rhs≥0, base_name="user[$i]")
        else
            @constraint(model, cons.lhs - cons.rhs≤0, base_name="user[$i]")
        end
    end
end

function add_initial_constraints!(model::InfiniteModel, u0, u0_idxs, ts)
    U = model[:U]
    @constraint(model, initial[i in u0_idxs], U[i](ts)==u0[i])
end

is_explicit(tableau) = tableau isa DiffEqDevTools.ExplicitRKTableau

function add_infopt_solve_constraints!(model::InfiniteModel, sys, pmap)
    iv = MTK.get_iv(sys)
    t = model[:t]
    U = model[:U]
    V = model[:V]

    stmap = Dict([v => U[i] for (i, v) in enumerate(unknowns(sys))])
    ctrlmap = Dict([v => V[i] for (i, v) in enumerate(controls(sys))])
    submap = merge(stmap, ctrlmap, Dict(pmap))
    @show submap

    # Differential equations
    diff_eqs = diff_equations(sys)
    D = Differential(iv)
    diffsubmap = Dict([D(U[i]) => ∂(U[i], t) for i in 1:length(U)])
    for u in unknowns(sys)
        diff_eqs = map(e -> Symbolics.substitute(e, submap), diff_eqs)
        diff_eqs = map(e -> Symbolics.substitute(e, diffsubmap), diff_eqs)
    end
    @constraint(model, D[i = 1:length(diff_eqs)], diff_eqs[i].lhs==diff_eqs[i].rhs)

    # Algebraic equations
    alg_eqs = alg_equations(sys)
    alg_eqs = map(e -> Symbolics.substitute(e, submap), alg_eqs)
    @constraint(model, A[i = 1:length(alg_eqs)], alg_eqs[i].lhs==alg_eqs[i].rhs)
end

function add_jump_solve_constraints!(prob, tableau)
    A = tableau.A
    α = tableau.α
    c = tableau.c
    model = prob.model
    f = prob.f
    p = prob.p
    t = model[:t]
    tsteps = supports(model[:t])
    pop!(tsteps)
    dt = tsteps[2] - tsteps[1]

    U = model[:U]
    nᵤ = length(U)
    if is_explicit(tableau)
        K = Any[]
        for τ in tsteps
            for (i, h) in enumerate(c)
                ΔU = sum([A[i, j] * K[j] for j in 1:(i - 1)], init = zeros(nᵤ))
                Uₙ = [U[i](τ) + ΔU[i] * dt for i in 1:nᵤ]
                Kₙ = f(Uₙ, p, τ + h * dt)
                push!(K, Kₙ)
            end
            ΔU = dt * sum([α[i] * K[i] for i in 1:length(α)])
            @constraint(model, [n = 1:nᵤ], U[n](τ) + ΔU[n]==U[n](τ + dt),
                base_name="solve_time_$τ")
            empty!(K)
        end
    else
        @variable(model, K[1:length(α), 1:nᵤ], Infinite(t), start=tsteps[1])
        for τ in tsteps
            ΔUs = A * K
            for (i, h) in enumerate(c)
                ΔU = ΔUs[i, :]
                Uₙ = [U[j] + ΔU[j] * dt for j in 1:nᵤ]
                @constraint(model, [j in 1:nᵤ], K[i, j]==f(Uₙ, p, τ + h * dt)[j],
                    DomainRestrictions(t => τ), base_name="solve_K($τ)")
            end
            ΔU = dt * sum([α[i] * K[i, :] for i in 1:length(α)])
            @constraint(model, [n = 1:nᵤ], U[n] + ΔU[n]==U[n](τ + dt),
                DomainRestrictions(t => τ), base_name="solve_U($τ)")
        end
    end
end

"""
Solve JuMPControlProblem. Arguments:
- prob: a JumpControlProblem
- jump_solver: a LP solver such as HiGHS
- ode_solver: Takes in a symbol representing the solver. Acceptable solvers may be found at https://docs.sciml.ai/DiffEqDevDocs/stable/internals/tableaus/. Note that the symbol may be different than the typical name of the solver, e.g. :Tsitouras5 rather than Tsit5.
- silent: set the model silent (suppress model output)

Returns a JuMPControlSolution, which contains both the model and the ODE solution.
"""
function DiffEqBase.solve(
        prob::JuMPControlProblem, jump_solver, ode_solver::Symbol; silent = false)
    model = prob.model
    tableau_getter = Symbol(:construct, ode_solver)
    @eval tableau = $tableau_getter()
    silent && set_silent(model)

    # Unregister current solver constraints
    for con in all_constraints(model)
        if occursin("solve", JuMP.name(con))
            unregister(model, Symbol(JuMP.name(con)))
            delete(model, con)
        end
    end
    unregister(model, :K)
    for var in all_variables(model)
        if occursin("K", JuMP.name(var))
            delete(model, var)
        end
    end
    add_jump_solve_constraints!(prob, tableau)
    _solve(prob, jump_solver, ode_solver)
end

"""
`derivative_method` kwarg refers to the method used by InfiniteOpt to compute derivatives. The list of possible options can be found at https://infiniteopt.github.io/InfiniteOpt.jl/stable/guide/derivative/. Defaults to FiniteDifference(Backward()).
"""
function DiffEqBase.solve(prob::InfiniteOptControlProblem, jump_solver;
        derivative_method = InfiniteOpt.FiniteDifference(Backward()))
    silent && set_silent(model)
    set_derivative_method(prob.model[:t], derivative_method)
    _solve(prob, jump_solver, derivative_method)
end

function _solve(prob::AbstractOptimalControlProblem, jump_solver, solver)
    model = prob.model
    set_optimizer(model, jump_solver)
    optimize!(model)

    tstatus = termination_status(model)
    pstatus = primal_status(model)
    !has_values(model) &&
        error("Model not solvable; please report this to github.com/SciML/ModelingToolkit.jl.")

    ts = supports(model[:t])
    U_vals = value.(model[:U])
    U_vals = [[U_vals[i][j] for i in 1:length(U_vals)] for j in 1:length(ts)]
    sol = DiffEqBase.build_solution(prob, solver, ts, U_vals)

    input_sol = nothing
    if !isempty(model[:V])
        V_vals = value.(model[:V])
        V_vals = [[V_vals[i][j] for i in 1:length(V_vals)] for j in 1:length(ts)]
        input_sol = DiffEqBase.build_solution(prob, solver, ts, V_vals)
    end

    if !(pstatus === FEASIBLE_POINT &&
         (tstatus === OPTIMAL || tstatus === LOCALLY_SOLVED || tstatus === ALMOST_OPTIMAL ||
          tstatus === ALMOST_LOCALLY_SOLVED))
        sol = SciMLBase.solution_new_retcode(sol, SciMLBase.ReturnCode.ConvergenceFailure)
        !isnothing(input_sol) && (input_sol = SciMLBase.solution_new_retcode(
            input_sol, SciMLBase.ReturnCode.ConvergenceFailure))
    end

    OptimalControlSolution(model, sol, input_sol)
end

end
