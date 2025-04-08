module MTKJuMPControlExt
using ModelingToolkit
using JuMP, InfiniteOpt
using DiffEqDevTools, DiffEqBase, SciMLBase
using LinearAlgebra
const MTK = ModelingToolkit

struct JuMPControlProblem{uType, tType, isinplace, P, F, K} <: SciMLBase.AbstractODEProblem{uType, tType, isinplace}
       f::F
       u0::uType
       tspan::tType
       p::P
       model::InfiniteModel
       kwargs::K

       function JuMPControlProblem(f, u0, tspan, p, model; kwargs...) 
           new{typeof(u0), typeof(tspan), SciMLBase.isinplace(f), typeof(p), typeof(f), typeof(kwargs)}(f, u0, tspan, p, model, kwargs)
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
function MTK.JuMPControlProblem(sys::ODESystem, u0map, tspan, pmap; dt = error("dt must be provided for JuMPControlProblem."), kwargs...)
    ts = tspan[1]
    te = tspan[2]
    steps = ts:dt:te
    ctrls = controls(sys)
    states = unknowns(sys)
    constraintsys = MTK.get_constraintsystem(sys)

    if !isnothing(constraintsys)
        (length(constraints(constraintsys)) + length(u0map) > length(sts)) &&
            @warn "The JuMPControlProblem is overdetermined. The total number of conditions (# constraints + # fixed initial values given by u0map) exceeds the total number of states. The solvers will default to doing a nonlinear least-squares optimization."
    end

    f, u0, p = MTK.process_SciMLProblem(ODEFunction, sys, u0map, pmap;
        t = tspan !== nothing ? tspan[1] : tspan, kwargs...)

    model = InfiniteModel()
    @infinite_parameter(model, t in [ts, te], num_supports = length(steps))
    @variable(model, U[i = 1:length(states)], Infinite(t))
    @variable(model, V[1:length(ctrls)], Infinite(t))

    add_jump_cost_function!(model, sys)
    add_user_constraints!(model, sys)

    stidxmap = Dict([v => i for (i, v) in enumerate(states)])
    u0_idxs = has_alg_eqs(sys) ? collect(1:length(states)) : [stidxmap[k] for (k, v) in u0map]
    add_initial_constraints!(model, u0, u0_idxs, tspan)

    JuMPControlProblem(f, u0, tspan, p, model, kwargs...)
end

function add_jump_cost_function!(model, sys)
    jcosts = MTK.get_costs(sys)
    consolidate = MTK.get_consolidate(sys)
    if isnothing(consolidate)
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
        jcosts = Symbolics.substitute(jcosts, Dict(x(t) => subval))
    end
    
    for ct in controls(sys)
        p = operation(ct)
        t = only(arguments(ct))
        idx = cidxmap[p(iv)]
        subval = isequal(t, iv) ? model[:V][idx] : model[:V][idx](t)
        jcosts = Symbolics.substitute(jcosts, Dict(x(t) => subval))
    end
    
    @objective(model, Min, consolidate(jcosts))
end

function add_user_constraints!(model, sys)
    jconstraints = if !(csys = MTK.get_constraintsystem(sys) isa Nothing)
        MTK.get_constraints(csys)
    else
        nothing
    end
    isnothing(jconstraints) && return nothing

    iv = MTK.get_iv(sys)
    stidxmap = Dict([v => i for (i, v) in enumerate(unknowns(sys))])
    cidxmap = Dict([v => i for (i, v) in enumerate(controls(sys))])

    for st in unknowns(sys)
        x = operation(st)
        t = only(arguments(st))
        idx = stidxmap[x(iv)]
        subval = isequal(t, iv) ? model[:U][idx] : model[:U][idx](t)
        jconstraints = Symbolics.substitute(jconstraints, Dict(x(t) => subval))
    end

    for ct in controls(sys)
        p = operation(ct)
        t = only(arguments(ct))
        idx = cidxmap[p(iv)]
        subval = isequal(t, iv) ? model[:V][idx] : model[:V][idx](t)
        jconstraints = Symbolics.substitute(jconstraints, Dict(p(t) => subval))
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
    U = model[:U]
    @constraint(model, initial[i in u0_idxs], U[i](ts) == u0[i])
end

is_explicit(tableau) = tableau isa DiffEqDevTools.ExplicitRKTableau

function add_solve_constraints!(prob, tableau)
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
                ΔU = sum([A[i, j] * K[j] for j in 1:i-1], init = zeros(nᵤ))
                Uₙ = [U[i](τ) + ΔU[i]*dt for i in 1:nᵤ]
                Kₙ = f(Uₙ, p, τ + h*dt)
                push!(K, Kₙ)
            end
            ΔU = dt*sum([α[i] * K[i] for i in 1:length(α)])
            @constraint(model, [n = 1:nᵤ], U[n](τ) + ΔU[n] == U[n](τ + dt), base_name = "solve_time_$τ")
            empty!(K)
        end
    else
        @variable(model, K[1:length(α), 1:nᵤ], Infinite(t), start = tsteps[1])
        for τ in tsteps
            ΔUs = A * K
            for (i, h) in enumerate(c)
                ΔU = ΔUs[i, :]
                Uₙ = [U[j] + ΔU[j]*dt for j in 1:nᵤ]
                @constraint(model, [j in 1:nᵤ], K[i, j] == f(Uₙ, p, τ + h*dt)[j], DomainRestrictions(t => τ), base_name = "solve_K($τ)")
            end
            ΔU = dt*sum([α[i] * K[i, :] for i in 1:length(α)])
            @constraint(model, [n = 1:nᵤ], U[n] + ΔU[n] == U[n](τ + dt), DomainRestrictions(t => τ), base_name = "solve_U($τ)")
        end
    end
end

"""
"""
struct JuMPControlSolution
    model::InfiniteModel
    sol::ODESolution
end

"""
Solve JuMPControlProblem. Arguments:
- prob: a JumpControlProblem
- jump_solver: a LP solver such as HiGHS
- ode_solver: Takes in a symbol representing the solver. Acceptable solvers may be found at https://docs.sciml.ai/DiffEqDevDocs/stable/internals/tableaus/. Note that the symbol may be different than the typical name of the solver, e.g. :Tsitouras5 rather than Tsit5.

Returns a JuMPControlSolution, which contains both the model and the ODE solution.
"""
function DiffEqBase.solve(prob::JuMPControlProblem, jump_solver, ode_solver::Symbol)
    model = prob.model
    tableau_getter = Symbol(:construct, ode_solver)
    @eval tableau = $tableau_getter()
    ts = supports(model[:t])

    # Unregister current solver constraints
    for con in all_constraints(model)
        if occursin("solve", JuMP.name(con))
            unregister(model, Symbol(JuMP.name(con)))
            delete(model, con)
        end
    end
    for var in all_variables(model)
        @show JuMP.name(var)
        if occursin("K", JuMP.name(var))
            unregister(model, Symbol(JuMP.name(var)))
            delete(model, var)
        end
    end
    add_solve_constraints!(prob, tableau)

    set_optimizer(model, jump_solver)
    optimize!(model)

    tstatus = termination_status(model)
    pstatus = primal_status(model)
    !has_values(model) && error("Model not solvable; please report this to github.com/SciML/ModelingToolkit.jl.")

    U_vals = value.(model[:U])
    U_vals = [[U_vals[i][j] for i in 1:length(U_vals)] for j in 1:length(ts)]
    sol = DiffEqBase.build_solution(prob, ode_solver, ts, U_vals)

    if !(pstatus === FEASIBLE_POINT && (tstatus === OPTIMAL || tstatus === LOCALLY_SOLVED || tstatus === ALMOST_OPTIMAL || tstatus === ALMOST_LOCALLY_SOLVED))
        sol = SciMLBase.solution_new_retcode(sol, SciMLBase.ReturnCode.ConvergenceFailure)
    end
    JuMPControlSolution(model, sol)
end

end
