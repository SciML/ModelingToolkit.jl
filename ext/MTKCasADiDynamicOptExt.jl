module MTKCasADiDynamicOptExt
using ModelingToolkit
using CasADi
using DiffEqBase
using UnPack
const MTK = ModelingToolkit

# Default linear interpolation for MX objects, likely to change down the line when we support interpolation with the collocation polynomial.
struct MXLinearInterpolation
    u::MX
    t::Vector{Float64}
    dt::Float64
end

struct CasADiModel
    opti::Opti
    U::MXLinearInterpolation
    V::MXLinearInterpolation
    tₛ::MX
end

struct CasADiDynamicOptProblem{uType, tType, isinplace, P, F, K} <:
       AbstractDynamicOptProblem{uType, tType, isinplace}
    f::F
    u0::uType
    tspan::tType
    p::P
    model::CasADiModel
    kwargs::K

    function CasADiDynamicOptProblem(f, u0, tspan, p, model, kwargs...)
        new{typeof(u0), typeof(tspan), SciMLBase.isinplace(f, 5),
            typeof(p), typeof(f), typeof(kwargs)}(f, u0, tspan, p, model, kwargs)
    end
end

function (M::MXLinearInterpolation)(τ) 
    nt = (τ - M.t[1]) / M.dt
    i = 1 + floor(Int, nt)
    Δ = nt - i + 1

    (i > length(M.t) || i < 1) && error("Cannot extrapolate past the tspan.")
    M.u[:, i] + Δ*(M.u[:, i + 1] - M.u[:, i])
end

"""
    CasADiDynamicOptProblem(sys::ODESystem, u0, tspan, p; dt, steps)

Convert an ODESystem representing an optimal control system into a CasADi model
for solving using optimization. Must provide either `dt`, the timestep between collocation 
points (which, along with the timespan, determines the number of points), or directly 
provide the number of points as `steps`.

The optimization variables:
- a vector-of-vectors U representing the unknowns as an interpolation array
- a vector-of-vectors V representing the controls as an interpolation array

The constraints are:
- The set of user constraints passed to the ODESystem via `constraints`
- The solver constraints that encode the time-stepping used by the solver
"""
function MTK.CasADiDynamicOptProblem(sys::ODESystem, u0map, tspan, pmap;
        dt = nothing,
        steps = nothing,
        guesses = Dict(), kwargs...) 
    MTK.warn_overdetermined(sys, u0map)
    _u0map = has_alg_eqs(sys) ? u0map : merge(Dict(u0map), Dict(guesses))
    f, u0, p = MTK.process_SciMLProblem(ODEInputFunction, sys, _u0map, pmap;
        t = tspan !== nothing ? tspan[1] : tspan, output_type = MX, kwargs...)

    pmap = Dict{Any, Any}(pmap)
    steps, is_free_t = MTK.process_tspan(tspan, dt, steps)
    model = init_model(sys, tspan, steps, u0map, pmap, u0; is_free_t)

    CasADiDynamicOptProblem(f, u0, tspan, p, model, kwargs...)
end

function init_model(sys, tspan, steps, u0map, pmap, u0; is_free_t = false)
    ctrls = MTK.unbound_inputs(sys)
    states = unknowns(sys)
    opti = CasADi.Opti()

    if is_free_t
        (ts_sym, te_sym) = tspan
        MTK.symbolic_type(ts_sym) !== MTK.NotSymbolic() &&
            error("Free initial time problems are not currently supported.")
        tₛ = variable!(opti)
        tsteps = LinRange(0, 1, steps)
    else
        tₛ = MX(1)
        tsteps = LinRange(tspan[1], tspan[2], steps)
    end
    
    U = CasADi.variable!(opti, length(states), steps)
    V = CasADi.variable!(opti, length(ctrls), steps)
    U_interp = MXLinearInterpolation(U, tsteps, tsteps[2]-tsteps[1])
    V_interp = MXLinearInterpolation(V, tsteps, tsteps[2]-tsteps[1])

    model = CasADiModel(opti, U_interp, V_interp, tₛ)

    set_casadi_bounds!(model, sys, pmap)
    add_cost_function!(model, sys, (tspan[1], tspan[2]), pmap)
    add_user_constraints!(model, sys, pmap; is_free_t)

    stidxmap = Dict([v => i for (i, v) in enumerate(states)])
    u0map = Dict([MTK.default_toterm(MTK.value(k)) => v for (k, v) in u0map])
    u0_idxs = has_alg_eqs(sys) ? collect(1:length(states)) :
              [stidxmap[MTK.default_toterm(k)] for (k, v) in u0map]
    add_initial_constraints!(model, u0, u0_idxs)

    model
end

function set_casadi_bounds!(model, sys, pmap)
    @unpack opti, U, V = model
    for (i, u) in enumerate(unknowns(sys))
        if MTK.hasbounds(u)
            lo, hi = MTK.getbounds(u)
            subject_to!(opti, lo <= U[i, :] <= hi)
        end
    end
    for (i, v) in enumerate(MTK.unbound_inputs(sys))
        if MTK.hasbounds(v)
            lo, hi = MTK.getbounds(v)
            subject_to!(opti, lo <= V[i, :] <= hi)
        end
    end
end

function add_initial_constraints!(model::CasADiModel, u0, u0_idxs)
    @unpack opti, U = model
    for i in u0_idxs
        subject_to!(opti, U.u[i, 1] == u0[i])
    end
end

function add_user_constraints!(model::CasADiModel, sys, pmap; is_free_t = false)
    @unpack opti, U, V, tₛ = model

    iv = MTK.get_iv(sys)
    conssys = MTK.get_constraintsystem(sys)
    jconstraints = isnothing(conssys) ? nothing : MTK.get_constraints(conssys)
    (isnothing(jconstraints) || isempty(jconstraints)) && return nothing

    stidxmap = Dict([v => i for (i, v) in enumerate(unknowns(sys))])
    cons_unknowns = map(MTK.default_toterm, unknowns(conssys))
    for st in cons_unknowns
        x = MTK.operation(st)
        t = only(MTK.arguments(st))
        idx = stidxmap[x(iv)]
        @show t
        MTK.symbolic_type(t) === MTK.NotSymbolic() || continue
        jconstraints = map(c -> Symbolics.substitute(c, Dict(x(t) => U(t)[idx])), jconstraints)
    end
    jconstraints = substitute_casadi_vars(model, sys, pmap, jconstraints)

    for (i, cons) in enumerate(jconstraints)
        if cons isa Equation
            subject_to!(opti, cons.lhs - cons.rhs==0)
        elseif cons.relational_op === Symbolics.geq
            subject_to!(opti, cons.lhs - cons.rhs≥0)
        else
            subject_to!(opti, cons.lhs - cons.rhs≤0)
        end
    end
end

function add_cost_function!(model::CasADiModel, sys, tspan, pmap)
    @unpack opti, U, V, tₛ = model
    jcosts = MTK.get_costs(sys)
    consolidate = MTK.get_consolidate(sys)

    if isnothing(jcosts) || isempty(jcosts)
        minimize!(opti, MX(0))
        return
    end
    stidxmap = Dict([v => i for (i, v) in enumerate(sts)])
    pidxmap = Dict([v => i for (i, v) in enumerate(ps)])

    for i in 1:length(jcosts)
        vars = vars(jcosts[i])
        for st in vars
            x = operation(st)
            t = only(arguments(st))
            t isa Union{Num, MTK.Symbolic} && continue
            idx = stidxmap[x(iv)]
            jcosts[i] = Symbolics.substitute(jcosts[i], Dict(x(t) => U(t)[idx]))
        end
    end
    jcosts = substitute_casadi_vars(model::CasADiModel, sys, pmap, jcosts; auxmap)

    dt = U.t[2] - U.t[1]
    intmap = Dict()
    for int in MTK.collect_applied_operators(jcosts, Symbolics.Integral)
        op = MTK.operation(int)
        arg = only(arguments(MTK.value(int)))
        lo, hi = (op.domain.domain.left, op.domain.domain.right)
        (lo, hi) !== tspan && error("Non-whole interval bounds for integrals are not currently supported.")
        intmap[int] = dt * tₛ * sum(arg)
    end
    jcosts = map(c -> Symbolics.substitute(c, intmap), jcosts)
    minimize!(opti, consolidate(jcosts))
end

function substitute_casadi_vars(model::CasADiModel, sys, pmap, exprs; auxmap = Dict())
    @unpack opti, U, V = model
    iv = MTK.get_iv(sys)
    sts = unknowns(sys)
    cts = MTK.unbound_inputs(sys)

    x_ops = [MTK.operation(MTK.unwrap(st)) for st in sts]
    c_ops = [MTK.operation(MTK.unwrap(ct)) for ct in cts]

    exprs = map(c -> Symbolics.fixpoint_sub(c, auxmap), exprs)
    exprs = map(c -> Symbolics.fixpoint_sub(c, Dict(pmap)), exprs)

    # for variables like x(t)
    whole_interval_map = Dict([[v => U.u[i, :] for (i, v) in enumerate(sts)];
                               [v => V.u[i, :] for (i, v) in enumerate(cts)]])
    exprs = map(c -> Symbolics.fixpoint_sub(c, whole_interval_map), exprs)
    exprs
end

function add_solve_constraints(prob, tableau; is_free_t = false)
    @unpack A, α, c = tableau
    @unpack model, f, p = prob 
    @unpack opti, U, V, tₛ = model
    solver_opti = copy(opti)

    tsteps = U.t 
    dt = tsteps[2] - tsteps[1]

    nᵤ = size(U.u, 1)
    nᵥ = size(V.u, 1)

    if MTK.is_explicit(tableau)
        K = MX[]
        for k in 1:length(tsteps)-1
            τ = tsteps[k]
            for (i, h) in enumerate(c)
                ΔU = sum([A[i, j] * K[j] for j in 1:(i - 1)], init = MX(zeros(nᵤ)))
                Uₙ = U.u[:, k] + ΔU*dt
                Vₙ = V.u[:, k]
                Kₙ = tₛ * f(Uₙ, Vₙ, p, τ + h * dt) # scale the time
                push!(K, Kₙ)
            end
            ΔU = dt * sum([α[i] * K[i] for i in 1:length(α)])
            subject_to!(solver_opti, U.u[:, k] + ΔU == U.u[:, k+1])
            empty!(K)
        end
    else
        for k in 1:length(tsteps)-1
            τ = tsteps[k]
            Kᵢ = variable!(solver_opti, nᵤ, length(α))
            ΔUs = A * Kᵢ' # the stepsize at each stage of the implicit method
            for (i, h) in enumerate(c)
                ΔU = ΔUs[i,:]'
                Uₙ = U.u[:,k] + ΔU*dt
                Vₙ = V.u[:,k]
                subject_to!(solver_opti, Kᵢ[:,i] == tₛ * f(Uₙ, Vₙ, p, τ + h*dt))
            end
            ΔU_tot = dt*(Kᵢ*α)
            subject_to!(solver_opti, U.u[:, k] + ΔU_tot == U.u[:,k+1])
        end
    end
    solver_opti
end

"""
    solve(prob::CasADiDynamicOptProblem, casadi_solver, ode_solver; plugin_options, solver_options, silent)

`plugin_options` and `solver_options` get propagated to the Opti object in CasADi.

NOTE: the solver should be passed in as a string to CasADi. "ipopt"
"""
function DiffEqBase.solve(prob::CasADiDynamicOptProblem, solver::Union{String, Symbol} = "ipopt", tableau_getter = MTK.constructDefault; plugin_options::Dict = Dict(), solver_options::Dict = Dict(), silent = false)
    @unpack model, u0, p, tspan, f = prob
    tableau = tableau_getter()
    @unpack opti, U, V, tₛ = model

    opti = add_solve_constraints(prob, tableau)
    solver!(opti, "$solver", plugin_options, solver_options)

    failed = false
    value_getter = nothing
    try
        sol = CasADi.solve!(opti)
        value_getter = x -> CasADi.value(sol, x)
    catch ErrorException
        value_getter = x -> CasADi.debug_value(opti, x)
        failed = true
    end

    ts = value_getter(tₛ) * U.t
    U_vals = value_getter(U.u)
    U_vals = [[U_vals[i, j] for i in 1:size(U_vals, 1)] for j in 1:length(ts)]
    sol = DiffEqBase.build_solution(prob, tableau_getter, ts, U_vals)

    input_sol = nothing
    if prod(size(V.u)) != 0
        V_vals = value_getter(V.u)
        V_vals = [[V_vals[i, j] for i in 1:size(V_vals, 1)] for j in 1:length(ts)]
        input_sol = DiffEqBase.build_solution(prob, tableau_getter, ts, V_vals)
    end

    if failed
        sol = SciMLBase.solution_new_retcode(sol, SciMLBase.ReturnCode.ConvergenceFailure)
        !isnothing(input_sol) && (input_sol = SciMLBase.solution_new_retcode(
            input_sol, SciMLBase.ReturnCode.ConvergenceFailure))
    end

    DynamicOptSolution(model, sol, input_sol)
end
end
