module MTKCasADiDynamicOptExt
using ModelingToolkit
using CasADi
using DiffEqBase
using UnPack
using NaNMath
const MTK = ModelingToolkit

for ff in [acos, log1p, acosh, log2, asin, tan, atanh, cos, log, sin, log10, sqrt]
    f = nameof(ff)
    @eval NaNMath.$f(x::CasadiSymbolicObject) = Base.$f(x)
end

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
    is_free_final::Bool
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
    if i < length(M.t)
        M.u[:, i] + Δ * (M.u[:, i + 1] - M.u[:, i])
    else
        M.u[:, i]
    end
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
    process_DynamicOptProblem(CasADiDynamicOptProblem, CasADiModel, sys, u0map, tspan, pmap; dt, steps, guesses, kwargs...)
end

MTK.generate_internal_model(::Type{CasADiModel}) = CasADi.opti()

function MTK.generate_state_variable(model::Opti, u0, ns, nt, tsteps)
    U = CasADi.variable!(model, ns, nt)
    set_initial!(opti, U, DM(repeat(u0, 1, steps)))
    MXLinearInterpolation(U, tsteps, tsteps[2] - tsteps[1])
end

function MTK.generate_input_variable(model::Opti, c0, nc, nt, tsteps)
    V = CasADi.variable!(model, nc, nt)
    !isempty(c0) && set_initial!(opti, V, DM(repeat(c0, 1, steps)))
    MXLinearInterpolation(V, tsteps, tsteps[2] - tsteps[1])
end

function MTK.generate_timescale(model::Opti, guess, is_free_t)
    if is_free_t
        tₛ = variable!(model)
        set_initial!(model, tₛ, guess)
        subject_to!(model, tₛ >= 0)
        tₛ
    else
        MX(1)
    end
end

function MTK.add_constraint!(model::CasADiModel, expr)
    @unpack opti = model
    if cons isa Equation
        subject_to!(opti, expr.lhs - expr.rhs == 0)
    elseif cons.relational_op === Symbolics.geq
        subject_to!(opti, expr.lhs - expr.rhs ≥ 0)
    else
        subject_to!(opti, expr.lhs - expr.rhs ≤ 0)
    end
end
MTK.set_objective!(model::CasADiModel, expr) = minimize!(model.opti, MX(expr))

function MTK.set_variable_bounds!(model, sys, pmap, tf)
    @unpack opti, U, V = model
    for (i, u) in enumerate(unknowns(sys))
        if MTK.hasbounds(u)
            lo, hi = MTK.getbounds(u)
            subject_to!(opti, Symbolics.fixpoint_sub(lo, pmap) <= U.u[i, :])
            subject_to!(opti, U.u[i, :] <= Symbolics.fixpoint_sub(hi, pmap))
        end
    end
    for (i, v) in enumerate(MTK.unbound_inputs(sys))
        if MTK.hasbounds(v)
            lo, hi = MTK.getbounds(v)
            subject_to!(opti, Symbolics.fixpoint_sub(lo, pmap) <= V.u[i, :])
            subject_to!(opti, V.u[i, :] <= Symbolics.fixpoint_sub(hi, pmap))
        end
    end
    if MTK.symbolic_type(tf) === MTK.ScalarSymbolic() && hasbounds(tf)
        lo, hi = MTK.getbounds(tf)
        subject_to!(opti, model.tₛ >= lo)
        subject_to!(opti, model.tₛ <= hi)
    end
end

function MTK.add_initial_constraints!(model::CasADiModel, u0, u0_idxs)
    @unpack opti, U = model
    for i in u0_idxs
        subject_to!(opti, U.u[i, 1] == u0[i])
    end
end

function MTK.substitute_model_vars(
        model::CasADiModel, sys, pmap, exprs; auxmap::Dict = Dict(), is_free_t)
    @unpack opti, U, V, tₛ = model
    iv = MTK.get_iv(sys)
    sts = unknowns(sys)
    cts = MTK.unbound_inputs(sys)

    x_ops = [MTK.operation(MTK.unwrap(st)) for st in sts]
    c_ops = [MTK.operation(MTK.unwrap(ct)) for ct in cts]

    exprs = map(c -> Symbolics.fast_substitute(c, auxmap), exprs)
    exprs = map(c -> Symbolics.fast_substitute(c, Dict(pmap)), exprs)
    # tf means different things in different contexts; a [tf] in a cost function
    # should be tₛ, while a x(tf) should translate to x[1]
    if is_free_t
        free_t_map = Dict([[x(tₛ) => U.u[i, end] for (i, x) in enumerate(x_ops)];
                           [c(tₛ) => V.u[i, end] for (i, c) in enumerate(c_ops)]])
        exprs = map(c -> Symbolics.fast_substitute(c, free_t_map), exprs)
    end

    exprs = substitute_fixed_t_vars(exprs)

    # for variables like x(t)
    whole_interval_map = Dict([[v => U.u[i, :] for (i, v) in enumerate(sts)];
                               [v => V.u[i, :] for (i, v) in enumerate(cts)]])
    exprs = map(c -> Symbolics.fast_substitute(c, whole_interval_map), exprs)
    exprs
end

function substitute_fixed_t_vars(exprs)
    for i in 1:length(exprs)
        subvars = MTK.vars(exprs[i])
        for st in subvars 
            MTK.iscall(st) || continue
            x = operation(st)
            t = only(arguments(st))
            MTK.symbolic_type(t) === MTK.NotSymbolic() || continue
            if haskey(stidxmap, x(iv))
                idx = stidxmap[x(iv)]
                cv = U
            else
                idx = ctidxmap[x(iv)]
                cv = V
            end
            exprs[i] = Symbolics.fast_substitute(exprs[i], Dict(x(t) => cv(t)[idx]))
        end
    end
end

MTK.substitute_differentials(model::CasADiModel, exprs, args...) = exprs

function MTK.substitute_integral(model::CasADiModel, exprs)
    @unpack U, opti = model
    dt = U.t[2] - U.t[1]
    intmap = Dict()
    for int in MTK.collect_applied_operators(exprs, Symbolics.Integral)
        op = MTK.operation(int)
        arg = only(arguments(MTK.value(int)))
        lo, hi = (op.domain.domain.left, op.domain.domain.right)
        !isequal((lo, hi), tspan) &&
            error("Non-whole interval bounds for integrals are not currently supported for CasADiDynamicOptProblem.")
        # Approximate integral as sum.
        intmap[int] = dt * tₛ * sum(arg)
    end
    exprs = map(c -> Symbolics.substitute(c, intmap), exprs)
    exprs = MTK.value.(exprs)
end

function add_solve_constraints!(prob, tableau)
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
        for k in 1:(length(tsteps) - 1)
            τ = tsteps[k]
            for (i, h) in enumerate(c)
                ΔU = sum([A[i, j] * K[j] for j in 1:(i - 1)], init = MX(zeros(nᵤ)))
                Uₙ = U.u[:, k] + ΔU * dt
                Vₙ = V.u[:, k]
                Kₙ = tₛ * f(Uₙ, Vₙ, p, τ + h * dt) # scale the time
                push!(K, Kₙ)
            end
            ΔU = dt * sum([α[i] * K[i] for i in 1:length(α)])
            subject_to!(solver_opti, U.u[:, k] + ΔU == U.u[:, k + 1])
            empty!(K)
        end
    else
        for k in 1:(length(tsteps) - 1)
            τ = tsteps[k]
            Kᵢ = variable!(solver_opti, nᵤ, length(α))
            ΔUs = A * Kᵢ'
            for (i, h) in enumerate(c)
                ΔU = ΔUs[i, :]'
                Uₙ = U.u[:, k] + ΔU * dt
                Vₙ = V.u[:, k]
                subject_to!(solver_opti, Kᵢ[:, i] == tₛ * f(Uₙ, Vₙ, p, τ + h * dt))
            end
            ΔU_tot = dt * (Kᵢ * α)
            subject_to!(solver_opti, U.u[:, k] + ΔU_tot == U.u[:, k + 1])
        end
    end
    solver_opti
end

function MTK.prepare_solver()
    opti = add_solve_constraints(prob, tableau)
    solver!(opti, "$solver", plugin_options, solver_options)
end
function MTK.get_U_values()
    U_vals = value_getter(U.u)
    size(U_vals, 2) == 1 && (U_vals = U_vals')
    U_vals = [[U_vals[i, j] for i in 1:size(U_vals, 1)] for j in 1:length(ts)]
end
function MTK.get_V_values()
end
function MTK.get_t_values()
    ts = value_getter(tₛ) * U.t
end

function MTK.optimize_model!()
    try
        sol = CasADi.solve!(opti)
        value_getter = x -> CasADi.value(sol, x)
    catch ErrorException
        value_getter = x -> CasADi.debug_value(opti, x)
        failed = true
    end
end
MTK.successful_solve() = true
end
