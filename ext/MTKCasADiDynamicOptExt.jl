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
Base.getindex(m::MXLinearInterpolation, i...) = length(i) == length(size(m.u)) ? m.u[i...] : m.u[i..., :]

mutable struct CasADiModel
    model::Opti
    U::MXLinearInterpolation
    V::MXLinearInterpolation
    tₛ::MX
    is_free_final::Bool
    solver_opti::Union{Nothing, Opti}

    function CasADiModel(opti, U, V, tₛ, is_free_final, solver_opti = nothing)
        new(opti, U, V, tₛ, is_free_final, solver_opti)
    end
end

struct CasADiDynamicOptProblem{uType, tType, isinplace, P, F, K} <:
       AbstractDynamicOptProblem{uType, tType, isinplace}
    f::F
    u0::uType
    tspan::tType
    p::P
    wrapped_model::CasADiModel
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
    colons = ntuple(_ -> (:), length(size(M.u)) - 1)
    if i < length(M.t)
        M.u[colons..., i] + Δ*(M.u[colons..., i+1] - M.u[colons..., i])
    else
        M.u[colons..., i]
    end
end

"""
    CasADiDynamicOptProblem(sys::System, u0, tspan, p; dt, steps)

Convert an System representing an optimal control system into a CasADi model
for solving using optimization. Must provide either `dt`, the timestep between collocation 
points (which, along with the timespan, determines the number of points), or directly 
provide the number of points as `steps`.

The optimization variables:
- a vector-of-vectors U representing the unknowns as an interpolation array
- a vector-of-vectors V representing the controls as an interpolation array

The constraints are:
- The set of user constraints passed to the System via `constraints`
- The solver constraints that encode the time-stepping used by the solver
"""
function MTK.CasADiDynamicOptProblem(sys::System, u0map, tspan, pmap;
        dt = nothing,
        steps = nothing,
        guesses = Dict(), kwargs...)
    MTK.process_DynamicOptProblem(CasADiDynamicOptProblem, CasADiModel, sys, u0map, tspan, pmap; dt, steps, guesses, kwargs...)
end

MTK.generate_internal_model(::Type{CasADiModel}) = CasADi.Opti()
MTK.generate_time_variable!(opti::Opti, args...) = nothing

function MTK.generate_state_variable!(model::Opti, u0, ns, tsteps)
    nt = length(tsteps)
    U = CasADi.variable!(model, ns, nt)
    set_initial!(model, U, DM(repeat(u0, 1, nt)))
    MXLinearInterpolation(U, tsteps, tsteps[2] - tsteps[1])
end

function MTK.generate_input_variable!(model::Opti, c0, nc, tsteps)
    nt = length(tsteps)
    V = CasADi.variable!(model, nc, nt)
    !isempty(c0) && set_initial!(model, V, DM(repeat(c0, 1, nt)))
    MXLinearInterpolation(V, tsteps, tsteps[2] - tsteps[1])
end

function MTK.generate_timescale!(model::Opti, guess, is_free_t)
    if is_free_t
        tₛ = variable!(model)
        set_initial!(model, tₛ, guess)
        subject_to!(model, tₛ >= 0)
        tₛ
    else
        MX(1)
    end
end

function MTK.add_constraint!(m::CasADiModel, expr)
    if expr isa Equation
        subject_to!(m.model, expr.lhs - expr.rhs == 0)
    elseif expr.relational_op === Symbolics.geq
        subject_to!(m.model, expr.lhs - expr.rhs ≥ 0)
    else
        subject_to!(m.model, expr.lhs - expr.rhs ≤ 0)
    end
end
MTK.set_objective!(m::CasADiModel, expr) = minimize!(m.model, MX(expr))

function MTK.add_initial_constraints!(m::CasADiModel, u0, u0_idxs, args...)
    @unpack model, U = m
    for i in u0_idxs
        subject_to!(model, U.u[i, 1] == u0[i])
    end
end

function MTK.substitute_model_vars(m::CasADiModel, sys, exprs, tspan)
    @unpack model, U, V, tₛ = m
    iv = MTK.get_iv(sys)
    sts = unknowns(sys)
    cts = MTK.unbound_inputs(sys)
    x_ops = [MTK.operation(MTK.unwrap(st)) for st in sts]
    c_ops = [MTK.operation(MTK.unwrap(ct)) for ct in cts]
    (ti, tf) = tspan
    if MTK.is_free_final(m)
        _tf = tₛ + ti
        exprs = map(c -> Symbolics.fast_substitute(c, Dict(tf => _tf)), exprs)
        free_t_map = Dict([[x(_tf) => U.u[i, end] for (i, x) in enumerate(x_ops)];
                           [c(_tf) => V.u[i, end] for (i, c) in enumerate(c_ops)]])
        exprs = map(c -> Symbolics.fast_substitute(c, free_t_map), exprs)
    end

    exprs = substitute_fixed_t_vars(m, sys, exprs)
    whole_interval_map = Dict([[v => U.u[i, :] for (i, v) in enumerate(sts)];
                               [v => V.u[i, :] for (i, v) in enumerate(cts)]])
    exprs = map(c -> Symbolics.fast_substitute(c, whole_interval_map), exprs)
end

function substitute_fixed_t_vars(model::CasADiModel, sys, exprs)
    stidxmap = Dict([v => i for (i, v) in enumerate(unknowns(sys))])
    ctidxmap = Dict([v => i for (i, v) in enumerate(MTK.unbound_inputs(sys))])
    iv = MTK.get_iv(sys)
    for i in 1:length(exprs)
        subvars = MTK.vars(exprs[i])
        for st in subvars 
            MTK.iscall(st) || continue
            x = operation(st)
            t = only(arguments(st))
            MTK.symbolic_type(t) === MTK.NotSymbolic() || continue
            if haskey(stidxmap, x(iv))
                idx = stidxmap[x(iv)]
                cv = model.U
            else
                idx = ctidxmap[x(iv)]
                cv = model.V
            end
            exprs[i] = Symbolics.fast_substitute(exprs[i], Dict(x(t) => cv(t)[idx]))
        end
        jcosts = Symbolics.substitute(jcosts, Dict(x(t) => cv(t)[idx]))
    end
    exprs
end

MTK.substitute_differentials(model::CasADiModel, sys, eqs) = exprs

function MTK.substitute_integral(m::CasADiModel, exprs, tspan)
    @unpack U, model, tₛ = m
    dt = U.t[2] - U.t[1]
    intmap = Dict()
    for int in MTK.collect_applied_operators(exprs, Symbolics.Integral)
        op = MTK.operation(int)
        arg = only(arguments(MTK.value(int)))
        lo, hi = MTK.value.((op.domain.domain.left, op.domain.domain.right))
        !isequal((lo, hi), tspan) &&
            error("Non-whole interval bounds for integrals are not currently supported for CasADiDynamicOptProblem.")
        # Approximate integral as sum.
        intmap[int] = dt * tₛ * sum(arg)
    end
    exprs = map(c -> Symbolics.substitute(c, intmap), exprs)
    exprs = MTK.value.(exprs)
end

function add_solve_constraints!(prob::CasADiDynamicOptProblem, tableau)
    @unpack A, α, c = tableau
    @unpack wrapped_model, f, p = prob
    @unpack model, U, V, tₛ = wrapped_model
    solver_opti = copy(model)

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

"""
CasADi Collocation solver.
- solver: an optimization solver such as Ipopt. Should be given as a string or symbol in all lowercase, e.g. "ipopt"
- tableau: An ODE RK tableau. Load a tableau by calling a function like `constructRK4` and may be found at https://docs.sciml.ai/DiffEqDevDocs/stable/internals/tableaus/. If this argument is not passed in, the solver will default to Radau second order.
"""
struct CasADiCollocation <: AbstractCollocation
    solver::Union{String, Symbol}
    tableau::DiffEqBase.ODERKTableau
end
MTK.CasADiCollocation(solver, tableau = MTK.constructDefault()) = CasADiCollocation(solver, tableau)

function MTK.prepare_and_optimize!(prob::CasADiDynamicOptProblem, solver::CasADiCollocation; verbose = false, solver_options = Dict(), plugin_options = Dict(), kwargs...)
    solver_opti = add_solve_constraints!(prob, solver.tableau)
    verbose || (solver_options["print_level"] = 0)
    solver!(solver_opti, "$(solver.solver)", plugin_options, solver_options)
    try
        CasADi.solve!(solver_opti)
    catch ErrorException
    end
    prob.wrapped_model.solver_opti = solver_opti
end

function MTK.get_U_values(model::CasADiModel)
    value_getter = MTK.successful_solve(model) ? CasADi.debug_value : CasADi.value
    (nu, nt) = size(model.U.u)
    U_vals = value_getter(model.solver_opti, model.U.u)
    size(U_vals, 2) == 1 && (U_vals = U_vals')
    U_vals = [[U_vals[i, j] for i in 1:nu] for j in 1:nt]
end

function MTK.get_V_values(model::CasADiModel)
    value_getter = MTK.successful_solve(model) ? CasADi.debug_value : CasADi.value
    (nu, nt) = size(model.V.u)
    if nu*nt != 0
        V_vals = value_getter(model.solver_opti, model.V.u)
        size(V_vals, 2) == 1 && (V_vals = V_vals')
        V_vals = [[V_vals[i, j] for i in 1:nu] for j in 1:nt]
    else
        nothing
    end
end

function MTK.get_t_values(model::CasADiModel)
    value_getter = MTK.successful_solve(model) ? CasADi.debug_value : CasADi.value
    ts = value_getter(model.solver_opti, model.tₛ) .* model.U.t
end

function MTK.successful_solve(m::CasADiModel) 
    isnothing(m.solver_opti) && return false
    retcode = CasADi.return_status(m.solver_opti)
    retcode == "Solve_Succeeded" || retcode == "Solved_To_Acceptable_Level"
end
end
