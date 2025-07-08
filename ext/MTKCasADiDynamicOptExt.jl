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
function Base.getindex(m::MXLinearInterpolation, i...)
    length(i) == length(size(m.u)) ? m.u[i...] : m.u[i..., :]
end

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
        M.u[colons..., i] + Δ*(M.u[colons..., i + 1] - M.u[colons..., i])
    else
        M.u[colons..., i]
    end
end

function MTK.CasADiDynamicOptProblem(sys::System, op, tspan;
        dt = nothing,
        steps = nothing,
        guesses = Dict(), kwargs...)
    prob,
    _ = MTK.process_DynamicOptProblem(
        CasADiDynamicOptProblem, CasADiModel, sys, op, tspan; dt, steps, guesses, kwargs...)
    prob
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

function MTK.lowered_var(m::CasADiModel, uv, i, t)
    X = getfield(m, uv)
    t isa Union{Num, Symbolics.Symbolic} ? X.u[i, :] : X(t)[i]
end

function MTK.lowered_integral(model::CasADiModel, expr, lo, hi)
    total = MX(0)
    dt = model.U.t[2] - model.U.t[1]
    for (i, t) in enumerate(model.U.t)
        if lo < t < hi
            Δt = min(dt, t - lo)
            total += (0.5*Δt*(expr[i] + expr[i - 1]))
        elseif t >= hi && (t - dt < hi)
            Δt = hi - t + dt
            total += (0.5*Δt*(expr[i] + expr[i - 1]))
        end
    end
    model.tₛ * total
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

struct CasADiCollocation <: AbstractCollocation
    solver::Union{String, Symbol}
    tableau::DiffEqBase.ODERKTableau
end

function MTK.CasADiCollocation(solver, tableau = MTK.constructDefault())
    CasADiCollocation(solver, tableau)
end

function MTK.prepare_and_optimize!(
        prob::CasADiDynamicOptProblem, solver::CasADiCollocation; verbose = false,
        solver_options = Dict(), plugin_options = Dict(), kwargs...)
    solver_opti = add_solve_constraints!(prob, solver.tableau)
    verbose || (solver_options["print_level"] = 0)
    solver!(solver_opti, "$(solver.solver)", plugin_options, solver_options)
    try
        CasADi.solve!(solver_opti)
    catch ErrorException
    end
    prob.wrapped_model.solver_opti = solver_opti
    prob.wrapped_model
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
function MTK.objective_value(model::CasADiModel)
    CasADi.pyconvert(Float64, model.solver_opti.py.value(model.solver_opti.py.f))
end

function MTK.successful_solve(m::CasADiModel)
    isnothing(m.solver_opti) && return false
    retcode = CasADi.return_status(m.solver_opti)
    retcode == "Solve_Succeeded" || retcode == "Solved_To_Acceptable_Level"
end
end
