module MTKInfiniteOptExt
using ModelingToolkit
using InfiniteOpt
using DiffEqBase
using LinearAlgebra
using StaticArrays
using UnPack
import SymbolicUtils
import NaNMath
const MTK = ModelingToolkit

struct InfiniteOptModel
    model::InfiniteModel
    U::Vector{<:AbstractVariableRef}
    V::Vector{<:AbstractVariableRef}
    tₛ::AbstractVariableRef
    is_free_final::Bool
end

struct JuMPDynamicOptProblem{uType, tType, isinplace, P, F, K} <:
       AbstractDynamicOptProblem{uType, tType, isinplace}
    f::F
    u0::uType
    tspan::tType
    p::P
    wrapped_model::InfiniteOptModel
    kwargs::K

    function JuMPDynamicOptProblem(f, u0, tspan, p, model, kwargs...)
        new{typeof(u0), typeof(tspan), SciMLBase.isinplace(f, 5),
            typeof(p), typeof(f), typeof(kwargs)}(f, u0, tspan, p, model, kwargs)
    end
end

struct InfiniteOptDynamicOptProblem{uType, tType, isinplace, P, F, K} <:
       AbstractDynamicOptProblem{uType, tType, isinplace}
    f::F
    u0::uType
    tspan::tType
    p::P
    wrapped_model::InfiniteOptModel
    kwargs::K

    function InfiniteOptDynamicOptProblem(f, u0, tspan, p, model, kwargs...)
        new{typeof(u0), typeof(tspan), SciMLBase.isinplace(f),
            typeof(p), typeof(f), typeof(kwargs)}(f, u0, tspan, p, model, kwargs)
    end
end

MTK.generate_internal_model(m::Type{InfiniteOptModel}) = InfiniteModel()
function MTK.generate_time_variable!(m::InfiniteModel, tspan, tsteps)
    @infinite_parameter(m, t in [tspan[1], tspan[2]], num_supports = length(tsteps))
end
function MTK.generate_state_variable!(m::InfiniteModel, u0::Vector, ns, ts)
    @variable(m, U[i = 1:ns], Infinite(m[:t]), start=u0[i])
end
function MTK.generate_input_variable!(m::InfiniteModel, c0, nc, ts)
    @variable(m, V[i = 1:nc], Infinite(m[:t]), start=c0[i])
end

function MTK.generate_timescale!(m::InfiniteModel, guess, is_free_t)
    @variable(m, tₛ ≥ 0, start = guess)
    if !is_free_t
        fix(tₛ, 1, force = true)
        set_start_value(tₛ, 1)
    end
    tₛ
end

function MTK.add_constraint!(m::InfiniteOptModel, expr::Union{Equation, Inequality})
    if expr isa Equation
        @constraint(m.model, expr.lhs - expr.rhs == 0)
    elseif expr.relational_op === Symbolics.geq
        @constraint(m.model, expr.lhs - expr.rhs ≥ 0)
    else
        @constraint(m.model, expr.lhs - expr.rhs ≤ 0)
    end
end
MTK.set_objective!(m::InfiniteOptModel, expr) = @objective(m.model, Min, expr)

function MTK.JuMPDynamicOptProblem(sys::System, op, tspan;
        dt = nothing,
        steps = nothing,
        guesses = Dict(), kwargs...)
    prob,
    _ = MTK.process_DynamicOptProblem(JuMPDynamicOptProblem, InfiniteOptModel, sys,
        op, tspan; dt, steps, guesses, kwargs...)
    prob
end

function MTK.InfiniteOptDynamicOptProblem(sys::System, op, tspan;
        dt = nothing,
        steps = nothing,
        guesses = Dict(), kwargs...)
    prob,
    pmap = MTK.process_DynamicOptProblem(InfiniteOptDynamicOptProblem, InfiniteOptModel,
        sys, op, tspan; dt, steps, guesses, kwargs...)
    MTK.add_equational_constraints!(prob.wrapped_model, sys, pmap, tspan)
    prob
end

function MTK.lowered_integral(model::InfiniteOptModel, expr, lo, hi)
    model.tₛ * InfiniteOpt.∫(expr, model.model[:t], lo, hi)
end
MTK.lowered_derivative(model::InfiniteOptModel, i) = ∂(model.U[i], model.model[:t])

function MTK.process_integral_bounds(model::InfiniteOptModel, integral_span, tspan)
    if MTK.is_free_final(model) && isequal(integral_span, tspan)
        integral_span = (0, 1)
    elseif MTK.is_free_final(model)
        error("Free final time problems cannot handle partial timespans.")
    else
        integral_span
    end
end

function MTK.add_initial_constraints!(m::InfiniteOptModel, u0, u0_idxs, ts)
    for i in u0_idxs
        fix(m.U[i](0), u0[i], force = true)
    end
end

function MTK.lowered_var(m::InfiniteOptModel, uv, i, t)
    X = getfield(m, uv)
    t isa Union{Num, Symbolics.Symbolic} ? X[i] : X[i](t)
end

function add_solve_constraints!(prob::JuMPDynamicOptProblem, tableau)
    @unpack A, α, c = tableau
    @unpack wrapped_model, f, p = prob
    @unpack tₛ, U, V, model = wrapped_model
    t = model[:t]
    tsteps = supports(t)
    dt = tsteps[2] - tsteps[1]

    nᵤ = length(U)
    nᵥ = length(V)
    if MTK.is_explicit(tableau)
        K = Any[]
        for τ in tsteps[1:(end - 1)]
            for (i, h) in enumerate(c)
                ΔU = sum([A[i, j] * K[j] for j in 1:(i - 1)], init = zeros(nᵤ))
                Uₙ = [U[i](τ) + ΔU[i] * dt for i in 1:nᵤ]
                Vₙ = [V[i](τ) for i in 1:nᵥ]
                Kₙ = tₛ * f(Uₙ, Vₙ, p, τ + h * dt)
                push!(K, Kₙ)
            end
            ΔU = dt * sum([α[i] * K[i] for i in 1:length(α)])
            @constraint(model, [n = 1:nᵤ], U[n](τ) + ΔU[n]==U[n](τ + dt),
                base_name="solve_time_$τ")
            empty!(K)
        end
    else
        K = @variable(model, K[1:length(α), 1:nᵤ], Infinite(model[:t]))
        ΔUs = A * K
        ΔU_tot = dt * (K' * α)
        for τ in tsteps[1:(end - 1)]
            for (i, h) in enumerate(c)
                ΔU = @view ΔUs[i, :]
                Uₙ = U + ΔU * dt
                @constraint(model, [j = 1:nᵤ], K[i, j]==(tₛ * f(Uₙ, V, p, τ + h * dt)[j]),
                    DomainRestrictions(t => τ), base_name="solve_K$i($τ)")
            end
            @constraint(model,
                [n = 1:nᵤ], U[n](τ) + ΔU_tot[n]==U[n](min(τ + dt, tsteps[end])),
                DomainRestrictions(t => τ), base_name="solve_U($τ)")
        end
    end
end

struct JuMPCollocation <: AbstractCollocation
    solver::Any
    tableau::DiffEqBase.ODERKTableau
end
function MTK.JuMPCollocation(solver, tableau = MTK.constructDefault())
    JuMPCollocation(solver, tableau)
end

struct InfiniteOptCollocation <: AbstractCollocation
    solver::Any
    derivative_method::InfiniteOpt.AbstractDerivativeMethod
end
function MTK.InfiniteOptCollocation(
        solver, derivative_method = InfiniteOpt.FiniteDifference(InfiniteOpt.Backward()))
    InfiniteOptCollocation(solver, derivative_method)
end

function MTK.prepare_and_optimize!(
        prob::JuMPDynamicOptProblem, solver::JuMPCollocation; verbose = false, kwargs...)
    model = prob.wrapped_model.model
    verbose || set_silent(model)
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
            unregister(model, Symbol(JuMP.name(var)))
            delete(model, var)
        end
    end
    add_solve_constraints!(prob, solver.tableau)
    set_optimizer(model, solver.solver)
    optimize!(model)
    model
end

function MTK.prepare_and_optimize!(prob::InfiniteOptDynamicOptProblem,
        solver::InfiniteOptCollocation; verbose = false, kwargs...)
    model = prob.wrapped_model.model
    verbose || set_silent(model)
    set_derivative_method(model[:t], solver.derivative_method)
    set_optimizer(model, solver.solver)
    optimize!(model)
    model
end

function MTK.get_V_values(m::InfiniteModel)
    nt = length(supports(m[:t]))
    if !isempty(m[:V])
        V_vals = value.(m[:V])
        V_vals = [[V_vals[i][j] for i in 1:length(V_vals)] for j in 1:nt]
    else
        nothing
    end
end
function MTK.get_U_values(m::InfiniteModel)
    nt = length(supports(m[:t]))
    U_vals = value.(m[:U])
    U_vals = [[U_vals[i][j] for i in 1:length(U_vals)] for j in 1:nt]
end
MTK.get_t_values(m::InfiniteModel) = value(m[:tₛ]) * supports(m[:t])
MTK.objective_value(m::InfiniteModel) = InfiniteOpt.objective_value(m)

function MTK.successful_solve(model::InfiniteModel)
    tstatus = termination_status(model)
    pstatus = primal_status(model)
    !has_values(model) &&
        error("Model not solvable; please report this to github.com/SciML/ModelingToolkit.jl with a MWE.")

    pstatus === FEASIBLE_POINT &&
        (tstatus === OPTIMAL || tstatus === LOCALLY_SOLVED || tstatus === ALMOST_OPTIMAL ||
         tstatus === ALMOST_LOCALLY_SOLVED)
end

import InfiniteOpt: JuMP, GeneralVariableRef

for ff in [acos, log1p, acosh, log2, asin, tan, atanh, cos, log, sin, log10, sqrt]
    f = nameof(ff)
    # These need to be defined so that JuMP can trace through functions built by Symbolics
    @eval NaNMath.$f(x::GeneralVariableRef) = Base.$f(x)
end

# JuMP variables and Symbolics variables never compare equal. When tracing through dynamics, a function argument can be either a JuMP variable or A Symbolics variable, it can never be both.
function Base.isequal(::SymbolicUtils.Symbolic,
        ::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr, InfiniteOpt.AbstractInfOptExpr})
    false
end
function Base.isequal(
        ::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr, InfiniteOpt.AbstractInfOptExpr},
        ::SymbolicUtils.Symbolic)
    false
end
end
