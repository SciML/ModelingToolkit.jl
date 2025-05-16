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
    model::InfiniteOptModel
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
    model::InfiniteOptModel
    kwargs::K

    function InfiniteOptDynamicOptProblem(f, u0, tspan, p, model, kwargs...)
        new{typeof(u0), typeof(tspan), SciMLBase.isinplace(f),
            typeof(p), typeof(f), typeof(kwargs)}(f, u0, tspan, p, model, kwargs)
    end
end

MTK.generate_internal_model(m::Type{InfiniteOptModel}) = InfiniteModel()
MTK.generate_time_variable!(m::InfiniteModel, tspan, steps) = @infinite_parameter(m, t in [tspan[1], tspan[2]], num_supports = steps)
MTK.generate_state_variable!(m::InfiniteModel, u0::Vector, ns, nt) = @variable(m, U[i = 1:ns], Infinite(m[:t]), start=u0[i])
MTK.generate_input_variable!(m::InfiniteModel, c0, nc, nt) = @variable(m, V[i = 1:nc], Infinite(m[:t]), start=c0[i])

function MTK.generate_timescale!(m::InfiniteModel, guess, is_free_t)
    @variable(m, tₛ ≥ 0, start = guess)
    if !is_free_t
        fix(tₛ, 1, force=true)
        set_start_value(tₛ, 1)
    end
    tₛ
end

function MTK.add_constraint!(m::InfiniteOptModel, expr::Union{Equation, Inequality}) 
    if expr isa Equation
        @constraint(m.model, expr.lhs - expr.rhs == 0)
    elseif expr.relational_op === Symbolics.geq
        @constraint(m.model, expr.lhs - eq.rhs ≥ 0)
    else
        @constraint(m.model, expr.lhs - eq.rhs ≤ 0)
    end
end
MTK.set_objective!(m::InfiniteOptModel, expr) = @objective(m.model, Min, expr)

"""
    JuMPDynamicOptProblem(sys::System, u0, tspan, p; dt)

Convert a System representing an optimal control system into a JuMP model
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
function MTK.JuMPDynamicOptProblem(sys::System, u0map, tspan, pmap;
        dt = nothing,
        steps = nothing,
        guesses = Dict(), kwargs...)
    MTK.process_DynamicOptProblem(JuMPDynamicOptProblem, InfiniteOptModel, sys, u0map, tspan, pmap; dt, steps, guesses, kwargs...)
end

"""
    InfiniteOptDynamicOptProblem(sys::System, u0map, tspan, pmap; dt)

Convert System representing an optimal control system into a InfiniteOpt model
for solving using optimization. Must provide `dt` for determining the length 
of the interpolation arrays.

Related to `JuMPDynamicOptProblem`, but directly adds the differential equations
of the system as derivative constraints, rather than using a solver tableau.
"""
function MTK.InfiniteOptDynamicOptProblem(sys::System, u0map, tspan, pmap;
        dt = nothing,
        steps = nothing,
        guesses = Dict(), kwargs...)
    MTK.process_DynamicOptProblem(InfiniteOptDynamicOptProblem, InfiniteOptModel, sys, u0map, tspan, pmap; dt, steps, guesses, kwargs...)
end

function MTK.set_variable_bounds!(model, sys, pmap, tf)
    for (i, u) in enumerate(unknowns(sys))
        if MTK.hasbounds(u)
            lo, hi = MTK.getbounds(u)
            set_lower_bound(model.U[i], Symbolics.fixpoint_sub(lo, pmap))
            set_upper_bound(model.U[i], Symbolics.fixpoint_sub(hi, pmap))
        end
    end

    for (i, v) in enumerate(MTK.unbound_inputs(sys))
        if MTK.hasbounds(v)
            lo, hi = MTK.getbounds(v)
            set_lower_bound(model.V[i], Symbolics.fixpoint_sub(lo, pmap))
            set_upper_bound(model.V[i], Symbolics.fixpoint_sub(hi, pmap))
        end
    end

    if MTK.symbolic_type(tf) === MTK.ScalarSymbolic() && hasbounds(tf)
        lo, hi = MTK.getbounds(tf)
        set_lower_bound(model.tₛ, lo)
        set_upper_bound(model.tₛ, hi)
    end
end

function MTK.substitute_integral(model, exprs)
    intmap = Dict()
    for int in MTK.collect_applied_operators(exprs, Symbolics.Integral)
        op = MTK.operation(int)
        arg = only(arguments(MTK.value(int)))
        lo, hi = MTK.value.((op.domain.domain.left, op.domain.domain.right))
        hi = (MTK.symbolic_type(hi) === MTK.ScalarSymbolic()) ? 1 : hi
        intmap[int] = model.tₛ * InfiniteOpt.∫(arg, model.model[:t], lo, hi)
    end
    exprs = map(c -> Symbolics.substitute(c, intmap), exprs)
end

function MTK.add_initial_constraints!(m::InfiniteOptModel, u0, u0_idxs, ts)
    @show m.U
    @constraint(m.model, initial[i in u0_idxs], m.U[i](ts)==u0[i])
end

function MTK.substitute_model_vars(model, sys, exprs; tf = nothing)
    whole_interval_map = Dict([[v => model.U[i] for (i, v) in enumerate(unknowns(sys))];
                               [v => model.V[i] for (i, v) in enumerate(MTK.unbound_inputs(sys))]])
    exprs = map(c -> Symbolics.fast_substitute(c, whole_interval_map), exprs)

    x_ops = [MTK.operation(MTK.unwrap(st)) for st in unknowns(sys)]
    c_ops = [MTK.operation(MTK.unwrap(ct)) for ct in MTK.unbound_inputs(sys)]

    if MTK.symbolic_type(tf) === MTK.ScalarSymbolic()
        free_t_map = Dict([[x(tf) => model.U[i](1) for (i, x) in enumerate(x_ops)];
                           [c(tf) => model.V[i](1) for (i, c) in enumerate(c_ops)]])
        exprs = map(c -> Symbolics.fast_substitute(c, free_t_map), exprs)
    end

    # for variables like x(1.0)
    fixed_t_map = Dict([[x_ops[i] => model.U[i] for i in 1:length(model.U)];
                        [c_ops[i] => model.V[i] for i in 1:length(model.V)]])
    exprs = map(c -> Symbolics.fast_substitute(c, fixed_t_map), exprs)
end

function MTK.substitute_differentials(model::InfiniteOptModel, eqs)
    U = model.U
    t = model.model[:t]
    D = Differential(MTK.get_iv(sys))
    diffsubmap = Dict([D(U[i]) => ∂(U[i], t) for i in 1:length(U)])
    map(e -> Symbolics.substitute(e, diffsubmap), diff_eqs)
end

function add_solve_constraints!(prob::JuMPDynamicOptProblem, tableau)
    @unpack A, α, c = tableau
    @unpack model, f, p = prob
    tsteps = supports(model.model[:t])
    dt = tsteps[2] - tsteps[1]

    tₛ = model.tₛ
    U = model.U
    V = model.V
    nᵤ = length(U)
    nᵥ = length(V)
    if MTK.is_explicit(tableau)
        K = Any[]
        for τ in tsteps[1:end-1]
            for (i, h) in enumerate(c)
                ΔU = sum([A[i, j] * K[j] for j in 1:(i - 1)], init = zeros(nᵤ))
                Uₙ = [U[i](τ) + ΔU[i] * dt for i in 1:nᵤ]
                Vₙ = [V[i](τ) for i in 1:nᵥ]
                Kₙ = tₛ * f(Uₙ, Vₙ, p, τ + h * dt)
                push!(K, Kₙ)
            end
            ΔU = dt * sum([α[i] * K[i] for i in 1:length(α)])
            @constraint(model.model, [n = 1:nᵤ], U[n](τ) + ΔU[n]==U[n](τ + dt),
                base_name="solve_time_$τ")
            empty!(K)
        end
    else
        @variable(model, K[1:length(α), 1:nᵤ], Infinite(t))
        ΔUs = A * K
        ΔU_tot = dt * (K' * α)
        for τ in tsteps[1:end-1]
            for (i, h) in enumerate(c)
                ΔU = @view ΔUs[i, :]
                Uₙ = U + ΔU * h * dt
                @constraint(model.model, [j = 1:nᵤ], K[i, j]==(tₛ * f(Uₙ, V, p, τ + h * dt)[j]),
                    DomainRestrictions(t => τ), base_name="solve_K$i($τ)")
            end
            @constraint(model.model, [n = 1:nᵤ], U[n](τ) + ΔU_tot[n]==U[n](min(τ + dt, tsteps[end])),
                DomainRestrictions(t => τ), base_name="solve_U($τ)")
        end
    end
end

"""
JuMP Collocation solver.
- solver: a optimization solver such as Ipopt
- tableau: An ODE RK tableau. Load a tableau by calling a function like `constructRK4` and may be found at https://docs.sciml.ai/DiffEqDevDocs/stable/internals/tableaus/. If this argument is not passed in, the solver will default to Radau second order.

Returns a DynamicOptSolution, which contains both the model and the ODE solution.
"""
struct JuMPCollocation <: AbstractCollocation
    solver::Any
    tableau::DiffEqBase.ODERKTableau
end
MTK.JuMPCollocation(solver, tableau = MTK.constructDefault()) = JuMPCollocation(solver, tableau)

"""
InfiniteOpt Collocation solver.
- solver: an optimization solver such as Ipopt
- `derivative_method` kwarg refers to the method used by InfiniteOpt to compute derivatives. The list of possible options can be found at https://infiniteopt.github.io/InfiniteOpt.jl/stable/guide/derivative/. Defaults to FiniteDifference(Backward()).
"""
struct InfiniteOptCollocation <: AbstractCollocation
    solver::Any
    derivative_method::InfiniteOpt.AbstractDerivativeMethod
end
MTK.InfiniteOptCollocation(solver, derivative_method = InfiniteOpt.FiniteDifference(InfiniteOpt.Backward())) = InfiniteOptCollocation(solver, derivative_method)

function MTK.prepare_solver!(prob::JuMPDynamicOptProblem, solver::JuMPCollocation; verbose = false, kwargs...)
    model = prob.model.model
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
end

function MTK.prepare_solver!(prob::InfiniteOptDynamicOptProblem, solver::InfiniteOptCollocation; verbose = false, kwargs...)
    model = prob.model.model
    verbose || set_silent(model)
    add_equational_constraints!(model, prob.f.sys, prob.tspan)
    set_derivative_method(model[:t], solver.derivative_method)
    set_optimizer(model, solver.solver)
end

function MTK.optimize_model!(prob::Union{InfiniteOptDynamicOptProblem, JuMPDynamicOptProblem}, solver) 
    optimize!(prob.model.model)
    prob.model
end

function MTK.get_V_values(m::InfiniteOptModel)
    nt = length(supports(m.model[:t]))
    if !isempty(m.V)
        V_vals = value.(m.V)
        V_vals = [[V_vals[i][j] for i in 1:length(V_vals)] for j in 1:nt]
    else
        nothing
    end
end
function MTK.get_U_values(m::InfiniteOptModel)
    nt = length(supports(m.model[:t]))
    U_vals = value.(m.U)
    U_vals = [[U_vals[i][j] for i in 1:length(U_vals)] for j in 1:nt]
end
MTK.get_t_values(model) = model.tₛ * supports(model.model[:t])

function MTK.successful_solve(m::InfiniteOptModel)
    model = m.model
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
