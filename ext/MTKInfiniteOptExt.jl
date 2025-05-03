module MTKInfiniteOptExt
using ModelingToolkit
using InfiniteOpt
using DiffEqBase
using LinearAlgebra
using StaticArrays
import SymbolicUtils
import NaNMath
const MTK = ModelingToolkit

struct JuMPDynamicOptProblem{uType, tType, isinplace, P, F, K} <:
       AbstractDynamicOptProblem{uType, tType, isinplace}
    f::F
    u0::uType
    tspan::tType
    p::P
    model::InfiniteModel
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
    model::InfiniteModel
    kwargs::K

    function InfiniteOptDynamicOptProblem(f, u0, tspan, p, model, kwargs...)
        new{typeof(u0), typeof(tspan), SciMLBase.isinplace(f),
            typeof(p), typeof(f), typeof(kwargs)}(f, u0, tspan, p, model, kwargs)
    end
end

"""
    JuMPDynamicOptProblem(sys::ODESystem, u0, tspan, p; dt)

Convert an ODESystem representing an optimal control system into a JuMP model
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
function MTK.JuMPDynamicOptProblem(sys::ODESystem, u0map, tspan, pmap;
        dt = nothing,
        steps = nothing,
        guesses = Dict(), kwargs...)
    MTK.warn_overdetermined(sys, u0map)
    _u0map = has_alg_eqs(sys) ? u0map : merge(Dict(u0map), Dict(guesses))
    f, u0, p = MTK.process_SciMLProblem(ODEInputFunction, sys, _u0map, pmap;
        t = tspan !== nothing ? tspan[1] : tspan, kwargs...)

    pmap = Dict{Any, Any}(pmap)
    steps, is_free_t = MTK.process_tspan(tspan, dt, steps)
    model = init_model(sys, tspan, steps, u0map, pmap, u0; is_free_t)

    JuMPDynamicOptProblem(f, u0, tspan, p, model, kwargs...)
end

"""
    InfiniteOptDynamicOptProblem(sys::ODESystem, u0map, tspan, pmap; dt)

Convert an ODESystem representing an optimal control system into a InfiniteOpt model
for solving using optimization. Must provide `dt` for determining the length 
of the interpolation arrays.

Related to `JuMPDynamicOptProblem`, but directly adds the differential equations
of the system as derivative constraints, rather than using a solver tableau.
"""
function MTK.InfiniteOptDynamicOptProblem(sys::ODESystem, u0map, tspan, pmap;
        dt = nothing,
        steps = nothing,
        guesses = Dict(), kwargs...)
    MTK.warn_overdetermined(sys, u0map)
    _u0map = has_alg_eqs(sys) ? u0map : merge(Dict(u0map), Dict(guesses))
    f, u0, p = MTK.process_SciMLProblem(ODEInputFunction, sys, _u0map, pmap;
        t = tspan !== nothing ? tspan[1] : tspan, kwargs...)

    pmap = Dict{Any, Any}(pmap)
    steps, is_free_t = MTK.process_tspan(tspan, dt, steps)
    model = init_model(sys, tspan, steps, u0map, pmap, u0; is_free_t)

    add_infopt_solve_constraints!(model, sys, pmap; is_free_t)
    InfiniteOptDynamicOptProblem(f, u0, tspan, p, model, kwargs...)
end

# Initialize InfiniteOpt model.
function init_model(sys, tspan, steps, u0map, pmap, u0; is_free_t = false)
    ctrls = MTK.unbound_inputs(sys)
    states = unknowns(sys)
    model = InfiniteModel()

    if is_free_t
        (ts_sym, te_sym) = tspan
        MTK.symbolic_type(ts_sym) !== MTK.NotSymbolic() &&
            error("Free initial time problems are not currently supported.")
        @variable(model, tf, start=pmap[te_sym])
        set_lower_bound(tf, ts_sym)
        hasbounds(te_sym) && begin
            lo, hi = getbounds(te_sym)
            set_lower_bound(tf, lo)
            set_upper_bound(tf, hi)
        end
        pmap[te_sym] = model[:tf]
        tspan = (0, 1)
    end

    @infinite_parameter(model, t in [tspan[1], tspan[2]], num_supports=steps)
    @variable(model, U[i = 1:length(states)], Infinite(t), start=u0[i])
    c0 = MTK.value.([pmap[c] for c in ctrls])
    @variable(model, V[i = 1:length(ctrls)], Infinite(t), start=c0[i])
    for (i, ct) in enumerate(ctrls)
        pmap[ct] = model[:V][i]
    end

    set_jump_bounds!(model, sys, pmap)
    add_jump_cost_function!(model, sys, (tspan[1], tspan[2]), pmap; is_free_t)
    add_user_constraints!(model, sys, pmap; is_free_t)

    stidxmap = Dict([v => i for (i, v) in enumerate(states)])
    u0map = Dict([MTK.default_toterm(MTK.value(k)) => v for (k, v) in u0map])
    u0_idxs = has_alg_eqs(sys) ? collect(1:length(states)) :
              [stidxmap[MTK.default_toterm(k)] for (k, v) in u0map]
    add_initial_constraints!(model, u0, u0_idxs, tspan[1])
    return model
end

"""
Modify the pmap by replacing controls with V[i](t), and tf with the model's final time variable for free final time problems.
"""
function modified_pmap(model, sys, pmap)
    pmap = Dict{Any, Any}(pmap)
end

function set_jump_bounds!(model, sys, pmap)
    U = model[:U]
    for (i, u) in enumerate(unknowns(sys))
        if MTK.hasbounds(u)
            lo, hi = MTK.getbounds(u)
            set_lower_bound(U[i], Symbolics.fixpoint_sub(lo, pmap))
            set_upper_bound(U[i], Symbolics.fixpoint_sub(hi, pmap))
        end
    end

    V = model[:V]
    for (i, v) in enumerate(MTK.unbound_inputs(sys))
        if MTK.hasbounds(v)
            lo, hi = MTK.getbounds(v)
            set_lower_bound(V[i], Symbolics.fixpoint_sub(lo, pmap))
            set_upper_bound(V[i], Symbolics.fixpoint_sub(hi, pmap))
        end
    end
end

function add_jump_cost_function!(model::InfiniteModel, sys, tspan, pmap; is_free_t = false)
    jcosts = MTK.get_costs(sys)
    consolidate = MTK.get_consolidate(sys)
    if isnothing(jcosts) || isempty(jcosts)
        @objective(model, Min, 0)
        return
    end
    jcosts = substitute_jump_vars(model, sys, pmap, jcosts; is_free_t)
    tₛ = is_free_t ? model[:tf] : 1

    # Substitute integral
    iv = MTK.get_iv(sys)

    intmap = Dict()
    for int in MTK.collect_applied_operators(jcosts, Symbolics.Integral)
        op = MTK.operation(int)
        arg = only(arguments(MTK.value(int)))
        lo, hi = (op.domain.domain.left, op.domain.domain.right)
        lo = MTK.value(lo)
        hi = haskey(pmap, hi) ? 1 : MTK.value(hi)
        intmap[int] = tₛ * InfiniteOpt.∫(arg, model[:t], lo, hi)
    end
    jcosts = map(c -> Symbolics.substitute(c, intmap), jcosts)
    @objective(model, Min, consolidate(jcosts))
end

function add_user_constraints!(model::InfiniteModel, sys, pmap; is_free_t = false)
    conssys = MTK.get_constraintsystem(sys)
    jconstraints = isnothing(conssys) ? nothing : MTK.get_constraints(conssys)
    (isnothing(jconstraints) || isempty(jconstraints)) && return nothing

    if is_free_t
        for u in MTK.get_unknowns(conssys)
            x = MTK.operation(u)
            t = only(arguments(u))
            if (MTK.symbolic_type(t) === MTK.NotSymbolic())
                error("Provided specific time constraint in a free final time problem. This is not supported by the JuMP/InfiniteOpt collocation solvers. The offending variable is $u. Specific-time constraints can only be specified at the beginning or end of the timespan.")
            end
        end
    end

    auxmap = Dict([u => MTK.default_toterm(MTK.value(u)) for u in unknowns(conssys)])
    jconstraints = substitute_jump_vars(model, sys, pmap, jconstraints; auxmap, is_free_t)

    # Substitute to-term'd variables
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

function substitute_jump_vars(model, sys, pmap, exprs; auxmap = Dict(), is_free_t = false)
    iv = MTK.get_iv(sys)
    sts = unknowns(sys)
    cts = MTK.unbound_inputs(sys)
    U = model[:U]
    V = model[:V]
    x_ops = [MTK.operation(MTK.unwrap(st)) for st in sts]
    c_ops = [MTK.operation(MTK.unwrap(ct)) for ct in cts]

    exprs = map(c -> Symbolics.fixpoint_sub(c, auxmap), exprs)
    exprs = map(c -> Symbolics.fixpoint_sub(c, Dict(pmap)), exprs)
    if is_free_t
        tf = model[:tf]
        free_t_map = Dict([[x(tf) => U[i](1) for (i, x) in enumerate(x_ops)];
                           [c(tf) => V[i](1) for (i, c) in enumerate(c_ops)]])
        exprs = map(c -> Symbolics.fixpoint_sub(c, free_t_map), exprs)
    end

    # for variables like x(t)
    whole_interval_map = Dict([[v => U[i] for (i, v) in enumerate(sts)];
                               [v => V[i] for (i, v) in enumerate(cts)]])
    exprs = map(c -> Symbolics.fixpoint_sub(c, whole_interval_map), exprs)

    # for variables like x(1.0)
    fixed_t_map = Dict([[x_ops[i] => U[i] for i in 1:length(U)];
                        [c_ops[i] => V[i] for i in 1:length(V)]])

    exprs = map(c -> Symbolics.fixpoint_sub(c, fixed_t_map), exprs)
    exprs
end

is_explicit(tableau) = tableau isa DiffEqBase.ExplicitRKTableau

function add_infopt_solve_constraints!(model::InfiniteModel, sys, pmap; is_free_t = false)
    # Differential equations
    U = model[:U]
    t = model[:t]
    D = Differential(MTK.get_iv(sys))
    diffsubmap = Dict([D(U[i]) => ∂(U[i], t) for i in 1:length(U)])
    tₛ = is_free_t ? model[:tf] : 1

    diff_eqs = substitute_jump_vars(model, sys, pmap, diff_equations(sys))
    diff_eqs = map(e -> Symbolics.substitute(e, diffsubmap), diff_eqs)
    @constraint(model, D[i = 1:length(diff_eqs)], diff_eqs[i].lhs==tₛ * diff_eqs[i].rhs)

    # Algebraic equations
    alg_eqs = substitute_jump_vars(model, sys, pmap, alg_equations(sys))
    @constraint(model, A[i = 1:length(alg_eqs)], alg_eqs[i].lhs==alg_eqs[i].rhs)
end

function add_jump_solve_constraints!(prob, tableau; is_free_t = false)
    A = tableau.A
    α = tableau.α
    c = tableau.c
    model = prob.model
    f = prob.f
    p = prob.p

    t = model[:t]
    tsteps = supports(t)
    tmax = tsteps[end]
    pop!(tsteps)
    tₛ = is_free_t ? model[:tf] : 1
    dt = tsteps[2] - tsteps[1]

    U = model[:U]
    V = model[:V]
    nᵤ = length(U)
    nᵥ = length(V)
    if is_explicit(tableau)
        K = Any[]
        for τ in tsteps
            for (i, h) in enumerate(c)
                ΔU = sum([A[i, j] * K[j] for j in 1:(i - 1)], init = zeros(nᵤ))
                Uₙ = [U[i](τ) + ΔU[i] * dt for i in 1:nᵤ]
                Vₙ = [V[i](τ) for i in 1:nᵥ]
                Kₙ = tₛ * f(Uₙ, Vₙ, p, τ + h * dt) # scale the time
                push!(K, Kₙ)
            end
            ΔU = dt * sum([α[i] * K[i] for i in 1:length(α)])
            @constraint(model, [n = 1:nᵤ], U[n](τ) + ΔU[n]==U[n](τ + dt),
                base_name="solve_time_$τ")
            empty!(K)
        end
    else
        @variable(model, K[1:length(α), 1:nᵤ], Infinite(t))
        ΔUs = A * K
        ΔU_tot = dt * (K' * α)
        for τ in tsteps
            for (i, h) in enumerate(c)
                ΔU = @view ΔUs[i, :]
                Uₙ = U + ΔU * h * dt
                @constraint(model, [j = 1:nᵤ], K[i, j]==(tₛ * f(Uₙ, V, p, τ + h * dt)[j]),
                    DomainRestrictions(t => τ), base_name="solve_K$i($τ)")
            end
            @constraint(model, [n = 1:nᵤ], U[n](τ) + ΔU_tot[n]==U[n](min(τ + dt, tmax)),
                DomainRestrictions(t => τ), base_name="solve_U($τ)")
        end
    end
end

"""
Default ODE Tableau: RadauIIA5
"""
function constructDefault(T::Type = Float64)
    sq6 = sqrt(6)
    A = [11 // 45-7sq6 / 360 37 // 225-169sq6 / 1800 -2 // 225+sq6 / 75
         37 // 225+169sq6 / 1800 11 // 45+7sq6 / 360 -2 // 225-sq6 / 75
         4 // 9-sq6 / 36 4 // 9+sq6 / 36 1//9]
    c = [2 // 5 - sq6 / 10; 2 / 5 + sq6 / 10; 1]
    α = [4 // 9 - sq6 / 36; 4 // 9 + sq6 / 36; 1 // 9]
    A = map(T, A)
    α = map(T, α)
    c = map(T, c)
    
    DiffEqBase.ImplicitRKTableau(A, c, α, 5)
end

"""
Solve JuMPDynamicOptProblem. Arguments:
- prob: a JumpDynamicOptProblem
- jump_solver: a LP solver such as HiGHS
- ode_solver: Takes in a symbol representing the solver. Acceptable solvers may be found at https://docs.sciml.ai/DiffEqDevDocs/stable/internals/tableaus/. Note that the symbol may be different than the typical name of the solver, e.g. :Tsitouras5 rather than Tsit5.
- silent: set the model silent (suppress model output)

Returns a DynamicOptSolution, which contains both the model and the ODE solution.
"""
function DiffEqBase.solve(
        prob::JuMPDynamicOptProblem, jump_solver, ode_solver::Symbol = :Default; silent = false)
    model = prob.model
    tableau_getter = Symbol(:construct, ode_solver)
    if ode_solver == :Default
        @eval tableau = $tableau_getter()
    else
        @eval tableau = Main.$tableau_getter()
    end
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
            unregister(model, Symbol(JuMP.name(var)))
            delete(model, var)
        end
    end
    add_jump_solve_constraints!(prob, tableau; is_free_t = haskey(model, :tf))
    _solve(prob, jump_solver, ode_solver)
end

"""
`derivative_method` kwarg refers to the method used by InfiniteOpt to compute derivatives. The list of possible options can be found at https://infiniteopt.github.io/InfiniteOpt.jl/stable/guide/derivative/. Defaults to FiniteDifference(Backward()).
"""
function DiffEqBase.solve(prob::InfiniteOptDynamicOptProblem, jump_solver;
        derivative_method = InfiniteOpt.FiniteDifference(Backward()), silent = false)
    model = prob.model
    silent && set_silent(model)
    set_derivative_method(model[:t], derivative_method)
    _solve(prob, jump_solver, derivative_method)
end

function _solve(prob::AbstractDynamicOptProblem, jump_solver, solver)
    model = prob.model
    set_optimizer(model, jump_solver)
    optimize!(model)

    tstatus = termination_status(model)
    pstatus = primal_status(model)
    !has_values(model) &&
        error("Model not solvable; please report this to github.com/SciML/ModelingToolkit.jl with a MWE.")

    tf = haskey(model, :tf) ? value(model[:tf]) : 1
    ts = tf * supports(model[:t])
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

    DynamicOptSolution(model, sol, input_sol)
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
