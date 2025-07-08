abstract type AbstractDynamicOptProblem{uType, tType, isinplace} <:
              SciMLBase.AbstractODEProblem{uType, tType, isinplace} end

abstract type AbstractCollocation end

struct DynamicOptSolution
    model::Any
    sol::ODESolution
    input_sol::Union{Nothing, ODESolution}
end

function Base.show(io::IO, sol::DynamicOptSolution)
    println("retcode: ", sol.sol.retcode, "\n")

    println("Optimal control solution for following model:\n")
    show(sol.model)

    print("\n\nPlease query the model using sol.model, the solution trajectory for the system using sol.sol, or the solution trajectory for the controllers using sol.input_sol.")
end

"""
    JuMPDynamicOptProblem(sys::System, op, tspan; dt, steps, guesses, kwargs...)

Convert an System representing an optimal control system into a JuMP model
for solving using optimization. Must provide either `dt`, the timestep between collocation 
points (which, along with the timespan, determines the number of points), or directly 
provide the number of points as `steps`.

To construct the problem, please load InfiniteOpt along with ModelingToolkit.
"""
function JuMPDynamicOptProblem end
"""
    InfiniteOptDynamicOptProblem(sys::System, op, tspan; dt)

Convert an System representing an optimal control system into a InfiniteOpt model
for solving using optimization. Must provide `dt` for determining the length 
of the interpolation arrays.

Related to `JuMPDynamicOptProblem`, but directly adds the differential equations
of the system as derivative constraints, rather than using a solver tableau.

To construct the problem, please load InfiniteOpt along with ModelingToolkit.
"""
function InfiniteOptDynamicOptProblem end
"""
    CasADiDynamicOptProblem(sys::System, op, tspan; dt, steps, guesses, kwargs...)

Convert an System representing an optimal control system into a CasADi model
for solving using optimization. Must provide either `dt`, the timestep between collocation 
points (which, along with the timespan, determines the number of points), or directly 
provide the number of points as `steps`.

To construct the problem, please load CasADi along with ModelingToolkit.
"""
function CasADiDynamicOptProblem end
"""
    PyomoDynamicOptProblem(sys::System, op, tspan; dt, steps)

Convert an System representing an optimal control system into a Pyomo model
for solving using optimization. Must provide either `dt`, the timestep between collocation 
points (which, along with the timespan, determines the number of points), or directly 
provide the number of points as `steps`.

To construct the problem, please load Pyomo along with ModelingToolkit.
"""
function PyomoDynamicOptProblem end

### Collocations
"""
JuMP Collocation solver. Takes two arguments:
- `solver`: a optimization solver such as Ipopt
- `tableau`: An ODE RK tableau. Load a tableau by calling a function like `constructRK4` and may be found at https://docs.sciml.ai/DiffEqDevDocs/stable/internals/tableaus/. If this argument is not passed in, the solver will default to Radau second order.
"""
function JuMPCollocation end
"""
InfiniteOpt Collocation solver.
- `solver`: an optimization solver such as Ipopt
- `derivative_method`: the method used by InfiniteOpt to compute derivatives. The list of possible options can be found at https://infiniteopt.github.io/InfiniteOpt.jl/stable/guide/derivative/. Defaults to FiniteDifference(Backward()).
"""
function InfiniteOptCollocation end
"""
CasADi Collocation solver.
- `solver`: an optimization solver such as Ipopt. Should be given as a string or symbol in all lowercase, e.g. "ipopt"
- `tableau`: An ODE RK tableau. Load a tableau by calling a function like `constructRK4` and may be found at https://docs.sciml.ai/DiffEqDevDocs/stable/internals/tableaus/. If this argument is not passed in, the solver will default to Radau second order.
"""
function CasADiCollocation end
"""
Pyomo Collocation solver.
- `solver`: an optimization solver such as Ipopt. Should be given as a string or symbol in all lowercase, e.g. "ipopt"
- `derivative_method`: a derivative method from Pyomo. The choices here are ForwardEuler, BackwardEuler, MidpointEuler, LagrangeRadau, or LagrangeLegendre. The last two should additionally have a number indicating the number of collocation points per timestep, e.g. PyomoCollocation("ipopt", LagrangeRadau(3)). Defaults to LagrangeRadau(5).
"""
function PyomoCollocation end

function warn_overdetermined(sys, op)
    cstrs = constraints(sys)
    init_conds = filter(x -> value(x) ∈ Set(unknowns(sys)), [k for (k, v) in op])
    if !isempty(cstrs)
        (length(cstrs) + length(init_conds) > length(unknowns(sys))) &&
            @warn "The control problem is overdetermined. The total number of conditions (# constraints + # fixed initial values given by op) exceeds the total number of states. The solvers will default to doing a nonlinear least-squares optimization."
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

is_explicit(tableau) = tableau isa DiffEqBase.ExplicitRKTableau

@fallback_iip_specialize function SciMLBase.ODEInputFunction{iip, specialize}(sys::System;
        inputs = unbound_inputs(sys),
        disturbance_inputs = disturbances(sys),
        u0 = nothing, tgrad = false,
        jac = false, controljac = false,
        p = nothing, t = nothing,
        eval_expression = false,
        sparse = false, simplify = false,
        eval_module = @__MODULE__,
        steady_state = false,
        checkbounds = false,
        sparsity = false,
        analytic = nothing,
        initialization_data = nothing,
        cse = true,
        kwargs...) where {iip, specialize}
    f, _,
    _ = generate_control_function(
        sys, inputs, disturbance_inputs; eval_module, cse, kwargs...)
    f = f[1]

    if tgrad
        _tgrad = generate_tgrad(sys;
            simplify = simplify,
            expression = Val{true},
            wrap_gfw = Val{true},
            expression_module = eval_module, cse,
            checkbounds = checkbounds, kwargs...)
    else
        _tgrad = nothing
    end

    if jac
        _jac = generate_jacobian(sys;
            simplify = simplify, sparse = sparse,
            expression = Val{true},
            wrap_gfw = Val{true},
            expression_module = eval_module, cse,
            checkbounds = checkbounds, kwargs...)
    else
        _jac = nothing
    end

    if controljac
        _cjac = generate_control_jacobian(sys;
            simplify = simplify, sparse = sparse,
            expression = Val{true}, wrap_gfw = Val{true},
            expression_module = eval_module, cse,
            checkbounds = checkbounds, kwargs...)
    else
        _cjac = nothing
    end

    M = calculate_massmatrix(sys)
    _M = concrete_massmatrix(M; sparse, u0)

    observedfun = ObservedFunctionCache(
        sys; steady_state, eval_expression, eval_module, checkbounds, cse)

    _W_sparsity = W_sparsity(sys)
    W_prototype = calculate_W_prototype(_W_sparsity; u0, sparse)
    if sparse
        uElType = u0 === nothing ? Float64 : eltype(u0)
        controljac_prototype = similar(calculate_control_jacobian(sys), uElType)
    else
        controljac_prototype = nothing
    end

    ODEInputFunction{iip, specialize}(f;
        sys = sys,
        jac = _jac === nothing ? nothing : _jac,
        controljac = _cjac === nothing ? nothing : _cjac,
        tgrad = _tgrad === nothing ? nothing : _tgrad,
        mass_matrix = _M,
        jac_prototype = W_prototype,
        controljac_prototype = controljac_prototype,
        observed = observedfun,
        sparsity = sparsity ? _W_sparsity : nothing,
        analytic = analytic,
        initialization_data)
end

# returns the JuMP timespan, the number of steps, and whether it is a free time problem.
function process_tspan(tspan, dt, steps)
    is_free_time = false
    symbolic_type(tspan[1]) !== NotSymbolic() &&
        error("Free initial time problems are not currently supported by the collocation solvers.")

    if isnothing(dt) && isnothing(steps)
        error("Must provide either the dt or the number of intervals to the collocation solvers (JuMP, InfiniteOpt, CasADi).")
    elseif symbolic_type(tspan[1]) === ScalarSymbolic() ||
           symbolic_type(tspan[2]) === ScalarSymbolic()
        isnothing(steps) &&
            error("Free final time problems require specifying the number of steps using the keyword arg `steps`, rather than dt.")
        isnothing(dt) ||
            @warn "Specified dt for free final time problem. This will be ignored; dt will be determined by the number of timesteps."

        return (0, 1), steps, true
    else
        isnothing(steps) ||
            @warn "Specified number of steps for problem with concrete tspan. This will be ignored; number of steps will be determined by dt."

        return tspan, length(tspan[1]:dt:tspan[2]), false
    end
end

##########################
### MODEL CONSTRUCTION ###
##########################
function process_DynamicOptProblem(
        prob_type::Type{<:AbstractDynamicOptProblem}, model_type, sys::System, op, tspan;
        dt = nothing,
        steps = nothing,
        guesses = Dict(), kwargs...)
    warn_overdetermined(sys, op)
    ctrls = unbound_inputs(sys)
    states = unknowns(sys)

    stidxmap = Dict([v => i for (i, v) in enumerate(states)])
    op = Dict([default_toterm(value(k)) => v for (k, v) in op])
    u0_idxs = has_alg_eqs(sys) ? collect(1:length(states)) :
              [stidxmap[default_toterm(k)] for (k, v) in op if haskey(stidxmap, k)]

    _op = has_alg_eqs(sys) ? op : merge(Dict(op), Dict(guesses))
    f, u0,
    p = process_SciMLProblem(ODEInputFunction, sys, _op;
        t = tspan !== nothing ? tspan[1] : tspan, kwargs...)
    model_tspan, steps, is_free_t = process_tspan(tspan, dt, steps)
    warn_overdetermined(sys, op)

    pmap = filter(p -> (first(p) ∉ Set(unknowns(sys))), op)
    pmap = recursive_unwrap(AnyDict(pmap))
    evaluate_varmap!(pmap, keys(pmap))
    c0 = value.([pmap[c] for c in ctrls])

    tsteps = LinRange(model_tspan[1], model_tspan[2], steps)
    model = generate_internal_model(model_type)
    generate_time_variable!(model, model_tspan, tsteps)
    U = generate_state_variable!(model, u0, length(states), tsteps)
    V = generate_input_variable!(model, c0, length(ctrls), tsteps)
    tₛ = generate_timescale!(model, get(pmap, tspan[2], tspan[2]), is_free_t)
    fullmodel = model_type(model, U, V, tₛ, is_free_t)

    set_variable_bounds!(fullmodel, sys, pmap, tspan[2])
    add_cost_function!(fullmodel, sys, tspan, pmap)
    add_user_constraints!(fullmodel, sys, tspan, pmap)
    add_initial_constraints!(fullmodel, u0, u0_idxs, model_tspan[1])

    prob_type(f, u0, tspan, p, fullmodel, kwargs...), pmap
end

function generate_time_variable! end
function generate_internal_model end
function generate_state_variable! end
function generate_input_variable! end
function generate_timescale! end
function add_initial_constraints! end
function add_constraint! end

function set_variable_bounds!(m, sys, pmap, tf)
    @unpack model, U, V, tₛ = m
    t = get_iv(sys)
    for (i, u) in enumerate(unknowns(sys))
        var = lowered_var(m, :U, i, t)
        if hasbounds(u)
            lo, hi = getbounds(u)
            add_constraint!(m, var ≳ Symbolics.fixpoint_sub(lo, pmap))
            add_constraint!(m, var ≲ Symbolics.fixpoint_sub(hi, pmap))
        end
    end
    for (i, v) in enumerate(unbound_inputs(sys))
        var = lowered_var(m, :V, i, t)
        if hasbounds(v)
            lo, hi = getbounds(v)
            add_constraint!(m, var ≳ Symbolics.fixpoint_sub(lo, pmap))
            add_constraint!(m, var ≲ Symbolics.fixpoint_sub(hi, pmap))
        end
    end
    if symbolic_type(tf) === ScalarSymbolic() && hasbounds(tf)
        lo, hi = getbounds(tf)
        set_lower_bound(tₛ, Symbolics.fixpoint_sub(lo, pmap))
        set_upper_bound(tₛ, Symbolics.fixpoint_sub(hi, pmap))
    end
end

is_free_final(model) = model.is_free_final

function add_cost_function!(model, sys, tspan, pmap)
    jcosts = cost(sys)
    if Symbolics._iszero(jcosts)
        set_objective!(model, 0)
        return
    end

    jcosts = substitute_model_vars(model, sys, [jcosts], tspan)
    jcosts = substitute_params(pmap, jcosts)
    jcosts = substitute_integral(model, only(jcosts), tspan)
    set_objective!(model, value(jcosts))
end

"""
Substitute integrals. For an integral from (ts, te):
- Free final time problems should transcribe this to (0, 1) in the case that (ts, te) is the original timespan. Free final time problems cannot handle partial timespans.
- CasADi cannot handle partial timespans, even for non-free-final time problems.
time problems and unchanged otherwise.
"""
function substitute_integral(model, expr, tspan)
    intmap = Dict()
    for int in collect_applied_operators(expr, Symbolics.Integral)
        op = operation(int)
        arg = only(arguments(value(int)))
        lo, hi = value.((op.domain.domain.left, op.domain.domain.right))
        lo, hi = process_integral_bounds(model, (lo, hi), tspan)
        intmap[int] = lowered_integral(model, arg, lo, hi)
    end
    Symbolics.substitute(expr, intmap)
end

function process_integral_bounds(model, integral_span, tspan)
    if is_free_final(model) && isequal(integral_span, tspan)
        integral_span = (0, 1)
    elseif is_free_final(model)
        error("Free final time problems cannot handle partial timespans.")
    else
        (lo, hi) = integral_span
        (lo < tspan[1] || hi > tspan[2]) &&
            error("Integral bounds are beyond the timespan.")
        integral_span
    end
end

"""Substitute variables like x(1.5), x(t), etc. with the corresponding model variables."""
function substitute_model_vars(model, sys, exprs, tspan)
    x_ops = [operation(unwrap(st)) for st in unknowns(sys)]
    c_ops = [operation(unwrap(ct)) for ct in unbound_inputs(sys)]
    t = get_iv(sys)

    exprs = map(
        c -> Symbolics.fast_substitute(c, whole_t_map(model, t, x_ops, c_ops)), exprs)

    (ti, tf) = tspan
    if symbolic_type(tf) === ScalarSymbolic()
        _tf = model.tₛ + ti
        exprs = map(
            c -> Symbolics.fast_substitute(c, free_t_map(model, tf, x_ops, c_ops)), exprs)
        exprs = map(c -> Symbolics.fast_substitute(c, Dict(tf => _tf)), exprs)
    end
    exprs = map(c -> Symbolics.fast_substitute(c, fixed_t_map(model, x_ops, c_ops)), exprs)
    exprs
end

"""Mappings for variables that depend on the final time parameter, x(tf)."""
function free_t_map(m, tf, x_ops, c_ops)
    Dict([[x(tf) => lowered_var(m, :U, i, 1) for (i, x) in enumerate(x_ops)];
          [c(tf) => lowered_var(m, :V, i, 1) for (i, c) in enumerate(c_ops)]])
end

"""Mappings for variables that cover the whole timespan, x(t)."""
function whole_t_map(m, t, x_ops, c_ops)
    Dict([[v(t) => lowered_var(m, :U, i, t) for (i, v) in enumerate(x_ops)];
          [v(t) => lowered_var(m, :V, i, t) for (i, v) in enumerate(c_ops)]])
end

"""Mappings for variables that cover the whole timespan, x(t)."""
function fixed_t_map(m, x_ops, c_ops)
    Dict([[v => (t -> lowered_var(m, :U, i, t)) for (i, v) in enumerate(x_ops)];
          [v => (t -> lowered_var(m, :V, i, t)) for (i, v) in enumerate(c_ops)]])
end

function process_integral_bounds end
function lowered_integral end
function lowered_derivative end
function lowered_var end
function fixed_t_map end

function add_user_constraints!(model, sys, tspan, pmap)
    jconstraints = get_constraints(sys)
    (isnothing(jconstraints) || isempty(jconstraints)) && return nothing
    cons_dvs,
    cons_ps = process_constraint_system(
        jconstraints, Set(unknowns(sys)), parameters(sys), get_iv(sys); validate = false)

    is_free_final(model) && check_constraint_vars(cons_dvs)

    jconstraints = substitute_toterm(cons_dvs, jconstraints)
    jconstraints = substitute_model_vars(model, sys, jconstraints, tspan)
    jconstraints = substitute_params(pmap, jconstraints)

    for c in jconstraints
        add_constraint!(model, c)
    end
end

function add_equational_constraints!(model, sys, pmap, tspan)
    diff_eqs = substitute_model_vars(model, sys, diff_equations(sys), tspan)
    diff_eqs = substitute_params(pmap, diff_eqs)
    diff_eqs = substitute_differentials(model, sys, diff_eqs)
    for eq in diff_eqs
        add_constraint!(model, eq.lhs ~ eq.rhs * model.tₛ)
    end

    alg_eqs = substitute_model_vars(model, sys, alg_equations(sys), tspan)
    alg_eqs = substitute_params(pmap, alg_eqs)
    for eq in alg_eqs
        add_constraint!(model, eq.lhs ~ eq.rhs)
    end
end

function set_objective! end
objective_value(sol::DynamicOptSolution) = objective_value(sol.model)

function substitute_differentials(model, sys, eqs)
    t = get_iv(sys)
    D = Differential(t)
    diffsubmap = Dict([D(lowered_var(model, :U, i, t)) => lowered_derivative(model, i)
                       for i in 1:length(unknowns(sys))])
    eqs = map(c -> Symbolics.substitute(c, diffsubmap), eqs)
end

function substitute_toterm(vars, exprs)
    toterm_map = Dict([u => default_toterm(value(u)) for u in vars])
    exprs = map(c -> Symbolics.fast_substitute(c, toterm_map), exprs)
end

function substitute_params(pmap, exprs)
    exprs = map(c -> Symbolics.fixpoint_sub(c, Dict(pmap)), exprs)
end

function check_constraint_vars(vars)
    for u in vars
        x = operation(u)
        t = only(arguments(u))
        if (symbolic_type(t) === NotSymbolic())
            error("Provided specific time constraint in a free final time problem. This is not supported by the collocation solvers at the moment. The offending variable is $u. Specific-time user constraints can only be specified at the end of the timespan.")
        end
    end
end

########################
### SOLVER UTILITIES ###
########################
"""
Add the solve constraints, set the solver (Ipopt, e.g.) and solver options, optimize the model.
"""
function prepare_and_optimize! end
function get_t_values end
function get_U_values end
function get_V_values end
function successful_solve end

"""
    solve(prob::AbstractDynamicOptProblem, solver::AbstractCollocation; verbose = false, kwargs...)

- kwargs are used for other options. For example, the `plugin_options` and `solver_options` will propagated to the Opti object in CasADi.
"""
function DiffEqBase.solve(prob::AbstractDynamicOptProblem,
        solver::AbstractCollocation; verbose = false, kwargs...)
    solved_model = prepare_and_optimize!(prob, solver; verbose, kwargs...)

    ts = get_t_values(solved_model)
    Us = get_U_values(solved_model)
    Vs = get_V_values(solved_model)
    is_free_final(prob.wrapped_model) && (ts .+ prob.tspan[1])

    ode_sol = DiffEqBase.build_solution(prob, solver, ts, Us)
    input_sol = isnothing(Vs) ? nothing : DiffEqBase.build_solution(prob, solver, ts, Vs)

    if !successful_solve(solved_model)
        ode_sol = SciMLBase.solution_new_retcode(
            ode_sol, SciMLBase.ReturnCode.ConvergenceFailure)
        !isnothing(input_sol) && (input_sol = SciMLBase.solution_new_retcode(
            input_sol, SciMLBase.ReturnCode.ConvergenceFailure))
    end
    DynamicOptSolution(solved_model, ode_sol, input_sol)
end
