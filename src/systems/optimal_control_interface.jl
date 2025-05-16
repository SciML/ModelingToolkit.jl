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

function JuMPDynamicOptProblem end
function InfiniteOptDynamicOptProblem end
function CasADiDynamicOptProblem end

function JuMPCollocation end
function InfiniteOptCollocation end
function CasADiCollocation end

function warn_overdetermined(sys, u0map)
    cstrs = constraints(sys)
    if !isempty(cstrs)
        (length(cstrs) + length(u0map) > length(unknowns(sys))) &&
            @warn "The control problem is overdetermined. The total number of conditions (# constraints + # fixed initial values given by u0map) exceeds the total number of states. The solvers will default to doing a nonlinear least-squares optimization."
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
    f, _, _ = generate_control_function(
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
function process_DynamicOptProblem(prob_type::Type{<:AbstractDynamicOptProblem}, model_type, sys::ODESystem, u0map, tspan, pmap;
    dt = nothing,
    steps = nothing,
    guesses = Dict(), kwargs...)
    warn_overdetermined(sys, u0map)
    ctrls = unbound_inputs(sys)
    states = unknowns(sys)

    _u0map = has_alg_eqs(sys) ? u0map : merge(Dict(u0map), Dict(guesses))
    stidxmap = Dict([v => i for (i, v) in enumerate(states)])
    u0map = Dict([default_toterm(value(k)) => v for (k, v) in u0map])
    u0_idxs = has_alg_eqs(sys) ? collect(1:length(states)) :
              [stidxmap[default_toterm(k)] for (k, v) in u0map]

    f, u0, p = process_SciMLProblem(ODEInputFunction, sys, _u0map, pmap;
        t = tspan !== nothing ? tspan[1] : tspan, kwargs...)
    model_tspan, steps, is_free_t = process_tspan(tspan, dt, steps)

    pmap = recursive_unwrap(AnyDict(pmap))
    evaluate_varmap!(pmap, keys(pmap))
    c0 = value.([pmap[c] for c in ctrls])

    tsteps = LinRange(model_tspan[1], model_tspan[2], steps)
    model = generate_internal_model(model_type)
    generate_time_variable!(model, model_tspan, steps)
    U = generate_state_variable!(model, u0, length(states), length(steps))
    V = generate_input_variable!(model, c0, length(ctrls), length(steps))
    tₛ = generate_timescale!(model, get(pmap, tspan[2], tspan[2]), is_free_t)
    fullmodel = model_type(model, U, V, tₛ, is_free_t)

    set_variable_bounds!(fullmodel, sys, pmap, tspan[2])
    add_cost_function!(fullmodel, sys, tspan, pmap)
    add_user_constraints!(fullmodel, sys, tspan, pmap)
    add_initial_constraints!(fullmodel, u0, u0_idxs, model_tspan[1])

    prob_type(f, u0, tspan, p, fullmodel, kwargs...)
end

function generate_time_variable! end
function generate_internal_model end
function generate_state_variable! end
function generate_input_variable! end
function generate_timescale! end
function set_variable_bounds! end
function add_initial_constraints! end
function add_constraint! end
is_free_final(model) = model.is_free_final 

function add_cost_function!(model, sys, tspan, pmap)
    jcosts = copy(get_costs(sys))
    consolidate = get_consolidate(sys)
    if isnothing(jcosts) || isempty(jcosts)
        set_objective!(model, 0)
        return
    end
    jcosts = substitute_model_vars(model, sys, jcosts; tf = tspan[2])
    jcosts = substitute_params(pmap, jcosts)
    jcosts = substitute_integral(model, jcosts)
    set_objective!(model, consolidate(jcosts))
end

function add_user_constraints!(model, sys, tspan, pmap)
    conssys = get_constraintsystem(sys)
    jconstraints = isnothing(conssys) ? nothing : get_constraints(conssys)
    (isnothing(jconstraints) || isempty(jconstraints)) && return nothing
    consvars = get_unknowns(conssys)
    is_free_final(model) && check_constraint_vars(consvars)

    jconstraints = substitute_toterm(consvars, jconstraints)
    jconstraints = substitute_params(pmap, jconstraints)
    jconstraints = substitute_model_vars(model, sys, jconstraints; tf = tspan[2])

    for c in jconstraints
        @show c
        add_constraint!(model, c)
    end
end

function add_equational_constraints!(model, sys, tspan)
    model = model.model
    diff_eqs = substitute_model_vars(model, sys, diff_equations(sys); tf = tspan[2])
    diff_eqs = substitute_differentials(model, sys, diff_eqs)
    for eq in diff_eqs
        add_constraint!(model, eq.lhs ~ eq.rhs * model.tₛ)
    end

    alg_eqs = substitute_model_vars(model, sys, alg_equations(sys); tf = tspan[2])
    for eq in alg_eqs
        add_constraint!(model, eq.lhs ~ eq.rhs * model.tₛ)
    end
end

function set_objective! end
"""Substitute variables like x(1.5) with the corresponding model variables."""
function substitute_model_vars end
function substitute_integral end
function substitute_differentials end

function substitute_toterm(vars, exprs)
    toterm_map = Dict([u => default_toterm(value(u)) for u in vars])
    exprs = map(c -> Symbolics.fast_substitute(c, toterm_map), exprs)
end

function substitute_params(pmap, exprs)
    exprs = map(c -> Symbolics.fast_substitute(c, Dict(pmap)), exprs)
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
Add the solve constraints, set the solver (Ipopt, e.g.) and solver options.
"""
function prepare_solver! end
function optimize_model! end
function get_t_values end
function get_U_values end
function get_V_values end
function successful_solve end

"""
    solve(prob::AbstractDynamicOptProblem, solver::AbstractCollocation; verbose = false, kwargs...)

- kwargs are used for other options. For example, the `plugin_options` and `solver_options` will propagated to the Opti object in CasADi.
"""
function DiffEqBase.solve(prob::AbstractDynamicOptProblem, solver::AbstractCollocation; verbose = false, kwargs...)
    solver = prepare_solver!(prob, solver; verbose, kwargs...)
    model = optimize_model!(prob, solver)
    
    ts = get_t_values(model)
    Us = get_U_values(model)
    Vs = get_V_values(model)
    is_free_final(model) && (ts .+ tspan[1])

    ode_sol = DiffEqBase.build_solution(prob, solver, ts, Us)
    input_sol = isnothing(Vs) ? nothing : DiffEqBase.build_solution(prob, solver, ts, Vs)

    if !successful_solve(model)
        ode_sol = SciMLBase.solution_new_retcode(
            ode_sol, SciMLBase.ReturnCode.ConvergenceFailure)
        !isnothing(input_sol) && (input_sol = SciMLBase.solution_new_retcode(
            input_sol, SciMLBase.ReturnCode.ConvergenceFailure))
    end
    DynamicOptSolution(model, ode_sol, input_sol)
end
