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

function warn_overdetermined(sys, u0map)
    constraintsys = get_constraintsystem(sys)
    if !isnothing(constraintsys)
        (length(constraints(constraintsys)) + length(u0map) > length(unknowns(sys))) &&
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

"""
Generate the control function f(x, u, p, t) from the ODESystem. 
Input variables are automatically inferred but can be manually specified.
"""
function SciMLBase.ODEInputFunction{iip, specialize}(sys::ODESystem,
        dvs = unknowns(sys),
        ps = parameters(sys), u0 = nothing,
        inputs = unbound_inputs(sys),
        disturbance_inputs = disturbances(sys);
        version = nothing, tgrad = false,
        jac = false, controljac = false,
        p = nothing, t = nothing,
        eval_expression = false,
        sparse = false, simplify = false,
        eval_module = @__MODULE__,
        steady_state = false,
        checkbounds = false,
        sparsity = false,
        analytic = nothing,
        split_idxs = nothing,
        initialization_data = nothing,
        cse = true,
        kwargs...) where {iip, specialize}
    f, _, _ = generate_control_function(
        sys, inputs, disturbance_inputs; eval_module, cse, kwargs...)
    f = f[1]

    if tgrad
        tgrad_gen = generate_tgrad(sys, dvs, ps;
            simplify = simplify,
            expression = Val{true},
            expression_module = eval_module, cse,
            checkbounds = checkbounds, kwargs...)
        tgrad_oop, tgrad_iip = eval_or_rgf.(tgrad_gen; eval_expression, eval_module)
        _tgrad = GeneratedFunctionWrapper{(2, 3, is_split(sys))}(tgrad_oop, tgrad_iip)
    else
        _tgrad = nothing
    end

    if jac
        jac_gen = generate_jacobian(sys, dvs, ps;
            simplify = simplify, sparse = sparse,
            expression = Val{true},
            expression_module = eval_module, cse,
            checkbounds = checkbounds, kwargs...)
        jac_oop, jac_iip = eval_or_rgf.(jac_gen; eval_expression, eval_module)

        _jac = GeneratedFunctionWrapper{(2, 3, is_split(sys))}(jac_oop, jac_iip)
    else
        _jac = nothing
    end

    if controljac
        cjac_gen = generate_control_jacobian(sys, dvs, ps;
            simplify = simplify, sparse = sparse,
            expression = Val{true},
            expression_module = eval_module, cse,
            checkbounds = checkbounds, kwargs...)
        cjac_oop, cjac_iip = eval_or_rgf.(cjac_gen; eval_expression, eval_module)

        _cjac = GeneratedFunctionWrapper{(2, 3, is_split(sys))}(cjac_oop, cjac_iip)
    else
        _cjac = nothing
    end

    M = calculate_massmatrix(sys)
    _M = if sparse && !(u0 === nothing || M === I)
        SparseArrays.sparse(M)
    elseif u0 === nothing || M === I
        M
    else
        ArrayInterface.restructure(u0 .* u0', M)
    end

    observedfun = ObservedFunctionCache(
        sys; steady_state, eval_expression, eval_module, checkbounds, cse)

    if sparse
        uElType = u0 === nothing ? Float64 : eltype(u0)
        W_prototype = similar(W_sparsity(sys), uElType)
        controljac_prototype = similar(calculate_control_jacobian(sys), uElType)
    else
        W_prototype = nothing
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
        sparsity = sparsity ? W_sparsity(sys) : nothing,
        analytic = analytic,
        initialization_data)
end

function SciMLBase.ODEInputFunction(sys::AbstractODESystem, args...; kwargs...)
    ODEInputFunction{true}(sys, args...; kwargs...)
end

function SciMLBase.ODEInputFunction{true}(sys::AbstractODESystem, args...;
        kwargs...)
    ODEInputFunction{true, SciMLBase.AutoSpecialize}(sys, args...; kwargs...)
end

function SciMLBase.ODEInputFunction{false}(sys::AbstractODESystem, args...;
        kwargs...)
    ODEInputFunction{false, SciMLBase.FullSpecialize}(sys, args...; kwargs...)
end

# returns the JuMP timespan, the number of steps, and whether it is a free time problem.
function process_tspan(tspan, dt, steps)
    is_free_time = false
    if isnothing(dt) && isnothing(steps)
        error("Must provide either the dt or the number of intervals to the collocation solvers (JuMP, InfiniteOpt, CasADi).")
    elseif symbolic_type(tspan[1]) === ScalarSymbolic() ||
           symbolic_type(tspan[2]) === ScalarSymbolic()
        isnothing(steps) &&
            error("Free final time problems require specifying the number of steps using the keyword arg `steps`, rather than dt.")
        isnothing(dt) ||
            @warn "Specified dt for free final time problem. This will be ignored; dt will be determined by the number of timesteps."

        return steps, true
    else
        isnothing(steps) ||
            @warn "Specified number of steps for problem with concrete tspan. This will be ignored; number of steps will be determined by dt."

        return length(tspan[1]:dt:tspan[2]), false
    end
end

function process_DynamicOptProblem(prob_type::AbstractDynamicOptProblem, model_type, sys::ODESystem, u0map, tspan, pmap;
    dt = nothing,
    steps = nothing,
    guesses = Dict(), kwargs...)

    MTK.warn_overdetermined(sys, u0map)
    _u0map = has_alg_eqs(sys) ? u0map : merge(Dict(u0map), Dict(guesses))
    f, u0, p = MTK.process_SciMLProblem(ODEInputFunction, sys, _u0map, pmap;
        t = tspan !== nothing ? tspan[1] : tspan, kwargs...)

    stidxmap = Dict([v => i for (i, v) in enumerate(states)])
    u0map = Dict([MTK.default_toterm(MTK.value(k)) => v for (k, v) in u0map])
    u0_idxs = has_alg_eqs(sys) ? collect(1:length(states)) :
              [stidxmap[MTK.default_toterm(k)] for (k, v) in u0map]
    pmap = Dict{Any, Any}(pmap)
    steps, is_free_t = MTK.process_tspan(tspan, dt, steps)

    ctrls = MTK.unbound_inputs(sys)
    states = unknowns(sys)

    model = generate_internal_model(model_type)
    U = generate_U(model, u0)
    V = generate_V()
    tₛ = generate_timescale()
    fullmodel = model_type(model, U, V, tₛ)

    set_variable_bounds!(fullmodel, sys, pmap)
    add_cost_function!(fullmodel, sys, tspan, pmap; is_free_t)
    add_user_constraints!(fullmodel, sys, tspan, pmap; is_free_t)
    add_initial_constraints!(fullmodel, u0, u0_idxs)

    prob_type(f, u0, tspan, p, fullmodel, kwargs...)
end

function add_cost_function!()
    jcosts = copy(MTK.get_costs(sys))
    consolidate = MTK.get_consolidate(sys)
    if isnothing(jcosts) || isempty(jcosts)
        minimize!(opti, MX(0))
        return
    end

    jcosts = substitute_model_vars(model, sys, pmap, jcosts; is_free_t)
    jcosts = substitute_free_final_vars(model, sys, pmap, jcosts; is_free_t)
    jcosts = substitute_fixed_t_vars(model, sys, pmap, jcosts; is_free_t)
    jcosts = substitute_integral()
end

function add_user_constraints!()
    conssys = MTK.get_constraintsystem(sys)
    jconstraints = isnothing(conssys) ? nothing : MTK.get_constraints(conssys)
    (isnothing(jconstraints) || isempty(jconstraints)) && return nothing

    auxmap = Dict([u => MTK.default_toterm(MTK.value(u)) for u in unknowns(conssys)])
    jconstraints = substitute_model_vars(model, sys, pmap, jconstraints; auxmap, is_free_t)

    for c in jconstraints
        if cons isa Equation
            add_constraint!()
        elseif cons.relational_op === Symbolics.geq
            add_constraint!()
        else
            add_constraint!()
        end
    end
end

function generate_U end
function generate_V end
function generate_timescale end

function add_initial_constraints! end
function add_constraint! end

function add_collocation_solve_constraints!(prob, tableau)
    nᵤ = size(U.u, 1)
    nᵥ = size(V.u, 1)

    if is_explicit(tableau)
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
            # Kᵢ = generate_K()
            Kᵢ = variable!(solver_opti, nᵤ, length(α))
            ΔUs = A * Kᵢ' # the stepsize at each stage of the implicit method
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
end

function add_equational_solve_constraints!()
    diff_eqs = substitute_differentials()
    add_constraint!()

    alg_eqs = substitute_model_vars()
    add_constraint!()
end

"""
Add the solve constraints, set the solver (Ipopt, e.g.)
"""
function prepare_solver end

function DiffEqBase.solve(prob::AbstractDynamicOptProblem, solver::AbstractCollocation)
    #add_solve_constraints!(prob, solver)
    solver = prepare_solver(prob, solver)
    sol = solve_prob(prob, solver)
    
    ts = get_t_values(sol)
    Us = get_U_values(sol)
    Vs = get_V_values(sol)

    ode_sol = DiffEqBase.build_solution(prob, solver, ts, Us)
    input_sol = DiffEqBase.build_solution(prob, solver, ts, Vs)

    if successful_solve(model) 
        ode_sol = SciMLBase.solution_new_retcode(
            ode_sol, SciMLBase.ReturnCode.ConvergenceFailure)
        !isnothing(input_sol) && (input_sol = SciMLBase.solution_new_retcode(
            input_sol, SciMLBase.ReturnCode.ConvergenceFailure))
    end
    DynamicOptSolution(model, ode_sol, input_sol)
end
