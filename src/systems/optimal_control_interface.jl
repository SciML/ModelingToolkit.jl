abstract type AbstractDynamicOptProblem{uType, tType, isinplace} <:
              SciMLBase.AbstractODEProblem{uType, tType, isinplace} end

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
