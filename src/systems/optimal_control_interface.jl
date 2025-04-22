abstract type AbstractOptimalControlProblem{uType, tType, isinplace} <:
              SciMLBase.AbstractODEProblem{uType, tType, isinplace} end

struct OptimalControlSolution
    model::Any
    sol::ODESolution
    input_sol::Union{Nothing, ODESolution}
end

function Base.show(io::IO, sol::OptimalControlSolution)
    println("retcode: ", sol.sol.retcode, "\n")

    println("Optimal control solution for following model:\n")
    show(sol.model)

    print("\n\nPlease query the model using sol.model, the solution trajectory for the system using sol.sol, or the solution trajectory for the controllers using sol.input_sol.")
end

function JuMPControlProblem end
function InfiniteOptControlProblem end
function CasADiControlProblem end
function PyomoControlProblem end

function warn_overdetermined(sys, u0map)
    constraintsys = get_constraintsystem(sys)
    if !isnothing(constraintsys)
        (length(constraints(constraintsys)) + length(u0map) > length(unknowns(sys))) &&
            @warn "The control problem is overdetermined. The total number of conditions (# constraints + # fixed initial values given by u0map) exceeds the total number of states. The solvers will default to doing a nonlinear least-squares optimization."
    end
end

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
    (f), _, _ = generate_control_function(
        sys, inputs, disturbance_inputs; eval_module, cse, kwargs...)

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

"""
Integral operator. When applied to an expression in a cost
function, assumes that the integration variable is the
iv of the system, and assumes that the bounds are the
tspan.
Equivalent to Integral(t in tspan) in Symbolics.
"""
struct ∫ <: Symbolics.Operator end
∫(x) = ∫()(x)
Base.show(io::IO, x::∫) = print(io, "∫")
Base.nameof(::∫) = :∫

function (I::∫)(x)
    Term{symtype(x)}(I, Any[x])
end

function (I::∫)(x::Num)
    v = value(x)
    Num(I(v))
end

SymbolicUtils.promote_symtype(::Int, t) = t

# returns the JuMP timespan, the number of steps, and whether it is a free time problem.
function process_tspan(tspan, dt, steps)
    is_free_time = false
    if isnothing(dt) && isnothing(steps)
        error("Must provide either the dt or the number of intervals to the collocation solvers (JuMP, InfiniteOpt, CasADi).")
    elseif symbolic_type(tspan[1]) === ScalarSymbolic() ||
           symbolic_type(tspan[2]) === ScalarSymbolic()
        isnothing(steps) &&
            error("Free final time problems require specifying the number of steps, rather than dt.")
        isnothing(dt) ||
            @warn "Specified dt for free final time problem. This will be ignored; dt will be determined by the number of timesteps."

        return steps, true
    else
        isnothing(steps) ||
            @warn "Specified number of steps for problem with concrete tspan. This will be ignored; number of steps will be determined by dt."

        return length(tspan[1]:dt:tspan[2]), false
    end
end

#"""
#$(SIGNATURES)
#
#Define one or more inputs.
#
#See also [`@independent_variables`](@ref), [`@variables`](@ref) and [`@constants`](@ref).
#"""
#macro inputs(xs...)
#    Symbolics._parse_vars(:inputs,
#        Real,
#        xs,
#        toparam) |> esc
#end
