function SciMLBase.OptimizationFunction(sys::System, args...; kwargs...)
    return OptimizationFunction{true}(sys, args...; kwargs...)
end

function SciMLBase.OptimizationFunction{iip}(sys::System;
        u0 = nothing, p = nothing, grad = false, hess = false,
        sparse = false, cons_j = false, cons_h = false, cons_sparse = false,
        linenumbers = true, eval_expression = false, eval_module = @__MODULE__,
        simplify = false, check_compatibility = true, checkbounds = false, cse = true,
        expression = Val{false}, kwargs...) where {iip}
    check_complete(sys, OptimizationFunction)
    check_compatibility && check_compatible_system(OptimizationFunction, sys)

    cstr = constraints(sys)

    f = generate_cost(sys; expression, wrap_gfw = Val{true}, eval_expression,
        eval_module, checkbounds, cse, kwargs...)

    if grad
        _grad = generate_cost_gradient(sys; expression, wrap_gfw = Val{true},
            eval_expression, eval_module, checkbounds, cse, kwargs...)
    else
        _grad = nothing
    end
    if hess
        _hess,
        hess_prototype = generate_cost_hessian(
            sys; expression, wrap_gfw = Val{true}, eval_expression,
            eval_module, checkbounds, cse, sparse, simplify, return_sparsity = true,
            kwargs...)
    else
        _hess = hess_prototype = nothing
        if sparse
            hess_prototype = cost_hessian_sparsity(sys)
        end
    end
    if isempty(cstr)
        cons = _cons_j = cons_jac_prototype = _cons_h = nothing
        cons_hess_prototype = cons_expr = nothing
    else
        cons = generate_cons(sys; expression, wrap_gfw = Val{true},
            eval_expression, eval_module, checkbounds, cse, kwargs...)
        if cons_j
            _cons_j,
            cons_jac_prototype = generate_constraint_jacobian(
                sys; expression, wrap_gfw = Val{true}, eval_expression,
                eval_module, checkbounds, cse, simplify, sparse = cons_sparse,
                return_sparsity = true, kwargs...)
        else
            _cons_j = cons_jac_prototype = nothing
        end
        if cons_h
            _cons_h,
            cons_hess_prototype = generate_constraint_hessian(
                sys; expression, wrap_gfw = Val{true}, eval_expression,
                eval_module, checkbounds, cse, simplify, sparse = cons_sparse,
                return_sparsity = true, kwargs...)
        else
            _cons_h = cons_hess_prototype = nothing
        end
        cons_expr = Code.toexpr.(expand.([eq.lhs for eq in Symbolics.canonical_form.(cstr)]))
    end

    obj_expr = Code.toexpr(expand(cost(sys)))

    observedfun = ObservedFunctionCache(
        sys; expression, eval_expression, eval_module, checkbounds, cse)

    args = (; f, ad = SciMLBase.NoAD())
    kwargs = (;
        sys = sys,
        grad = _grad,
        hess = _hess,
        hess_prototype = hess_prototype,
        cons = cons,
        cons_j = _cons_j,
        cons_jac_prototype = cons_jac_prototype,
        cons_h = _cons_h,
        cons_hess_prototype = cons_hess_prototype,
        cons_expr = cons_expr,
        expr = obj_expr,
        observed = observedfun)

    return maybe_codegen_scimlfn(expression, OptimizationFunction{iip}, args; kwargs...)
end

function SciMLBase.OptimizationProblem(sys::System, args...; kwargs...)
    return OptimizationProblem{true}(sys, args...; kwargs...)
end

function SciMLBase.OptimizationProblem{iip}(
        sys::System, op; lb = nothing,
        ub = nothing, check_compatibility = true, expression = Val{false},
        kwargs...) where {iip}
    check_complete(sys, OptimizationProblem)
    check_compatibility && check_compatible_system(OptimizationProblem, sys)

    f, u0,
    p = process_SciMLProblem(OptimizationFunction{iip}, sys, op;
        check_compatibility, tofloat = false, check_length = false, expression, kwargs...)

    dvs = unknowns(sys)
    int = symtype.(unwrap.(dvs)) .<: Integer
    if lb === nothing && ub === nothing
        lb = first.(getbounds.(dvs))
        ub = last.(getbounds.(dvs))
        isboolean = symtype.(unwrap.(dvs)) .<: Bool
        lb[isboolean] .= 0
        ub[isboolean] .= 1
    else
        xor(isnothing(lb), isnothing(ub)) &&
            throw(ArgumentError("Expected both `lb` and `ub` to be supplied"))
        !isnothing(lb) && length(lb) != length(dvs) &&
            throw(ArgumentError("Expected both `lb` to be of the same length as the vector of optimization variables"))
        !isnothing(ub) && length(ub) != length(dvs) &&
            throw(ArgumentError("Expected both `ub` to be of the same length as the vector of optimization variables"))
    end

    ps = parameters(sys)
    defs = defaults(sys)
    op = to_varmap(op, dvs)
    lbmap = merge(op, AnyDict(dvs .=> lb))
    _, _ = build_operating_point!(sys, lbmap, Dict(), Dict(), defs, dvs, ps)
    lb = varmap_to_vars(lbmap, dvs; tofloat = false)
    ubmap = merge(op, AnyDict(dvs .=> ub))
    _, _ = build_operating_point!(sys, ubmap, Dict(), Dict(), defs, dvs, ps)
    ub = varmap_to_vars(ubmap, dvs; tofloat = false)

    if !isnothing(lb) && all(lb .== -Inf) && !isnothing(ub) && all(ub .== Inf)
        lb = nothing
        ub = nothing
    end

    cstr = constraints(sys)
    if isempty(cstr)
        lcons = ucons = nothing
    else
        lcons = fill(-Inf, length(cstr))
        ucons = zeros(length(cstr))
        lcons[findall(Base.Fix2(isa, Equation), cstr)] .= 0.0
    end

    kwargs = process_kwargs(sys; kwargs...)
    kwargs = (; lb, ub, int, lcons, ucons, kwargs...)
    args = (; f, u0, p)
    return maybe_codegen_scimlproblem(expression, OptimizationProblem{iip}, args; kwargs...)
end

function check_compatible_system(
        T::Union{Type{OptimizationFunction}, Type{OptimizationProblem}}, sys::System)
    check_time_independent(sys, T)
    check_not_dde(sys)
    check_has_cost(sys, T)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
    check_no_equations(sys, T)
end
