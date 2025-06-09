@fallback_iip_specialize function SciMLBase.DiscreteFunction{iip, spec}(
        sys::System; u0 = nothing, p = nothing, t = nothing,
        eval_expression = false, eval_module = @__MODULE__, expression = Val{false},
        checkbounds = false, analytic = nothing, simplify = false, cse = true,
        initialization_data = nothing, check_compatibility = true,
        kwargs...) where {iip, spec}
    check_complete(sys, DiscreteFunction)
    check_compatibility && check_compatible_system(DiscreteFunction, sys)

    f = generate_rhs(sys; expression, wrap_gfw = Val{true},
        eval_expression, eval_module, checkbounds = checkbounds, cse,
        kwargs...)

    if spec === SciMLBase.FunctionWrapperSpecialize && iip
        if u0 === nothing || p === nothing || t === nothing
            error("u0, p, and t must be specified for FunctionWrapperSpecialize on DiscreteFunction.")
        end
        if expression == Val{true}
            f = :($(SciMLBase.wrapfun_iip)($f, ($u0, $u0, $p, $t)))
        else
            f = SciMLBase.wrapfun_iip(f, (u0, u0, p, t))
        end
    end

    observedfun = ObservedFunctionCache(
        sys; steady_state = false, expression, eval_expression, eval_module, checkbounds,
        cse)

    kwargs = (;
        sys = sys,
        observed = observedfun,
        analytic = analytic,
        initialization_data)
    args = (; f)

    return maybe_codegen_scimlfn(expression, DiscreteFunction{iip, spec}, args; kwargs...)
end

@fallback_iip_specialize function SciMLBase.DiscreteProblem{iip, spec}(
        sys::System, op, tspan;
        check_compatibility = true, expression = Val{false}, kwargs...) where {iip, spec}
    check_complete(sys, DiscreteProblem)
    check_compatibility && check_compatible_system(DiscreteProblem, sys)

    dvs = unknowns(sys)
    op = to_varmap(op, dvs)
    add_toterms!(op; replace = true)
    f, u0,
    p = process_SciMLProblem(DiscreteFunction{iip, spec}, sys, op;
        t = tspan !== nothing ? tspan[1] : tspan, check_compatibility, expression,
        kwargs...)

    if expression == Val{true}
        u0 = :(f($u0, p, tspan[1]))
    else
        u0 = f(u0, p, tspan[1])
    end

    kwargs = process_kwargs(sys; kwargs...)
    args = (; f, u0, tspan, p)

    return maybe_codegen_scimlproblem(expression, DiscreteProblem{iip}, args; kwargs...)
end

function check_compatible_system(
        T::Union{Type{DiscreteFunction}, Type{DiscreteProblem}}, sys::System)
    check_time_dependent(sys, T)
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
    check_is_discrete(sys, T)
    check_is_explicit(sys, T, ImplicitDiscreteProblem)
end

function shift_u0map_forward(sys::System, u0map, defs)
    iv = get_iv(sys)
    updated = AnyDict()
    for k in collect(keys(u0map))
        v = u0map[k]
        if !((op = operation(k)) isa Shift)
            isnothing(getunshifted(k)) &&
                error("Initial conditions must be for the past state of the unknowns. Instead of providing the condition for $k, provide the condition for $(Shift(iv, -1)(k)).")

            updated[Shift(iv, 1)(k)] = v
        elseif op.steps > 0
            error("Initial conditions must be for the past state of the unknowns. Instead of providing the condition for $k, provide the condition for $(Shift(iv, -1)(only(arguments(k)))).")
        else
            updated[Shift(iv, op.steps + 1)(only(arguments(k)))] = v
        end
    end
    for var in unknowns(sys)
        op = operation(var)
        root = getunshifted(var)
        shift = getshift(var)
        isnothing(root) && continue
        (haskey(updated, Shift(iv, shift)(root)) || haskey(updated, var)) && continue
        haskey(defs, root) || error("Initial condition for $var not provided.")
        updated[var] = defs[root]
    end
    return updated
end
