@fallback_iip_specialize function SciMLBase.DiscreteFunction{iip, spec}(
        sys::System, _d = nothing, u0 = nothing, p = nothing;
        t = nothing, eval_expression = false, eval_module = @__MODULE__,
        checkbounds = false, analytic = nothing, simplify = false, cse = true,
        initialization_data = nothing, check_compatibility = true, kwargs...) where {
        iip, spec}
    check_complete(sys, DiscreteFunction)
    check_compatibility && check_compatible_system(DiscreteFunction, sys)

    dvs = unknowns(sys)
    ps = parameters(sys)
    f = generate_rhs(sys, dvs, ps; expression = Val{false},
        eval_expression, eval_module, checkbounds = checkbounds, cse,
        kwargs...)

    if spec === SciMLBase.FunctionWrapperSpecialize && iip
        if u0 === nothing || p === nothing || t === nothing
            error("u0, p, and t must be specified for FunctionWrapperSpecialize on DiscreteFunction.")
        end
        f = SciMLBase.wrapfun_iip(f, (u0, u0, p, t))
    end

    observedfun = ObservedFunctionCache(
        sys; steady_state = false, eval_expression, eval_module, checkbounds, cse)

    DiscreteFunction{iip, spec}(f;
        sys = sys,
        observed = observedfun,
        analytic = analytic,
        initialization_data)
end

@fallback_iip_specialize function SciMLBase.DiscreteProblem{iip, spec}(
        sys::System, u0map, tspan, parammap = SciMLBase.NullParameters();
        check_compatibility = true, kwargs...) where {iip, spec}
    check_complete(sys, DiscreteProblem)
    check_compatibility && check_compatible_system(DiscreteProblem, sys)

    dvs = unknowns(sys)
    u0map = to_varmap(u0map, dvs)
    add_toterms!(u0map; replace = true)
    f, u0, p = process_SciMLProblem(DiscreteFunction{iip, spec}, sys, u0map, parammap;
        t = tspan !== nothing ? tspan[1] : tspan, check_compatibility, kwargs...)
    u0 = f(u0, p, tspan[1])

    kwargs = process_kwargs(sys; kwargs...)
    # Call `remake` so it runs initialization if it is trivial
    return remake(DiscreteProblem{iip}(f, u0, tspan, p; kwargs...))
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
