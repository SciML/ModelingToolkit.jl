@fallback_iip_specialize function SciMLBase.BVProblem{iip, spec}(
        sys::System, u0map, tspan, parammap = SciMLBase.NullParameters();
        check_compatibility = true, cse = true, checkbounds = false, eval_expression = false,
        eval_module = @__MODULE__, guesses = Dict(), kwargs...) where {iip, spec}
    check_complete(sys, BVProblem)
    check_compatibility && check_compatible_system(BVProblem, sys)

    # ODESystems without algebraic equations should use both fixed values + guesses
    # for initialization.
    _u0map = has_alg_eqs(sys) ? u0map : merge(Dict(u0map), Dict(guesses))
    fode, u0, p = process_SciMLProblem(
        ODEFunction{iip, spec}, sys, _u0map, parammap; guesses,
        t = tspan !== nothing ? tspan[1] : tspan, check_compatibility = false, cse, checkbounds,
        time_dependent_init = false, kwargs...)

    dvs = unknowns(sys)
    stidxmap = Dict([v => i for (i, v) in enumerate(dvs)])
    u0_idxs = has_alg_eqs(sys) ? collect(1:length(dvs)) : [stidxmap[k] for (k, v) in u0map]
    fbc = generate_boundary_conditions(
        sys, u0, u0_idxs, tspan; expression = Val{false}, cse, checkbounds)
    kwargs = process_kwargs(sys; kwargs...)
    # Call `remake` so it runs initialization if it is trivial
    return remake(BVProblem{iip}(fode, fbc, u0, tspan[1], p; kwargs...))
end

function check_compatible_system(T::Union{Type{BVPFunction}, Type{BVProblem}}, sys::System)
    check_time_dependent(sys, T)
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
    check_is_continuous(sys, T)
end
