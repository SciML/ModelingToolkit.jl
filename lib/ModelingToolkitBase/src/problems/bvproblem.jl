@fallback_iip_specialize function SciMLBase.BVProblem{iip, spec}(
        sys::System, op, tspan;
        check_compatibility = true, cse = true,
        checkbounds = false, eval_expression = false, eval_module = @__MODULE__,
        expression = Val{false}, guesses = Dict(), callback = nothing,
        kwargs...
    ) where {iip, spec}
    check_complete(sys, BVProblem)
    check_compatibility && check_compatible_system(BVProblem, sys)
    isnothing(callback) || error("BVP solvers do not support callbacks.")

    dvs = unknowns(sys)
    op = to_varmap(op, dvs)
    # Systems without algebraic equations should use both fixed values + guesses
    # for initialization.
    _op = has_alg_eqs(sys) ? op : merge(Dict(op), Dict(guesses))

    fode, u0,
        p = process_SciMLProblem(
        ODEFunction{iip, spec}, sys, _op; guesses,
        t = tspan !== nothing ? tspan[1] : tspan, check_compatibility = false, cse,
        checkbounds, time_dependent_init = false, expression, kwargs...
    )

    fcost = generate_bvp_cost(sys; expression = Val{false}, wrap_gfw = Val{false},
        eval_expression, eval_module, cse, checkbounds)

    stidxmap = Dict([v => i for (i, v) in enumerate(dvs)])
    u0_idxs = has_alg_eqs(sys) ? collect(1:length(dvs)) :
        [stidxmap[k] for (k, v) in op if haskey(stidxmap, k)]
    fbc = generate_boundary_conditions(
        sys, u0, u0_idxs, tspan[1]; expression = Val{false},
        wrap_gfw = Val{true}, cse, checkbounds
    )

    n_controls = length(unbound_inputs(sys))
    f_prototype = n_controls > 0 ? zeros(eltype(u0), length(dvs) - n_controls) : nothing
    bcresid_prototype = zeros(eltype(u0), length(u0_idxs) + length(constraints(sys)))

    bvpfn = BVPFunction{iip}(fode, fbc; cost = fcost, f_prototype, bcresid_prototype)

    if (length(constraints(sys)) + length(op) > length(dvs))
        @warn "The BVProblem is overdetermined. The total number of conditions (# constraints + # fixed initial values given by op) exceeds the total number of states. The BVP solvers will default to doing a nonlinear least-squares optimization."
    end

    kwargs = process_kwargs(sys; expression, kwargs...)
    args = (; bvpfn, u0, tspan, p)

    return maybe_codegen_scimlproblem(expression, BVProblem{iip}, args; kwargs...)
end

function check_compatible_system(T::Type{BVProblem}, sys::System)
    check_time_dependent(sys, T)
    check_not_dde(sys)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
    check_is_continuous(sys, T)

    return if !isempty(discrete_events(sys)) || !isempty(continuous_events(sys))
        throw(SystemCompatibilityError("BVP solvers do not support events."))
    end
end
