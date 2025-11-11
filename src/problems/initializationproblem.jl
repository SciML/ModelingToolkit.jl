struct InitializationProblem{iip, specialization} end

@doc """
    InitializationProblem(sys::AbstractSystem, t, op = Dict(); kwargs...)
    InitializationProblem{iip}(sys::AbstractSystem, t, op = Dict(); kwargs...)
    InitializationProblem{iip, specialize}(sys::AbstractSystem, t, op = Dict(); kwargs...)

Generate a `NonlinearProblem`, `SCCNonlinearProblem` or `NonlinearLeastSquaresProblem` to
represent a consistent initialization of `sys` given the initial time `t` and operating
point `op`. The initial time can be `nothing` for time-independent systems.

# Keyword arguments

$INITIALIZEPROB_KWARGS
$INTERNAL_INITIALIZEPROB_KWARGS

All other keyword arguments are forwarded to the wrapped nonlinear problem constructor.
""" InitializationProblem

@fallback_iip_specialize function InitializationProblem{iip, specialize}(
        sys::AbstractSystem,
        t, op = Dict();
        fast_path = false,
        guesses = [],
        check_length = true,
        warn_initialize_determined = true,
        initialization_eqs = [],
        fully_determined = nothing,
        check_units = true,
        use_scc = true,
        allow_incomplete = false,
        algebraic_only = false,
        time_dependent_init = is_time_dependent(sys),
        kwargs...) where {iip, specialize}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `mtkcompile` on the system before creating an `ODEProblem`")
    end
    has_u0_ics = false
    for k in keys(op)
        has_u0_ics |= is_variable(sys, k) || isdifferential(k)
    end
    if !has_u0_ics && get_initializesystem(sys) !== nothing
        isys = get_initializesystem(sys; initialization_eqs, check_units)
        simplify_system = false
    elseif !has_u0_ics && get_initializesystem(sys) === nothing
        isys = generate_initializesystem(
            sys; initialization_eqs, check_units, op, guesses, algebraic_only,
            fast_path)
        simplify_system = true
    else
        isys = generate_initializesystem(
            sys; op, initialization_eqs, check_units, time_dependent_init,
            guesses, algebraic_only, fast_path)
        simplify_system = true
    end

    # useful for `SteadyStateProblem` since `f` has to be autonomous and the
    # initialization should be too
    if !time_dependent_init
        idx = findfirst(isequal(get_iv(sys)), get_ps(isys))
        idx === nothing || deleteat!(get_ps(isys), idx)
    end

    if simplify_system
        isys = mtkcompile(isys; fully_determined, split = is_split(sys))
    end

    ts = get_tearing_state(isys)
    unassigned_vars = StructuralTransformations.singular_check(ts)
    if warn_initialize_determined && !isempty(unassigned_vars)
        errmsg = """
        The initialization system is structurally singular. Guess values may \
        significantly affect the initial values of the ODE. The problematic variables \
        are $unassigned_vars.

        Note that the identification of problematic variables is a best-effort heuristic.
        """
        @warn errmsg
    end

    uninit = as_atomic_array_set(unknowns(sys))
    setdiff!(uninit, as_atomic_array_set(unknowns(isys)))
    setdiff!(uninit, as_atomic_array_set(observables(isys)))

    if time_dependent_init && !isempty(uninit)
        allow_incomplete || throw(IncompleteInitializationError(uninit, sys))
        # for incomplete initialization, we will add the missing variables as parameters.
        # they will be updated by `update_initializeprob!` and `initializeprobmap` will
        # use them to construct the new `u0`.
        new_ps = copy(get_ps(isys))
        append!(new_ps, uninit)
        @set! isys.ps = new_ps
        isys = complete(isys)
    end

    neqs = length(equations(isys))
    nunknown = length(unknowns(isys))

    if use_scc
        scc_message = "`SCCNonlinearProblem` can only be used for initialization of fully determined systems and hence will not be used here. "
    else
        scc_message = ""
    end

    if warn_initialize_determined && neqs > nunknown
        @warn "Initialization system is overdetermined. $neqs equations for $nunknown unknowns. Initialization will default to using least squares. $(scc_message)To suppress this warning pass warn_initialize_determined = false. To make this warning into an error, pass fully_determined = true"
    end
    if warn_initialize_determined && neqs < nunknown
        @warn "Initialization system is underdetermined. $neqs equations for $nunknown unknowns. Initialization will default to using least squares. $(scc_message)To suppress this warning pass warn_initialize_determined = false. To make this warning into an error, pass fully_determined = true"
    end

    if t !== nothing
        op[get_iv(sys)] = t
    end
    filter!(!Base.Fix2(===, COMMON_MISSING) ∘ last, op)
    TProb = if neqs == nunknown && isempty(unassigned_vars)
        if use_scc && neqs > 0
            if is_split(isys)
                SCCNonlinearProblem
            else
                @warn "`SCCNonlinearProblem` can only be used with `split = true` systems. Simplify your `ODESystem` with `split = true` or pass `use_scc = false` to disable this warning"
                NonlinearProblem
            end
        else
            NonlinearProblem
        end
    else
        NonlinearLeastSquaresProblem
    end
    TProb{iip}(isys, op; kwargs..., build_initializeprob = false, is_initializeprob = true)
end

const INCOMPLETE_INITIALIZATION_MESSAGE = """
                                Initialization incomplete. Not all of the state variables of the
                                DAE system can be determined by the initialization. Missing
                                variables:
                                """

struct IncompleteInitializationError <: Exception
    uninit::Any
    sys::Any
end

function Base.showerror(io::IO, e::IncompleteInitializationError)
    println(io, INCOMPLETE_INITIALIZATION_MESSAGE)
    println(io, underscore_to_D(collect(e.uninit), e.sys))
end
