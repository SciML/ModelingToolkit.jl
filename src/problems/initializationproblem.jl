struct InitializationProblem{iip, specialization} end

"""
```julia
InitializationProblem{iip}(sys::AbstractSystem, t, u0map,
                           parammap = DiffEqBase.NullParameters();
                           version = nothing, tgrad = false,
                           jac = false,
                           checkbounds = false, sparse = false,
                           simplify = false,
                           linenumbers = true, parallel = SerialForm(),
                           initialization_eqs = [],
                           fully_determined = false,
                           kwargs...) where {iip}
```

Generates a NonlinearProblem or NonlinearLeastSquaresProblem from a System
which represents the initialization, i.e. the calculation of the consistent
initial conditions for the given DAE.
"""
@fallback_iip_specialize function InitializationProblem{iip, specialize}(
        sys::AbstractSystem,
        t, u0map = [],
        parammap = DiffEqBase.NullParameters();
        guesses = [],
        check_length = true,
        warn_initialize_determined = true,
        initialization_eqs = [],
        fully_determined = nothing,
        check_units = true,
        use_scc = true,
        allow_incomplete = false,
        force_time_independent = false,
        algebraic_only = false,
        time_dependent_init = is_time_dependent(sys),
        kwargs...) where {iip, specialize}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system before creating an `ODEProblem`")
    end
    if isempty(u0map) && get_initializesystem(sys) !== nothing
        isys = get_initializesystem(sys; initialization_eqs, check_units)
        simplify_system = false
    elseif isempty(u0map) && get_initializesystem(sys) === nothing
        isys = generate_initializesystem(
            sys; initialization_eqs, check_units, pmap = parammap,
            guesses, algebraic_only)
        simplify_system = true
    else
        isys = generate_initializesystem(
            sys; u0map, initialization_eqs, check_units, time_dependent_init,
            pmap = parammap, guesses, algebraic_only)
        simplify_system = true
    end

    # useful for `SteadyStateProblem` since `f` has to be autonomous and the
    # initialization should be too
    if force_time_independent
        idx = findfirst(isequal(get_iv(sys)), get_ps(isys))
        idx === nothing || deleteat!(get_ps(isys), idx)
    end

    if simplify_system
        isys = structural_simplify(isys; fully_determined)
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

    uninit = setdiff(unknowns(sys), [unknowns(isys); observables(isys)])

    # TODO: throw on uninitialized arrays
    filter!(x -> !(x isa Symbolics.Arr), uninit)
    if time_dependent_init && !isempty(uninit)
        allow_incomplete || throw(IncompleteInitializationError(uninit))
        # for incomplete initialization, we will add the missing variables as parameters.
        # they will be updated by `update_initializeprob!` and `initializeprobmap` will
        # use them to construct the new `u0`.
        newparams = map(toparam, uninit)
        append!(get_ps(isys), newparams)
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

    parammap = recursive_unwrap(anydict(parammap))
    if t !== nothing
        parammap[get_iv(sys)] = t
    end
    filter!(kvp -> kvp[2] !== missing, parammap)

    u0map = to_varmap(u0map, unknowns(sys))
    if isempty(guesses)
        guesses = Dict()
    end

    filter_missing_values!(u0map)
    filter_missing_values!(parammap)
    u0map = merge(ModelingToolkit.guesses(sys), todict(guesses), u0map)

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
    TProb(isys, u0map, parammap; kwargs...,
        build_initializeprob = false, is_initializeprob = true)
end

const INCOMPLETE_INITIALIZATION_MESSAGE = """
                                Initialization incomplete. Not all of the state variables of the
                                DAE system can be determined by the initialization. Missing
                                variables:
                                """

struct IncompleteInitializationError <: Exception
    uninit::Any
end

function Base.showerror(io::IO, e::IncompleteInitializationError)
    println(io, INCOMPLETE_INITIALIZATION_MESSAGE)
    println(io, e.uninit)
end
