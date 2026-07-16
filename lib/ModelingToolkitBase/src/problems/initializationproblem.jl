struct InitializationProblem{iip, specialization} end

@doc """
    InitializationProblem(sys::AbstractSystem, t, op = Dict(); kwargs...)
    InitializationProblem{iip}(sys::AbstractSystem, t, op = Dict(); kwargs...)
    InitializationProblem{iip, specialize}(sys::AbstractSystem, t, op = Dict(); kwargs...)

Generate a `LinearProblem`, `NonlinearProblem`, `SCCNonlinearProblem`,
`NonlinearLeastSquaresProblem` or `SciMLBase.HomotopyProblem` to represent a consistent
initialization of `sys` given the initial time `t` and operating point `op`. The initial
time can be `nothing` for time-independent systems. A `LinearProblem` is used when the
initialization system is linear (affine). A `SciMLBase.HomotopyProblem` is used when the
(square) initialization system contains Modelica `homotopy(actual, simplified)` nodes, so
the initialization is solved by continuation from the `simplified` form (see
[`homotopy`](@ref)).

# Keyword arguments

$INITIALIZEPROB_KWARGS
$INTERNAL_INITIALIZEPROB_KWARGS

All other keyword arguments are forwarded to the wrapped problem constructor.
""" InitializationProblem

@fallback_iip_specialize function InitializationProblem{iip, specialize}(
        sys::AbstractSystem,
        t, op = Dict();
        fast_path = false,
        guesses = [],
        # `check_length` defaults to `nothing`, meaning "no opinion — let the underlying
        # problem constructor apply its own default". It is only forwarded below when the
        # caller explicitly sets it. This matters because the initialization problem types
        # have different defaults (`NonlinearProblem` uses `true` for square systems,
        # `NonlinearLeastSquaresProblem` uses `false` for non-square); forwarding a single
        # value unconditionally would override and break one of those cases.
        check_length = nothing,
        warn_initialize_determined = true,
        initialization_eqs = [],
        fully_determined = nothing,
        check_units = true,
        allow_incomplete = false,
        algebraic_only = false,
        time_dependent_init = is_time_dependent(sys),
        initsys_mtkcompile_kwargs = (;),
        is_steadystateprob = false,
        kwargs...
    ) where {iip, specialize}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `mtkcompile` on the system before creating an `ODEProblem`")
    end
    _iip = resolve_iip(iip, op)
    if !fast_path
        op = build_operating_point(sys, op)
        fast_path = true
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
            fast_path
        )
        simplify_system = true
    else
        isys = generate_initializesystem(
            sys; op, initialization_eqs, check_units, time_dependent_init,
            guesses, algebraic_only, fast_path
        )
        simplify_system = true
    end

    # useful for `SteadyStateProblem` since `f` has to be autonomous and the
    # initialization should be too
    if !time_dependent_init
        idx = findfirst(isequal(get_iv(sys)), get_ps(isys))
        idx === nothing || deleteat!(get_ps(isys), idx)
    end

    if !is_split(sys)
        @set! isys.ps = mapreduce(collect, vcat, get_ps(isys))
    end
    if is_steadystateprob && time_dependent_init
        @set! isys.ps = filter(!isequal(get_iv(sys)::SymbolicT), get_ps(isys))
        binds = copy(parent(bindings(isys)))
        # Steady state problems can assume `t0 = 0`
        binds[get_iv(sys)::SymbolicT] = Symbolics.COMMON_ZERO
        @set! isys.bindings = ROSymmapT(binds)
    end
    pareqs, resteqs = find_all_parameter_equations(isys)
    @set! isys.eqs = resteqs
    if simplify_system
        isys = mtkcompile(isys; fully_determined, split = is_split(sys), initsys_mtkcompile_kwargs...)
    end
    for i in eachindex(pareqs)
        eq = pareqs[i]
        pareqs[i] = Symbolics.COMMON_ZERO ~ (eq.rhs - eq.lhs)
    end
    @set! isys.eqs = [equations(isys); pareqs]

    ts = get_tearing_state(isys)
    unassigned_vars = singular_check(ts)
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
    for (k, v) in bindings(sys)
        v === COMMON_MISSING || continue
        push!(uninit, k)
    end
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
        gs = ModelingToolkitBase.initial_conditions(isys)
        sys_gs = ModelingToolkitBase.guesses(sys)
        for k in uninit
            haskey(gs, k) && continue
            haskey(sys_gs, k) || continue
            gs[k] = sys_gs[k]
        end
        isys = complete(isys)
    end

    if t !== nothing
        op = copy(op)
        op[get_iv(sys)] = t
    end
    # Observed of `sys` aren't present in `isys` anymore, so this enables guess-propagation
    # to work properly.
    add_observed!(sys, ModelingToolkitBase.initial_conditions(isys))
    filter!(!Base.Fix2(===, COMMON_MISSING) ∘ last, op)
    TProb = get_initialization_problem_type(
        sys, isys; warn_initialize_determined,
        kwargs...
    )
    # Only forward `check_length` when the caller explicitly set it; otherwise let the
    # underlying problem type apply its own default (see the keyword's definition above).
    check_length_kw = check_length === nothing ? (;) : (; check_length)
    if TProb === SciMLBase.HomotopyProblem
        # The init system contains Modelica `homotopy` nodes: build a `HomotopyProblem`
        # so the initialization is solved by continuation from the `simplified` form.
        # Pass the resolved `_iip` explicitly (like the sibling `TProb{_iip}` branch): the
        # bare `HomotopyProblem(sys, op)` would re-derive in-place-ness from `op`, which is
        # a varmap and so always reads as in-place — wrong for an out-of-place system.
        SciMLBase.HomotopyProblem{_iip}(
            isys, op; kwargs..., check_length_kw...,
            build_initializeprob = false, is_initializeprob = true
        )
    else
        TProb{_iip}(
            isys, op; kwargs..., check_length_kw...,
            build_initializeprob = false, is_initializeprob = true
        )
    end
end

function overdetermined_initialization_message(neqs::Integer, nunknown::Integer, extra::AbstractString)
    return """
    Initialization system is overdetermined. $neqs equations for $nunknown unknowns. \
    Initialization will default to using least squares. $(extra)

    Call `analyze_initialization_jacobian(prob)` on the constructed problem to see which \
    equations are redundant (and which unknowns, if any, remain underdetermined).

    To suppress this warning, pass `warn_initialize_determined = false`. To turn this \
    warning into an error, pass `fully_determined = true`.
    """
end

function underdetermined_initialization_message(neqs::Integer, nunknown::Integer, extra::AbstractString)
    return """
    Initialization system is underdetermined. $neqs equations for $nunknown unknowns. \
    Initialization will default to using least squares. $(extra)

    Call `analyze_initialization_jacobian(prob)` on the constructed problem to see which \
    unknowns are underdetermined (and which equations, if any, are redundant).

    To suppress this warning, pass `warn_initialize_determined = false`. To turn this \
    warning into an error, pass `fully_determined = true`.
    """
end

"""
    $TYPEDSIGNATURES

Get the type of the initialization problem to use, given the system `sys`, initialization
system `isys` and arbitrary keyword arguments. Returns `LinearInitializationProblem` for an
affine init system, `NonlinearLeastSquaresProblem` for a non-square one, and — for a square
one — the type selected by `get_nonlinear_problem_type(isys)`: `SciMLBase.HomotopyProblem`
when `isys` contains Modelica `homotopy(actual, simplified)` nodes, otherwise
`NonlinearProblem`.
"""
function get_initialization_problem_type(
        sys::AbstractSystem, isys::AbstractSystem;
        warn_initialize_determined = true, kwargs...
    )
    neqs = length(equations(isys))
    nunknown = length(unknowns(isys))

    if warn_initialize_determined && neqs > nunknown
        @warn overdetermined_initialization_message(neqs, nunknown, "")
    end
    if warn_initialize_determined && neqs < nunknown
        @warn underdetermined_initialization_message(neqs, nunknown, "")
    end

    # Avoid using this for underdetermined systems
    return if isys isa System && nunknown > 0 && nunknown <= neqs && calculate_A_b(isys; throw = false) !== nothing
        LinearInitializationProblem
    elseif neqs == nunknown
        # Square nonlinear init system: select the concrete problem type with the same
        # dispatch as `AbstractNonlinearProblem(isys, op)`. When the init system carries
        # Modelica `homotopy(actual, simplified)` nodes that depend on the unknowns this
        # returns `SciMLBase.HomotopyProblem` (solved by continuation from the `simplified`
        # form), otherwise `NonlinearProblem`. (A `homotopy` node whose arguments are
        # constant in the unknowns collapses to an affine term and is correctly caught by
        # the `LinearInitializationProblem` branch above instead — there is nothing to
        # continue. `get_nonlinear_problem_type` is only defined for `System`; other
        # `AbstractSystem`s keep the previous `NonlinearProblem` behavior.)
        isys isa System ? get_nonlinear_problem_type(isys) : NonlinearProblem
    else
        NonlinearLeastSquaresProblem
    end
end

"""
    $TYPEDSIGNATURES

Return a list of possibly singular variables, given `get_tearing_state(sys)`.
"""
singular_check(::Nothing) = SymbolicT[]

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
    return println(io, underscore_to_D(collect(e.uninit), e.sys))
end

struct LinearInitializationProblem{iip} end

function LinearInitializationProblem{iip}(
        sys::AbstractSystem, op; u0_constructor = identity, kwargs...
    ) where {iip}
    # check_length = false allows using this for non-square systems
    linprob = LinearProblem{iip}(sys, op; u0_constructor, check_length = false, kwargs...)
    # Required for filling missing parameter values when this is an initialization
    # problem
    if state_values(linprob) === nothing
        linprob = remake(
            linprob;
            u0 = u0_constructor(ones(eltype(linprob.A), size(linprob.A, 2)))
        )
    end
    return SCCNonlinearProblem((linprob,), (Returns(nothing),), parameter_values(linprob), true; sys)
end
