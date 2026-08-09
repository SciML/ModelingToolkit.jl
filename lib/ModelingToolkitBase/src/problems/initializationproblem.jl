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

"""
    $(TYPEDSIGNATURES)

Equivalent to the keyword-argument-based `InitializationProblem`, given a pre-assembled
[`SciMLProblemOptions`](@ref). Public entry point for callers (namely
`maybe_build_initialization_problem`) that already hold an options struct, mirroring the
`(sys, opts::SciMLFunctionOptions)` methods on the `*Function` constructors. This is the
primary implementation; the keyword-argument-based method builds `opts` and delegates here.

`opts.check_length` is intentionally *not* used here: unlike the rest of `opts`, this
function's own `check_length` keyword is a nullable sentinel ("no opinion — let the
underlying problem constructor pick its default"), which has no equivalent in
`SciMLProblemOptions` (there, `check_length` is a concrete, always-set `Bool` describing a
different check on the outer system). Callers that want to override it must still pass it
as an explicit keyword.
"""
function InitializationProblem{iip}(
        sys::AbstractSystem, t, op, opts::SciMLProblemOptions; kwargs...
    ) where {iip}
    return InitializationProblem{iip, SciMLBase.AutoSpecialize}(sys, t, op, opts; kwargs...)
end

function InitializationProblem{iip, specialize}(
        sys::AbstractSystem, t, op, opts::SciMLProblemOptions;
        fast_path = false, guesses = [], check_length = nothing, kwargs...
    ) where {iip, specialize}
    (;
        warn_initialize_determined, initialization_eqs, fully_determined,
        check_initialization_units, allow_incomplete, algebraic_only,
        time_dependent_init, initsys_mtkcompile_kwargs, is_steadystateprob,
        u0_constructor, p_constructor, use_scc, missing_guess_value,
        warn_cyclic_dependency, circular_dependency_max_cycle_length,
        circular_dependency_max_cycles, init_compiler_options,
    ) = opts
    (; eval_expression, eval_module) = opts.fn_opts.codegen
    check_units = check_initialization_units

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
        isys = mtkcompile_initialization_system(
            isys, sys; fully_determined, initsys_mtkcompile_kwargs...
        )
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
    TProb = get_initialization_problem_type(sys, isys; warn_initialize_determined, use_scc)
    # Only forward `check_length` when the caller explicitly set it; otherwise let the
    # underlying problem type apply its own default (see the keyword's definition above).
    check_length_kw = check_length === nothing ? (;) : (; check_length)
    return TProb{_iip}(
        isys, op; kwargs..., check_length_kw...,
        u0_constructor, p_constructor, missing_guess_value,
        eval_expression, eval_module, warn_cyclic_dependency,
        circular_dependency_max_cycle_length, circular_dependency_max_cycles,
        compiler_options = init_compiler_options,
        build_initializeprob = false, is_initializeprob = true
    )
end

"""
    $TYPEDSIGNATURES

`mtkcompile` the initialization system `isys` of `sys`, rewriting structural errors so
they say that it is the initialization system which is unbalanced.

Structural analysis has no way of knowing that the system handed to it is an
initialization system, so on its own it reports the imbalance in terms which read as if
the model itself were at fault. Since this is almost always hit through a problem
constructor rather than a direct `mtkcompile` call, the message is rewritten to name the
initialization system, count the missing equations and list the ways of supplying them.
"""
function mtkcompile_initialization_system(
        isys::AbstractSystem, sys::AbstractSystem; fully_determined, kwargs...
    )
    try
        return mtkcompile(isys; fully_determined, split = is_split(sys), kwargs...)
    catch err
        newerr = with_initialization_context(err, isys, sys; kwargs...)
        newerr === err && rethrow()
        rethrow(newerr)
    end
end

"""
Structural analysis exceptions which report an unbalanced or singular system. These are
defined both here and in StateSelection.jl, which ModelingToolkitBase does not depend on,
so they are identified by name to cover whichever of the two was thrown.
"""
const UNBALANCED_SYSTEM_EXCEPTIONS = (
    :ExtraVariablesSystemException, :ExtraEquationsSystemException, :InvalidSystemException,
)

"""
    $TYPEDSIGNATURES

Given an error `err` thrown while compiling the initialization system `isys` of `sys`,
return an equivalent error whose message identifies it as an initialization failure.
Returns `err` unchanged if it is not an error about the balance of a system.
"""
function with_initialization_context(
        err, isys::AbstractSystem, sys::AbstractSystem; kwargs...
    )
    T = typeof(err)
    kind = nameof(T)
    kind in UNBALANCED_SYSTEM_EXCEPTIONS || return err
    fieldcount(T) == 1 && fieldtype(T, 1) <: AbstractString || return err
    msg = getfield(err, 1)
    # `InvalidSystemException` also reports problems which have nothing to do with the
    # balance of the system, and for which the guidance below would be wrong.
    kind === :InvalidSystemException && !occursin("singular", msg) && return err
    return T(
        initialization_structure_error_message(
            kind, msg, initialization_system_size(isys, sys; kwargs...)
        )
    )
end

"""
    $TYPEDSIGNATURES

Return `(neqs, nunknowns)` of the compiled initialization system `isys` of `sys`, or
`nothing` if it cannot be compiled without the balance check. Used to report how many
equations an unbalanced initialization system is missing, in the same terms as
`underdetermined_initialization_message`.
"""
function initialization_system_size(isys::AbstractSystem, sys::AbstractSystem; kwargs...)
    return try
        compiled = mtkcompile(
            isys; fully_determined = false, split = is_split(sys), kwargs...
        )
        (length(equations(compiled)), length(unknowns(compiled)))
    catch
        nothing
    end
end

const INITIALIZATION_DOCS_URL = "https://docs.sciml.ai/ModelingToolkit/stable/tutorials/initialization/"

const INITIALIZATION_CONTEXT_MESSAGE = """
This is an error in the initialization system, not in the equations being integrated. \
The initialization system solves for the value of every unknown of the model, and of \
their derivatives, at the initial time. It needs as many equations as there are such \
values to solve for. Initial values given for the model become equations of this system.
"""

const INITIALIZATION_UNDERDETERMINED_MESSAGE = """
Each missing equation is supplied by giving one more initial value, or one more equation \
relating initial values. Guesses are starting values for the initialization solve rather \
than constraints on it, so adding a guess does not supply one.
"""

const INITIALIZATION_UNDERDETERMINED_API_MESSAGE = """
In ModelingToolkit, any of the following supplies one missing equation:

  * Give an initial value in the problem constructor, as in
    `ODEProblem(sys, [x => 1.0], tspan)`.
  * Give one in the model, via `initial_conditions` or `bindings`.
  * Add an equation relating initial values, via the `initialization_eqs` keyword argument \
of the system or of the problem constructor.

To solve an underdetermined initialization in a least squares sense instead, pass \
`fully_determined = false` to the problem constructor; the guesses then decide which of \
the possible initial states is found.
"""

const INITIALIZATION_OVERDETERMINED_MESSAGE = """
Remove initial values or equations relating them until as many remain as there are \
unknowns to solve for.
"""

const INITIALIZATION_OVERDETERMINED_API_MESSAGE = """
In ModelingToolkit, these come from `initial_conditions`, `bindings` and \
`initialization_eqs`, and from the operating point passed to the problem constructor. To \
solve an overdetermined initialization in a least squares sense instead, pass \
`fully_determined = false` to the problem constructor. That only finds an initial state \
if the extra equations are consistent with the rest.
"""

const INITIALIZATION_SINGULAR_MESSAGE = """
There are as many equations as unknowns, but they do not determine a unique initial state: \
at least one equation is redundant given the others, leaving some unknown undetermined. \
This usually means two of the given initial conditions are related by an equation of the \
model, so that one of them carries no new information. Giving an initial condition for a \
different variable, or for a derivative, in place of the redundant one resolves it.
"""

const INITIALIZATION_NOTATION_MESSAGE = """
Note that variables named with a `ˍt` suffix are derivatives at the initial time: `xˍt` \
is `D(x)`.
"""

"""
    $TYPEDSIGNATURES

Build the message of a structural error from the initialization system. `kind` is the name
of the exception type thrown by structural analysis, `msg` its original message, and `size`
the `(neqs, nunknowns)` of the initialization system, if known.
"""
function initialization_structure_error_message(
        kind::Symbol, msg::AbstractString, size::Union{Nothing, Tuple{Int, Int}}
    )
    io = IOBuffer()
    if kind === :ExtraVariablesSystemException
        headline = "Initialization system is underdetermined."
        remedy = INITIALIZATION_UNDERDETERMINED_MESSAGE
        api_remedy = INITIALIZATION_UNDERDETERMINED_API_MESSAGE
    elseif kind === :ExtraEquationsSystemException
        headline = "Initialization system is overdetermined."
        remedy = INITIALIZATION_OVERDETERMINED_MESSAGE
        api_remedy = INITIALIZATION_OVERDETERMINED_API_MESSAGE
    else
        headline = "Initialization system is structurally singular."
        remedy = INITIALIZATION_SINGULAR_MESSAGE
        api_remedy = nothing
    end
    println(io, headline, '\n')
    println(io, INITIALIZATION_CONTEXT_MESSAGE)
    println(io, msg, '\n')
    deficit = initialization_deficit_message(kind, size)
    deficit === nothing || println(io, deficit)
    println(io, remedy)
    if show_api_guidance()
        api_remedy === nothing || println(io, api_remedy)
        println(io, "See $INITIALIZATION_DOCS_URL for more information.")
    end
    print(io, INITIALIZATION_NOTATION_MESSAGE)
    return String(take!(io))
end

"""
    $TYPEDSIGNATURES

Return a sentence stating by how many equations an initialization system of size
`size = (neqs, nunknowns)` is out of balance, or `nothing` if `size` is unknown or does
not corroborate the imbalance that structural analysis reported.

Only the difference is reported, not the counts themselves: `size` is measured on the
system compiled without the balance check, which simplifies further than the compilation
that failed, so its counts do not line up with the ones structural analysis reports.
"""
function initialization_deficit_message(kind::Symbol, size::Union{Nothing, Tuple{Int, Int}})
    size === nothing && return nothing
    neqs, nunknown = size
    if kind === :ExtraVariablesSystemException && nunknown > neqs
        n = nunknown - neqs
        verb = n == 1 ? "is" : "are"
        return "$n more $(pluralize(n, "equation")) $verb needed to determine the initial state.\n"
    elseif kind === :ExtraEquationsSystemException && neqs > nunknown
        n = neqs - nunknown
        return "There $(n == 1 ? "is" : "are") $n $(pluralize(n, "equation")) too many.\n"
    end
    return nothing
end

pluralize(n::Integer, word::AbstractString) = n == 1 ? word : word * "s"

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
        u0_constructor = identity, p_constructor = identity, use_scc = true,
        missing_guess_value = default_missing_guess_value(),
        eval_expression = false, eval_module = @__MODULE__,
        warn_cyclic_dependency = false,
        circular_dependency_max_cycle_length = length(all_symbols(sys)),
        circular_dependency_max_cycles = 10,
        compiler_options::CompilerOptions = CompilerOptions(),
        kwargs...
    ) where {iip, specialize}
    fn_opts = SciMLFunctionOptions(; eval_expression, eval_module)
    opts = SciMLProblemOptions(
        sys;
        fn_opts, warn_initialize_determined, initialization_eqs, fully_determined,
        check_initialization_units = check_units, allow_incomplete, algebraic_only,
        time_dependent_init, initsys_mtkcompile_kwargs, is_steadystateprob,
        u0_constructor, p_constructor, use_scc, missing_guess_value,
        warn_cyclic_dependency, circular_dependency_max_cycle_length,
        circular_dependency_max_cycles, init_compiler_options = compiler_options,
        build_initializeprob = false,
    )
    return InitializationProblem{iip, specialize}(
        sys, t, op, opts; fast_path, guesses, check_length, kwargs...
    )
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
        # continue.)
        get_nonlinear_problem_type(isys)
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
