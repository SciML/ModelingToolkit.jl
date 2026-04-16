abstract type AbstractCallback end

function has_functional_affect(cb)
    return affects(cb) isa ImperativeAffect
end

struct SymbolicAffect
    affect::Vector{Equation}
    discrete_parameters::Vector{SymbolicT}
end

@noinline function depwarn_alg_eqs()
    return Base.depwarn(
        """
        The `alg_eqs` keyword for callback affects is deprecated. The equations can \
        simply be passed as equations of the affect.
        """,
        :callback_alg_eqs
    )
end

function SymbolicAffect(
        affect::Vector{Equation}; alg_eqs = Equation[],
        discrete_parameters = SymbolicT[], kwargs...
    )
    if !isempty(alg_eqs)
        depwarn_alg_eqs()
        append!(affect, alg_eqs)
    end
    if symbolic_type(discrete_parameters) !== NotSymbolic()
        discrete_parameters = SymbolicT[unwrap(discrete_parameters)]
    elseif !(discrete_parameters isa Vector{SymbolicT})
        _discs = SymbolicT[]
        for p in discrete_parameters
            push!(_discs, unwrap(p))
        end
        discrete_parameters = _discs
    end
    return SymbolicAffect(affect, discrete_parameters)
end
function SymbolicAffect(affect::SymbolicAffect; kwargs...)
    return SymbolicAffect(
        affect.affect; discrete_parameters = affect.discrete_parameters, kwargs...
    )
end
SymbolicAffect(affect; kwargs...) = make_affect(affect; kwargs...)

"""
    AssignmentAffect(assignments)

A special shorthand affect usable in symbolic callbacks. `assignments` is an array of
`Pair`s, where the LHS of each pair is the variable/discrete to assign and the RHS is
the value to assign it. The RHS is evaluated using the values of variables and parameters
before the callback. In other words, it is evaluated as if the entire expression is enclosed
in [`Pre`](@ref). In addition, if the LHS is created using `@discretes` it will automatically
be added to the list of `discrete_parameters` typically provided when manually creating
events.
"""
function AssignmentAffect(affect::AbstractArray{T}) where {T <: Pair}
    discrete_parameters = SymbolicT[]
    new_affect = Equation[]
    for (lhs, rhs) in affect
        lhs = unwrap(lhs)::SymbolicT
        rhs = unwrap(rhs)
        push!(new_affect, lhs ~ Pre(rhs))
        if first(getmetadata(lhs, Symbolics.VariableSource)::NTuple{2, Symbol}) == :discretes
            push!(discrete_parameters, lhs)
        end
    end
    return SymbolicAffect(new_affect; discrete_parameters)
end

function SymbolicAffect(affect::Vector; discrete_parameters = SymbolicT[], kwargs...)
    assigns = Pair{SymbolicT, SymbolicT}[]
    eqs = Equation[]
    for aff in affect
        if aff isa Equation
            push!(eqs, aff)
        elseif aff isa Pair
            push!(assigns, aff)
        else
            throw(
                ArgumentError(
                    lazy"""
                    Unrecognized affect format $aff. Affects must be equations or `Pair`s.
                    """
                )
            )
        end
    end
    aff1 = SymbolicAffect(eqs; discrete_parameters, kwargs...)
    aff2 = AssignmentAffect(assigns; kwargs...)
    return SymbolicAffect([aff1.affect; aff2.affect], [aff1.discrete_parameters; aff2.discrete_parameters])
end

function (s::SymbolicUtils.Substituter)(aff::SymbolicAffect)
    return SymbolicAffect(s(aff.affect), s(aff.discrete_parameters))
end

discretes(affect::SymbolicAffect) = affect.discrete_parameters

struct AffectSystem
    """The internal implicit discrete system whose equations are solved to obtain values after the affect."""
    system::AbstractSystem
    """Unknowns of the parent ODESystem whose values are modified or accessed by the affect."""
    unknowns::Vector{SymbolicT}
    """Parameters of the parent ODESystem whose values are accessed by the affect."""
    parameters::Vector{SymbolicT}
    """Parameters of the parent ODESystem whose values are modified by the affect."""
    discretes::Vector{SymbolicT}
end

function (s::SymbolicUtils.Substituter)(aff::AffectSystem)
    sys = aff.system
    @set! sys.eqs = s(get_eqs(sys))
    @set! sys.parameter_dependencies = (get_parameter_dependencies(sys))
    @set! sys.defaults = Dict([k => s(v) for (k, v) in defaults(sys)])
    @set! sys.guesses = Dict([k => s(v) for (k, v) in guesses(sys)])
    @set! sys.unknowns = s(get_unknowns(sys))
    @set! sys.ps = s(get_ps(sys))
    return AffectSystem(sys, s(aff.unknowns), s(aff.parameters), s(aff.discretes))

end

function AffectSystem(spec::SymbolicAffect; iv = nothing, alg_eqs = Equation[], kwargs...)
    affect = spec.affect
    if !isempty(alg_eqs)
        depwarn_alg_eqs()
        affect = [affect; alg_eqs]
    end
    return AffectSystem(affect; iv, discrete_parameters = spec.discrete_parameters, kwargs...)
end

function AffectSystem(
        affect::Vector{Equation}; discrete_parameters = SymbolicT[],
        iv = nothing, extra_eqs = Equation[], kwargs...
    )
    isempty(affect) && return nothing
    if isnothing(iv)
        iv = t_nounits
        @warn "No independent variable specified. Defaulting to t_nounits."
    end
    affect = [affect; extra_eqs]

    discrete_parameters = SymbolicAffect(affect; discrete_parameters).discrete_parameters

    for p in discrete_parameters
        SU.query(isequal(unwrap(iv)), unwrap(p)) ||
            error("Non-time dependent parameter $p passed in as a discrete. Must be declared as @parameters $p(t).")
    end

    dvs = OrderedSet{SymbolicT}()
    params = OrderedSet{SymbolicT}()
    _varsbuf = Set{SymbolicT}()
    for eq in affect
        collect_vars!(dvs, params, eq, iv, Pre)
        empty!(_varsbuf)
        SU.search_variables!(_varsbuf, eq; is_atomic = OperatorIsAtomic{Pre}())
        filter!(x -> iscall(x) && operation(x) === Pre(), _varsbuf)
        union!(params, _varsbuf)
        diffvs = collect_applied_operators(eq, Differential)
        union!(dvs, diffvs)
    end
    pre_params = filter(haspre, params)
    sys_params = SymbolicT[]
    disc_ps_set = Set{SymbolicT}(discrete_parameters)
    disc_ps_set = gather_array_params(disc_ps_set)
    discrete_parameters = collect(disc_ps_set)
    for p in params
        p in disc_ps_set && continue
        p in pre_params && continue
        push!(sys_params, p)
    end
    discretes = discrete_parameters
    dvs = collect(dvs)
    _dvs = map(default_toterm, dvs)

    subs = Dict{SymbolicT, SymbolicT}(zip(dvs, _dvs))
    affect = substitute(affect, subs)

    ps = collect(gather_array_params(union(pre_params, sys_params)))

    @named affectsys = System(
        affect, iv, collect(union(_dvs, discretes)),
        ps; is_discrete = true
    )
    # This `@invokelatest` should not be necessary, but it works around the inference bug
    # in https://github.com/JuliaLang/julia/issues/59943. Remove it at your own risk, the
    # bug took weeks to reduce to an MWE.
    affectsys = (@invokelatest mtkcompile(affectsys; fully_determined = nothing))::System
    # get accessed parameters p from Pre(p) in the callback parameters
    accessed_params = Vector{SymbolicT}(filter(isparameter, map(unPre, collect(pre_params))))
    union!(accessed_params, sys_params)

    # add scalarized unknowns to the map.
    _obs = observed(unhack_system(affectsys))
    _dvs = vcat(unknowns(affectsys), map(eq -> eq.lhs, _obs))
    _dvs = __safe_scalarize_vars(_dvs)
    _discs = __safe_scalarize_vars(discretes)
    setdiff!(_dvs, _discs)
    return AffectSystem(affectsys, _dvs, accessed_params, discrete_parameters)
end

function __safe_scalarize_vars(vars::Vector{SymbolicT})
    _vars = SymbolicT[]
    for v in vars
        sh = SU.shape(v)::SU.ShapeVecT
        if isempty(sh)
            push!(_vars, v)
            continue
        end
        for i in SU.stable_eachindex(v)
            push!(_vars, v[i])
        end
    end
    return _vars
end

safe_vec(@nospecialize(x)) = x isa SymbolicT ? [x] : vec(x::Array{SymbolicT})

system(a::AffectSystem) = a.system
discretes(a::AffectSystem) = a.discretes
unknowns(a::AffectSystem) = a.unknowns
parameters(a::AffectSystem) = a.parameters
all_equations(a::AffectSystem) = vcat(equations(system(a)), observed(system(a)))

function Base.show(iio::IO, aff::AffectSystem)
    println(iio, "Affect system defined by equations:")
    eqs = all_equations(aff)
    return show(iio, eqs)
end

function Base.:(==)(a1::AffectSystem, a2::AffectSystem)
    return isequal(system(a1), system(a2)) &&
        isequal(discretes(a1), discretes(a2)) &&
        isequal(unknowns(a1), unknowns(a2)) &&
        isequal(parameters(a1), parameters(a2))
end

function Base.hash(a::AffectSystem, s::UInt)
    s = hash(system(a), s)
    s = hash(unknowns(a), s)
    s = hash(parameters(a), s)
    return hash(discretes(a), s)
end

function SU.search_variables!(vars, aff::AffectSystem; kwargs...)
    SU.search_variables!(vars, unknowns(aff); kwargs...)
    SU.search_variables!(vars, parameters(aff); kwargs...)
    return SU.search_variables!(vars, discretes(aff); kwargs...)
end

"""
    Pre(x)

The `Pre` operator. Used by the callback system to indicate the value of a parameter or variable
before the callback is triggered.
"""
struct Pre <: Symbolics.Operator end
Pre(x) = Pre()(x)
SymbolicUtils.promote_symtype(::Type{Pre}, T) = T
SymbolicUtils.isbinop(::Pre) = false
Base.nameof(::Pre) = :Pre
Base.show(io::IO, x::Pre) = print(io, "Pre")
unPre(x::Num) = unPre(unwrap(x))
unPre(x::Symbolics.Arr) = unPre(unwrap(x))
unPre(x::SymbolicT) = (iscall(x) && operation(x) isa Pre) ? only(arguments(x)) : x
distribute_shift_into_operator(::Pre) = false

(p::Pre)(x::Num) = Num(p(unwrap(x)))
(p::Pre)(x::Symbolics.Arr{T, N}) where {T, N} = Symbolics.Arr{T, N}(p(unwrap(x)))
(p::Pre)(x::Symbolics.SymStruct{T}) where {T} = Symbolics.SymStruct{T}(p(unwrap(x)))
(p::Pre)(x::Symbolics.CallAndWrap{T}) where {T} = Symbolics.CallAndWrap{T}(p(unwrap(x)))
function (p::Pre)(x::SymbolicT)
    iscall(x) || return x
    return Moshi.Match.@match x begin
        BSImpl.Term(; f) && if f isa Pre end => return x
        BSImpl.Term(; f) && if f isa Differential end => begin
            return p(default_toterm(x))
        end
        BSImpl.Term(; f, args, type, shape) && if f === getindex end => begin
            arrpre = p(args[1])
            Moshi.Match.@match arrpre begin
                BSImpl.Term(; f = f2) && if f2 isa Pre end => begin
                    newargs = SArgsT((x,))
                    return toparam(BSImpl.Term{VartypeT}(p, newargs; type, shape))
                end
                _ => begin
                    newargs = copy(parent(args))
                    newargs[1] = arrpre
                    return toparam(BSImpl.Term{VartypeT}(f, newargs; type, shape))
                end
            end
        end
        BSImpl.Term(; f, type, shape) && if f isa SymbolicT && !SU.is_function_symbolic(f) end => begin
            return toparam(BSImpl.Term{VartypeT}(p, SArgsT((x,)); type, shape))
        end
        _ => begin
            op = operation(x)
            args = map(p, arguments(x))
            return toparam(maketerm(SymbolicT, op, args, nothing; type = symtype(x)))
        end
    end
end
(::Pre)(x) = x
haspre(eq::Equation) = haspre(eq.lhs) || haspre(eq.rhs)
haspre(O) = recursive_hasoperator(Pre, O)

function validate_operator(op::Pre, args, iv; context = nothing)
end

###############################
###### Continuous events ######
###############################
const Affect = Union{AffectSystem, ImperativeAffect}

"""
    SymbolicContinuousCallback(eqs::Vector{Equation}, affect = nothing, iv = nothing; 
                               affect_neg = affect, initialize = nothing, finalize = nothing,
                               rootfind = SciMLBase.LeftRootFind, initialize_save_discretes = true)

A [`ContinuousCallback`](@ref SciMLBase.ContinuousCallback) specified symbolically. Takes a vector of equations `eq`
as well as the positive-edge `affect` and negative-edge `affect_neg` that apply when *any* of `eq` are satisfied.
By default `affect_neg = affect`; to only get rising edges specify `affect_neg = nothing`.

Assume without loss of generality that the equation is of the form `c(u,p,t) ~ 0`; we denote the integrator state as `i.u`.
For compactness, we define `prev_sign = sign(c(u[t-1], p[t-1], t-1))` and `cur_sign = sign(c(u[t], p[t], t))`.
A condition edge will be detected and the callback will be invoked iff `prev_sign * cur_sign <= 0`.
The positive edge `affect` will be triggered iff an edge is detected and if `prev_sign < 0`; similarly, `affect_neg` will be
triggered iff an edge is detected and `prev_sign > 0`.

Inter-sample condition activation is not guaranteed; for example if we use the dirac delta function as `c` to insert a
sharp discontinuity between integrator steps (which in this example would not normally be identified by adaptivity) then the condition is not
guaranteed to be triggered.

Once detected the integrator will "wind back" through a root-finding process to identify the point when the condition became active; the method used
is specified by `rootfind` from [`SciMLBase.RootfindOpt`](@ref). If we denote the time when the condition becomes active as `tc`,
the value in the integrator after windback will be:
* `u[tc-epsilon], p[tc-epsilon], tc` if `LeftRootFind` is used,
* `u[tc+epsilon], p[tc+epsilon], tc` if `RightRootFind` is used,
* or `u[t], p[t], t` if `NoRootFind` is used.
For example, if we want to detect when an unknown variable `x` satisfies `x > 0` using the condition `x ~ 0` on a positive edge (that is, `D(x) > 0`),
then left root finding will get us `x=-epsilon`, right root finding `x=epsilon` and no root finding will produce whatever the next step of the integrator was after
it passed through 0.

Multiple callbacks in the same system with different `rootfind` operations will be grouped
by their `rootfind` value into separate VectorContinuousCallbacks in the enumeration order of `SciMLBase.RootfindOpt`. This may cause some callbacks to not fire if several become
active at the same instant. See the `SciMLBase` documentation for more information on the semantic rules.

Affects (i.e. `affect` and `affect_neg`) can be specified as either:
* A list of equations that should be applied when the callback is triggered (e.g. `x ~ 3, y ~ 7`) which must be of the form `unknown ~ observed value` where each `unknown` appears only once. Equations will be applied in the order that they appear in the vector; parameters and state updates will become immediately visible to following equations. Instead of
equations, elements can also be `Pair`s. These are parsed according to the [`AssignmentAffect`](@ref) semantics
and combined with the remaining `Equation`s.
* An [`AssignmentAffect`](@ref).
* A tuple `(f!, unknowns, read_parameters, modified_parameters, ctx)`, where:
    + `f!` is a function with signature `(integ, u, p, ctx)` that is called with the integrator, a state *index* vector `u` derived from `unknowns`, a parameter *index* vector `p` derived from `read_parameters`, and the `ctx` that was given at construction time. Note that `ctx` is aliased between instances.
    + `unknowns` is a vector of symbolic unknown variables and optionally their aliases (e.g. if the model was defined with `@variables x(t)` then a valid value for `unknowns` would be `[x]`). A variable can be aliased with a pair `x => :y`. The indices of these `unknowns` will be passed to `f!` in `u` in a named tuple; in the earlier example, if we pass `[x]` as `unknowns` then `f!` can access `x` as `integ.u[u.x]`. If no alias is specified the name of the index will be the symbol version of the variable name.
    + `read_parameters` is a vector of the parameters that are *used* by `f!`. Their indices are passed to `f` in `p` similarly to the indices of `unknowns` passed in `u`.
    + `modified_parameters` is a vector of the parameters that are *modified* by `f!`. Note that a parameter will not appear in `p` if it only appears in `modified_parameters`; it must appear in both `parameters` and `modified_parameters` if it is used in the affect definition.
    + `ctx` is a user-defined context object passed to `f!` when invoked. This value is aliased for each problem.
* A [`ImperativeAffect`](@ref); refer to its documentation for details.

`reinitializealg` is used to set how the system will be reinitialized after the callback. 
- Symbolic affects have reinitialization built in. In this case the algorithm will default to SciMLBase.NoInit(), and should **not** be provided.
- Functional and imperative affects will default to SciMLBase.CheckInit(), which will error if the system is not properly reinitialized after the callback. If your system is a DAE, pass in an algorithm like SciMLBase.BrownBasicFullInit() to properly re-initialize.

`initialize_save_discretes` is a flag indicating whether the discrete variables modified by this
callback should be saved at the start of the integration (when the `initialize` runs).

Initial and final affects can also be specified identically to positive and negative edge affects. Initialization affects
will run as soon as the solver starts, while finalization affects will be executed after termination.
"""
struct SymbolicContinuousCallback <: AbstractCallback
    conditions::Vector{Equation}
    affect::Union{Affect, SymbolicAffect, Nothing}
    affect_neg::Union{Affect, SymbolicAffect, Nothing}
    initialize::Union{Affect, SymbolicAffect, Nothing}
    finalize::Union{Affect, SymbolicAffect, Nothing}
    rootfind::Union{Nothing, SciMLBase.RootfindOpt}
    reinitializealg::SciMLBase.DAEInitializationAlgorithm
    zero_crossing_id::Symbol
    initialize_save_discretes::Bool
end

function SymbolicContinuousCallback(
        conditions::Union{Equation, Vector{Equation}},
        affect = nothing;
        affect_neg = affect,
        initialize = nothing,
        finalize = nothing,
        rootfind = SciMLBase.LeftRootFind,
        reinitializealg = nothing,
        zero_crossing_id = gensym(),
        initialize_save_discretes = true,
        kwargs...
    )
    conditions = (conditions isa AbstractVector) ? conditions : [conditions]

    if isnothing(reinitializealg)
        if any(
                a -> a isa ImperativeAffect,
                [affect, affect_neg, initialize, finalize]
            )
            reinitializealg = SciMLBase.CheckInit()
        else
            reinitializealg = SciMLBase.NoInit()
        end
    end

    return SymbolicContinuousCallback(
        conditions, SymbolicAffect(affect; kwargs...),
        SymbolicAffect(affect_neg; kwargs...),
        SymbolicAffect(initialize; kwargs...), SymbolicAffect(
            finalize; kwargs...
        ),
        rootfind, reinitializealg, zero_crossing_id, initialize_save_discretes
    )
end # Default affect to nothing

function SymbolicContinuousCallback(p::Pair, args...; kwargs...)
    return SymbolicContinuousCallback(p[1], p[2], args...; kwargs...)
end
SymbolicContinuousCallback(cb::SymbolicContinuousCallback, args...; kwargs...) = cb
SymbolicContinuousCallback(cb::Nothing; kwargs...) = nothing
SymbolicContinuousCallback(cb::Nothing, args...; kwargs...) = nothing
function SymbolicContinuousCallback(cb::Tuple, args...; kwargs...)
    return if length(cb) == 2
        SymbolicContinuousCallback(cb[1]; kwargs..., cb[2]...)
    else
        error("Malformed tuple specifying callback. Should be a condition => affect pair, followed by a vector of kwargs.")
    end
end

function complete(cb::SymbolicContinuousCallback; kwargs...)
    return SymbolicContinuousCallback(
        cb.conditions, make_affect(cb.affect; kwargs...),
        make_affect(cb.affect_neg; kwargs...), make_affect(cb.initialize; kwargs...),
        make_affect(cb.finalize; kwargs...), cb.rootfind, cb.reinitializealg,
        cb.zero_crossing_id, cb.initialize_save_discretes
    )
end

make_affect(affect::SymbolicAffect; kwargs...) = AffectSystem(affect; kwargs...)
make_affect(affect::Nothing; kwargs...) = nothing
make_affect(affect::Tuple; kwargs...) = ImperativeAffect(affect...)
make_affect(affect::NamedTuple; kwargs...) = ImperativeAffect(; affect...)
make_affect(affect::Affect; kwargs...) = affect
make_affect(affect::Vector{Equation}; kwargs...) = AffectSystem(affect; kwargs...)

function make_affect(affect; kwargs...)
    error("Malformed affect $(affect). This should be a vector of equations or a tuple specifying a functional affect.")
end

function Base.show(io::IO, cb::AbstractCallback)
    indent = get(io, :indent, 0)
    iio = IOContext(io, :indent => indent + 1)
    is_discrete(cb) ? print(io, "SymbolicDiscreteCallback(") :
        print(io, "SymbolicContinuousCallback(")
    print(iio, "Conditions:")
    show(iio, equations(cb))
    print(iio, "; ")
    if affects(cb) != nothing
        print(iio, "Affect:")
        show(iio, affects(cb))
        print(iio, ", ")
    end
    if !is_discrete(cb) && affect_negs(cb) != nothing
        print(iio, "Negative-edge affect:")
        show(iio, affect_negs(cb))
        print(iio, ", ")
    end
    if initialize_affects(cb) != nothing
        print(iio, "Initialization affect:")
        show(iio, initialize_affects(cb))
        print(iio, ", ")
    end
    if finalize_affects(cb) != nothing
        print(iio, "Finalization affect:")
        show(iio, finalize_affects(cb))
    end
    return print(iio, ")")
end

function Base.show(io::IO, mime::MIME"text/plain", cb::AbstractCallback)
    indent = get(io, :indent, 0)
    iio = IOContext(io, :indent => indent + 1)
    is_discrete(cb) ? println(io, "SymbolicDiscreteCallback:") :
        println(io, "SymbolicContinuousCallback:")
    println(iio, "Conditions:")
    show(iio, mime, equations(cb))
    print(iio, "\n")
    if affects(cb) != nothing
        println(iio, "Affect:")
        show(iio, mime, affects(cb))
        print(iio, "\n")
    end
    if !is_discrete(cb) && affect_negs(cb) != nothing
        print(iio, "Negative-edge affect:\n")
        show(iio, mime, affect_negs(cb))
        print(iio, "\n")
    end
    if initialize_affects(cb) != nothing
        println(iio, "Initialization affect:")
        show(iio, mime, initialize_affects(cb))
        print(iio, "\n")
    end
    return if finalize_affects(cb) != nothing
        println(iio, "Finalization affect:")
        show(iio, mime, finalize_affects(cb))
        print(iio, "\n")
    end
end

function SU.search_variables!(vars, cb::AbstractCallback; kwargs...)
    if symbolic_type(conditions(cb)) isa NotSymbolic
        SU.search_variables!(vars, conditions(cb); kwargs...)
    end
    affs = affects(cb)
    affs === nothing || SU.search_variables!(vars, affs; kwargs...)
    affs = initialize_affects(cb)
    affs === nothing || SU.search_variables!(vars, affs; kwargs...)
    affs = finalize_affects(cb)
    affs === nothing || SU.search_variables!(vars, affs; kwargs...)
    return is_discrete(cb) || SU.search_variables!(vars, affect_negs(cb); kwargs...)
end

################################
######## Discrete events #######
################################
"""
    SymbolicDiscreteCallback(conditions::Vector{Equation}, affect = nothing, iv = nothing;
                             initialize = nothing, finalize = nothing)

A callback that triggers at the first timestep that the conditions are satisfied.

The condition can be one of: 
- Δt::Real              - periodic events with period Δt
- ts::Vector{Real}      - events trigger at these preset times given by `ts`
- eqs::Vector{SymbolicT} - events trigger when the condition evaluates to true
- A `SciMLBase.Clock(dt; phase)` - events trigger with period `dt` and phase `phase`.
  Note that this form will ignore the `initialize_save_discretes` keyword argument in the
  interest of correctness. The callback will trigger and save at `tspan[1]` if the clock
  would tick at `tspan[1]`.

Arguments: 
- iv: The independent variable of the system. This must be specified if the independent variable appears in one of the equations explicitly, as in x ~ t + 1.
"""
struct SymbolicDiscreteCallback <: AbstractCallback
    conditions::Union{Number, Vector{<:Number}, SymbolicT, SciMLBase.TimeDomain}
    affect::Union{Affect, SymbolicAffect, Nothing}
    initialize::Union{Affect, SymbolicAffect, Nothing}
    finalize::Union{Affect, SymbolicAffect, Nothing}
    reinitializealg::SciMLBase.DAEInitializationAlgorithm
    initialize_save_discretes::Bool
end

function SymbolicDiscreteCallback(
        condition::Union{SymbolicT, Number, Vector{<:Number}, SciMLBase.TimeDomain},
        affect = nothing; initialize = nothing, finalize = nothing,
        reinitializealg = nothing, initialize_save_discretes = true, kwargs...
    )
    # Manual error check (to prevent events like `[X < 5.0] => [X ~ Pre(X) + 10.0]` from being created).
    (condition isa Vector) && (eltype(condition) <: Num) &&
        error("Vectors of symbolic conditions are not allowed for `SymbolicDiscreteCallback`.")
    if condition isa SciMLBase.TimeDomain
        if !SciMLBase.isclock(condition)
            throw(
                ArgumentError("Clock given to `SymbolicDiscreteCallback` must be a `SciMLBase.Clock`.")
            )
        end
        c = condition
    else
        @assert !(condition isa SymbolicT && symtype(condition) != Bool)
        c = is_timed_condition(condition) ? condition : value(scalarize(condition))

        c = is_timed_condition(condition) ? condition : value(scalarize(condition))
    end
    if isnothing(reinitializealg)
        if any(
                a -> a isa ImperativeAffect,
                [affect, initialize, finalize]
            )
            reinitializealg = SciMLBase.CheckInit()
        else
            reinitializealg = SciMLBase.NoInit()
        end
    end
    return SymbolicDiscreteCallback(
        c, SymbolicAffect(affect; kwargs...),
        SymbolicAffect(initialize; kwargs...),
        SymbolicAffect(finalize; kwargs...), reinitializealg,
        initialize_save_discretes
    )
end # Default affect to nothing

function SymbolicDiscreteCallback(p::Pair, args...; kwargs...)
    return SymbolicDiscreteCallback(p[1], p[2], args...; kwargs...)
end
SymbolicDiscreteCallback(cb::SymbolicDiscreteCallback, args...; kwargs...) = cb
SymbolicDiscreteCallback(cb::Nothing, args...; kwargs...) = nothing
function SymbolicDiscreteCallback(cb::Tuple, args...; kwargs...)
    return if length(cb) == 2
        SymbolicDiscreteCallback(cb[1]; cb[2]...)
    else
        error("Malformed tuple specifying callback. Should be a condition => affect pair, followed by a vector of kwargs.")
    end
end

function complete(cb::SymbolicDiscreteCallback; kwargs...)
    return SymbolicDiscreteCallback(
        cb.conditions, make_affect(cb.affect; kwargs...),
        make_affect(cb.initialize; kwargs...),
        make_affect(cb.finalize; kwargs...), cb.reinitializealg,
        cb.initialize_save_discretes
    )
end

function is_timed_condition(condition::T) where {T}
    return if T === Num
        false
    elseif T <: Real
        true
    elseif T <: AbstractVector
        eltype(condition) <: Real
    elseif T <: SciMLBase.TimeDomain
        true
    else
        false
    end
end

to_cb_vector(cbs::Vector{<:AbstractCallback}; kwargs...) = cbs
to_cb_vector(cbs::Union{Nothing, Vector{Nothing}}; kwargs...) = AbstractCallback[]
to_cb_vector(cb::AbstractCallback; kwargs...) = [cb]
function to_cb_vector(cbs; CB_TYPE = SymbolicContinuousCallback, kwargs...)
    return if cbs isa Pair
        [CB_TYPE(cbs; kwargs...)]
    else
        cbs_filtered = filter!(c -> !isnothing(c), [CB_TYPE(cb; kwargs...) for cb in cbs])
        Vector{CB_TYPE}(cbs_filtered)
    end
end

############################################
########## Namespacing Utilities ###########
############################################
function namespace_affects(affect::AffectSystem, s)
    affsys = system(affect)
    old_ts = get_tearing_state(affsys)
    # if we just `renamespace` the system, it updates the name. However, this doesn't
    # namespace the returned values from `equations(affsys)`, etc. which we need. So we
    # need to manually namespace everything. This is done by renaming the system to the
    # namespace, putting it as a subsystem of an empty system called `affectsys`, and then
    # flatten the system. The resultant system has everything namespaced, and is still
    # called `affectsys` for further namespacing
    affsys = rename(affsys, nameof(s))
    affsys = toggle_namespacing(affsys, true)
    affsys = System(Equation[], get_iv(affsys); systems = [affsys], name = :affectsys)
    affsys = complete(affsys)
    @set! affsys.tearing_state = old_ts
    return AffectSystem(
        affsys,
        renamespace.((s,), unknowns(affect)),
        renamespace.((s,), parameters(affect)),
        renamespace.((s,), discretes(affect))
    )
end
function namespace_affects(affect::SymbolicAffect, s)
    return SymbolicAffect(
        namespace_equation.(affect.affect, (s,)),
        renamespace.((s,), affect.discrete_parameters)
    )
end

namespace_affects(af::Nothing, s) = nothing

function namespace_callback(cb::SymbolicContinuousCallback, s)::SymbolicContinuousCallback
    return SymbolicContinuousCallback(
        namespace_equation.(equations(cb), (s,)),
        namespace_affects(affects(cb), s),
        affect_neg = namespace_affects(affect_negs(cb), s),
        initialize = namespace_affects(initialize_affects(cb), s),
        finalize = namespace_affects(finalize_affects(cb), s),
        rootfind = cb.rootfind, reinitializealg = cb.reinitializealg,
        zero_crossing_id = cb.zero_crossing_id
    )
end

function namespace_conditions(condition, s)
    return is_timed_condition(condition) ? condition : namespace_expr(condition, s)
end

function namespace_callback(cb::SymbolicDiscreteCallback, s)::SymbolicDiscreteCallback
    return SymbolicDiscreteCallback(
        namespace_conditions(conditions(cb), s),
        namespace_affects(affects(cb), s),
        initialize = namespace_affects(initialize_affects(cb), s),
        finalize = namespace_affects(finalize_affects(cb), s), reinitializealg = cb.reinitializealg
    )
end

function Base.hash(cb::AbstractCallback, s::UInt)
    s = conditions(cb) isa AbstractVector ? foldr(hash, conditions(cb), init = s) :
        hash(conditions(cb), s)
    s = hash(affects(cb), s)
    !is_discrete(cb) && (s = hash(affect_negs(cb), s))
    s = hash(initialize_affects(cb), s)
    s = hash(finalize_affects(cb), s)
    !is_discrete(cb) && (s = hash(cb.rootfind, s))
    hash(cb.reinitializealg, s)
    !is_discrete(cb) && (s = hash(cb.zero_crossing_id, s))
    return s
end

###########################
######### Helpers #########
###########################

conditions(cb::AbstractCallback) = cb.conditions
function conditions(cbs::Vector{<:AbstractCallback})
    return reduce(vcat, conditions(cb) for cb in cbs; init = [])
end
function conditions(cbs::Vector{SymbolicContinuousCallback})
    return mapreduce(conditions, vcat, cbs; init = Equation[])
end
equations(cb::AbstractCallback) = conditions(cb)
equations(cb::Vector{<:AbstractCallback}) = conditions(cb)

affects(cb::AbstractCallback) = cb.affect
function affects(cbs::Vector{<:AbstractCallback})
    return reduce(vcat, affects(cb) for cb in cbs; init = [])
end

affect_negs(cb::SymbolicContinuousCallback) = cb.affect_neg
function affect_negs(cbs::Vector{SymbolicContinuousCallback})
    return reduce(vcat, affect_negs(cb) for cb in cbs; init = [])
end

initialize_affects(cb::AbstractCallback) = cb.initialize
function initialize_affects(cbs::Vector{<:AbstractCallback})
    return reduce(initialize_affects, vcat, cbs; init = [])
end

finalize_affects(cb::AbstractCallback) = cb.finalize
function finalize_affects(cbs::Vector{<:AbstractCallback})
    return reduce(finalize_affects, vcat, cbs; init = [])
end

function Base.:(==)(e1::AbstractCallback, e2::AbstractCallback)
    is_discrete(e1) === is_discrete(e2) || return false
    isequal(e1.conditions, e2.conditions) && isequal(e1.affect, e2.affect) || return false
    isequal(e1.initialize, e2.initialize) || return false
    isequal(e1.finalize, e2.finalize) || return false
    isequal(e1.reinitializealg, e2.reinitializealg) || return false
    if !is_discrete(e1)
        isequal(e1.affect_neg, e2.affect_neg) || return false
        isequal(e1.rootfind, e2.rootfind) || return false
        isequal(e1.zero_crossing_id, e2.zero_crossing_id) || return false
    end
    return true
end

Base.isempty(cb::AbstractCallback) = isempty(cb.conditions)

####################################
####### Compilation functions ######
####################################

struct CompiledCondition{IsDiscrete, F}
    f::F
end

function CompiledCondition{ID}(f::F) where {ID, F}
    return CompiledCondition{ID, F}(f)
end

function (cc::CompiledCondition)(out, u, t, integ)
    return cc.f(out, u, parameter_values(integ), t)
end

function (cc::CompiledCondition{false})(u, t, integ)
    return if DiffEqBase.isinplace(SciMLBase.get_sol(integ).prob)
        tmp, = DiffEqBase.get_tmp_cache(integ)
        cc.f(tmp, u, parameter_values(integ), t)
        tmp[1]
    else
        cc.f(u, parameter_values(integ), t)
    end
end

function (cc::CompiledCondition{true})(u, t, integ)
    return cc.f(u, parameter_values(integ), t)
end

####################################
#### Callable structs for affects ##
####################################

"""
    ExplicitAffect{DVS, PS, UF, PF}

Callable struct representing a compiled explicit affect (one with no algebraic equations).
Invokes `u_up!` to update state variables and `p_up!` to update discrete parameters, then
optionally resets aggregated jumps. Created by [`compile_explicit_affect`](@ref).

# Fields
- `dvs_to_update`: symbolic unknowns modified by this affect (emptiness is checked at call time)
- `ps_to_update`: symbolic discrete parameters modified by this affect
- `reset_jumps`: if `true`, call `reset_aggregated_jumps!` after the update
- `u_up!`: compiled in-place function that writes updated state into the integrator
- `p_up!`: compiled in-place function that writes updated parameters into the integrator
"""
struct ExplicitAffect{DVS, PS, UF, PF}
    dvs_to_update::DVS
    ps_to_update::PS
    reset_jumps::Bool
    u_up!::UF
    p_up!::PF
end

function (ea::ExplicitAffect)(integ)
    isempty(ea.dvs_to_update) || ea.u_up!(integ)
    isempty(ea.ps_to_update) || ea.p_up!(integ)
    return ea.reset_jumps && reset_aggregated_jumps!(integ)
end

"""
    ImplicitAffect{DVS, PS, AFFSYS, AFF, UG, AG, AUS, APS, US, PST, UGT, PG, PROB}

Callable struct representing a compiled implicit affect (one whose equations require solving
an `ImplicitDiscreteProblem` at each callback invocation). Created by
[`compile_implicit_affect`](@ref).

The `affprob` field is a plain (non-`Ref`) immutable field of type `PROB`. The in-place
setter functions `affu_setter!` and `affp_setter!` mutate the problem's internal mutable
arrays directly; `remake` is then called to produce a transient local copy with the current
`tspan` for the `init` call, without writing back to the struct.

# Fields
- `dvs_to_update`: unknowns of the parent system written back after solving
- `ps_to_update`: discrete parameters of the parent system written back after solving
- `affsys`: the compiled affect `System` (an `ImplicitDiscreteSystem`)
- `aff`: the `AffectSystem` descriptor (used for error reporting)
- `reset_jumps`: if `true`, call `reset_aggregated_jumps!` after the update
- `affu_getter`: reads current parent-system unknowns into the affect problem's `u0`
- `affp_getter`: reads current parent-system parameters into the affect problem's `p`
- `affu_setter!`: sets the affect problem's unknowns in-place
- `affp_setter!`: sets the affect problem's parameters in-place
- `u_setter!`: writes solved unknowns back into the parent integrator
- `p_setter!`: writes solved parameters back into the parent integrator
- `u_getter`: reads solved unknowns from the affect solution
- `p_getter`: reads solved parameters from the affect solution
- `affprob`: the pre-built `ImplicitDiscreteProblem` (mutated in-place each call)
"""
struct ImplicitAffect{DVS, PS, AFFSYS, AFF, UG, AG, AUS, APS, US, PST, UGT, PG, PROB}
    dvs_to_update::DVS
    ps_to_update::PS
    affsys::AFFSYS
    aff::AFF
    reset_jumps::Bool
    affu_getter::UG
    affp_getter::AG
    affu_setter!::AUS
    affp_setter!::APS
    u_setter!::US
    p_setter!::PST
    u_getter::UGT
    p_getter::PG
    affprob::PROB
end

function (ia::ImplicitAffect)(integ)
    ia.affu_setter!(ia.affprob, ia.affu_getter(integ))
    ia.affp_setter!(ia.affprob, ia.affp_getter(integ))
    # remake only updates tspan; result is a transient local, not stored back to struct
    affprob = remake(ia.affprob, tspan = (integ.t, integ.t))
    affsol = init(affprob, IDSolve())
    (check_error(affsol) === ReturnCode.InitialFailure) &&
        throw(UnsolvableCallbackError(all_equations(ia.aff)))
    ia.u_setter!(integ, ia.u_getter(affsol))
    ia.p_setter!(integ, ia.p_getter(affsol))
    return ia.reset_jumps && reset_aggregated_jumps!(integ)
end

"""
    VectorAffect{E2A, AFFS}

Callable struct for the positive-edge arm of a `VectorContinuousCallback`. Routes an
integrator call to the appropriate per-equation affect based on the equation index `idx`.
Created inside [`generate_callback`](@ref) for vectors of `SymbolicContinuousCallback`s.

# Fields
- `eq2affect`: maps condition equation index → affect index
- `affects`: vector of compiled affect callables, one per callback in the group
"""
struct VectorAffect{E2A, AFFS}
    eq2affect::E2A
    affects::AFFS
end

(va::VectorAffect)(integ, idx) = va.affects[va.eq2affect[idx]](integ)

"""
    VectorAffectNeg{E2A, AFFS}

Callable struct for the negative-edge arm of a `VectorContinuousCallback`. Like
[`VectorAffect`](@ref) but skips `nothing` entries (callbacks with no negative-edge affect).

# Fields
- `eq2affect`: maps condition equation index → affect index
- `affect_negs`: vector of compiled negative-edge affect callables (entries may be `nothing`)
"""
struct VectorAffectNeg{E2A, AFFS}
    eq2affect::E2A
    affect_negs::AFFS
end

function (va::VectorAffectNeg)(integ, idx)
    f = va.affect_negs[va.eq2affect[idx]]
    f === nothing && return
    return f(integ)
end

"""
    VectorOptionalAffect{FUNS}

Callable struct for the `initialize` or `finalize` slot of a `VectorContinuousCallback`.
Iterates over a heterogeneous vector of optional compiled affect callables, skipping
`nothing` entries. The four-argument signature `(cb, u, t, integ)` matches the SciMLBase
initialize/finalize protocol. Created by [`wrap_vector_optional_affect`](@ref).

# Fields
- `funs`: vector of compiled affect callables or `nothing`, one per callback in the group
"""
struct VectorOptionalAffect{FUNS}
    funs::FUNS
end

function (voa::VectorOptionalAffect)(cb, u, t, integ)
    for func in voa.funs
        func === nothing || func(integ)
    end
    return
end

"""
    InitFinalizeWrapper{F}

Callable struct that adapts a compiled single-integrator affect `f(integrator)` to the
four-argument `(cb, u, t, integrator)` protocol expected by SciMLBase's `initialize` and
`finalize` callback fields. Used in [`generate_callback`](@ref) for single callbacks whose
`initialize` or `finalize` field is non-`nothing`.

# Fields
- `f`: the inner compiled affect callable with signature `f(integrator)`
"""
struct InitFinalizeWrapper{F}
    f::F
end

(ifw::InitFinalizeWrapper)(c, u, t, integ) = ifw.f(integ)

"""
    compile_condition(cb::AbstractCallback, sys, dvs, ps; expression, kwargs...)

Returns a function `condition(u,t,integrator)`, condition(out,u,t,integrator)` returning the `condition(cb)`.
"""
Base.@nospecializeinfer function compile_condition(
        @nospecialize(cbs::Union{AbstractCallback, Vector{<:AbstractCallback}}),
        sys, @nospecialize(dvs), @nospecialize(ps);
        eval_expression = false, eval_module = @__MODULE__, kwargs...
    )
    u = map(value, dvs)
    p = map.(value, reorder_parameters(sys, ps))
    t = get_iv(sys)
    condit = conditions(cbs)

    if !is_discrete(cbs)
        condit = reduce(vcat, flatten_equations(Vector{Equation}(condit)))
        condit = condit isa AbstractVector ? [c.lhs - c.rhs for c in condit] :
            [condit.lhs - condit.rhs]
    end

    fs = build_function_wrapper(
        sys, condit, u, p..., t; kwargs..., cse = false
    )
    fs = GeneratedFunctionWrapper{(2, 3, is_split(sys))}(
        Val{false}, fs...; eval_expression, eval_module
    )
    return CompiledCondition{is_discrete(cbs)}(fs)
end

is_discrete(cb::AbstractCallback) = cb isa SymbolicDiscreteCallback
is_discrete(cb::Vector{<:AbstractCallback}) = eltype(cb) isa SymbolicDiscreteCallback

"""
    generate_continuous_callbacks(sys, dvs, ps; kwargs...)

Compile all continuous events of `sys` into SciMLBase callbacks. Groups the
`SymbolicContinuousCallback`s by their `(rootfind, reinitializealg, initialize_save_discretes)`
tuple and generates one `VectorContinuousCallback` (or `ContinuousCallback`) per group.

Returns a single callback, a `CallbackSet`, or `nothing` if the system has no continuous events.
"""
function generate_continuous_callbacks(
        sys::AbstractSystem, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true); kwargs...
    )
    cbs = continuous_events(sys)
    isempty(cbs) && return nothing
    cb_classes = Dict{
        Tuple{SciMLBase.RootfindOpt, SciMLBase.DAEInitializationAlgorithm, Bool},
        Vector{SymbolicContinuousCallback},
    }()

    # Sort the callbacks by their rootfinding method
    for cb in cbs
        _cbs = get!(
            () -> SymbolicContinuousCallback[],
            cb_classes, (cb.rootfind, cb.reinitializealg, cb.initialize_save_discretes)
        )
        push!(_cbs, cb)
    end
    sort!(OrderedDict(cb_classes), by = cb -> cb[1])
    # Use explicit loop to avoid Generator type inference overhead
    compiled_callbacks = Vector{Any}(undef, length(cb_classes))
    for (i, ((rf, reinit), cb)) in enumerate(cb_classes)
        compiled_callbacks[i] = generate_callback(cb, sys; kwargs...)
    end
    if length(compiled_callbacks) == 1
        return only(compiled_callbacks)
    else
        return CallbackSet(compiled_callbacks...)
    end
end

"""
    generate_discrete_callbacks(sys, dvs, ps; tspan, kwargs...)

Compile all discrete events of `sys` into a `Vector` of SciMLBase callbacks. Each
`SymbolicDiscreteCallback` becomes a `DiscreteCallback`, `PeriodicCallback`, or
`PresetTimeCallback` depending on its condition type.

Returns a non-empty `Vector`, or `nothing` if the system has no discrete events.
"""
function generate_discrete_callbacks(
        sys::AbstractSystem, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true); kwargs...
    )
    dbs = discrete_events(sys)
    isempty(dbs) && return nothing
    # Use explicit loop to avoid Generator type inference overhead
    result = Vector{Any}(undef, length(dbs))
    for (i, db) in enumerate(dbs)
        result[i] = generate_callback(db, sys; kwargs...)
    end
    return result
end

"""
    EMPTY_AFFECT(args...)

Sentinel affect that does nothing. Used as the default `affect` and `affect_neg` when a
callback specifies no affect, rather than passing `nothing` to SciMLBase.
"""
EMPTY_AFFECT(args...) = nothing

"""
    generate_callback(cbs::Vector{SymbolicContinuousCallback}, sys; kwargs...)

Generate a `VectorContinuousCallback` (or a single `ContinuousCallback` if there is only
one equation) from a homogeneous group of `SymbolicContinuousCallback`s that share a
rootfinding class. Delegates to the single-callback overload when `sum(num_eqs) == 1`.

Affect routing (from condition equation index to per-callback affect) is encoded in
[`VectorAffect`](@ref) and [`VectorAffectNeg`](@ref) callable structs.
Initialize/finalize are wrapped in [`VectorOptionalAffect`](@ref) via
[`wrap_vector_optional_affect`](@ref).
"""
function generate_callback(cbs::Vector{SymbolicContinuousCallback}, sys; kwargs...)
    eqs = map(cb -> flatten_equations(equations(cb)), cbs)
    num_eqs = length.(eqs)
    (isempty(eqs) || sum(num_eqs) == 0) && return nothing
    if sum(num_eqs) == 1
        cb_ind = findfirst(>(0), num_eqs)
        return generate_callback(cbs[cb_ind], sys; kwargs...)
    end

    if is_split(sys)
        ic = get_index_cache(sys)
    else
        ic = nothing
    end
    trigger = compile_condition(
        cbs, sys, unknowns(sys), parameters(sys; initial_parameters = true); kwargs...
    )
    compiled = compile_vector_callback_affects(cbs, sys, ic; kwargs...)

    # Build eq→affect index map: condition equation i maps to the affect for the callback
    # that owns it (a callback with k equations contributes k entries to this map).
    eq2affect = reduce(vcat, [fill(i, num_eqs[i]) for i in eachindex(compiled.affects)])
    eqs = reduce(vcat, eqs)

    affect     = VectorAffect(eq2affect, compiled.affects)
    affect_neg = VectorAffectNeg(eq2affect, compiled.affect_negs)
    initialize = wrap_vector_optional_affect(compiled.inits, SciMLBase.INITIALIZE_DEFAULT)
    finalize   = wrap_vector_optional_affect(compiled.finals, SciMLBase.FINALIZE_DEFAULT)

    return VectorContinuousCallback(
        trigger, affect, affect_neg, length(eqs); initialize, finalize,
        rootfind = cbs[1].rootfind, initializealg = cbs[1].reinitializealg,
        saved_clock_partitions = compiled.saved_clock_partitions,
        initialize_save_discretes = cbs[1].initialize_save_discretes
    )
end

"""
    compile_vector_callback_affects(cbs, sys, ic; kwargs...)

Compile the affect, affect_neg, initialize, and finalize callables for each callback in a
homogeneous group of `SymbolicContinuousCallback`s, and collect clock-partition save indices.

Returns a `NamedTuple` with fields:
- `affects`: compiled positive-edge affect callables
- `affect_negs`: compiled negative-edge affect callables (may contain `nothing`)
- `inits`: compiled initialize callables (may contain `nothing`)
- `finals`: compiled finalize callables (may contain `nothing`)
- `saved_clock_partitions`: `Vector{Vector{Int}}` of clock partition indices (empty for non-split systems)

# Arguments
- `cbs`: the `Vector{SymbolicContinuousCallback}` being compiled
- `sys`: the parent system providing compilation context
- `ic`: the index cache if `is_split(sys)`, else `nothing`
"""
function compile_vector_callback_affects(cbs, sys, ic; kwargs...)
    affects = []
    affect_negs = []
    inits = []
    finals = []
    saved_clock_partitions = Vector{Int}[]
    for cb in cbs
        affect = compile_affect(cb.affect, cb, sys; default = EMPTY_AFFECT, kwargs...)
        push!(affects, affect)
        affect_neg = (cb.affect_neg === cb.affect) ? affect :
            compile_affect(cb.affect_neg, cb, sys; default = EMPTY_AFFECT, kwargs...)
        push!(affect_negs, affect_neg)
        push!(inits, compile_affect(cb.initialize, cb, sys; default = nothing, kwargs...))
        push!(finals, compile_affect(cb.finalize, cb, sys; default = nothing, kwargs...))
        if ic !== nothing
            save_idxs = get(ic.callback_to_clocks, cb, Int[])
            for _ in conditions(cb)
                push!(saved_clock_partitions, save_idxs)
            end
        end
    end
    return (; affects, affect_negs, inits, finals, saved_clock_partitions)
end

"""
    generate_callback(cb::AbstractCallback, sys; tspan, kwargs...)

Generate a concrete SciMLBase callback from a single `SymbolicContinuousCallback` or
`SymbolicDiscreteCallback`. Dispatches on the condition type and discreteness:

- Timed + vector condition → `PresetTimeCallback`
- Timed + `SciMLBase.TimeDomain` condition → `PeriodicCallback`
- Timed + scalar → `PeriodicCallback`
- Untimed + discrete → `DiscreteCallback`
- Untimed + continuous → `ContinuousCallback`

Non-trivial `initialize`/`finalize` affects are wrapped in [`InitFinalizeWrapper`](@ref)
to adapt the `f(integrator)` signature to the `(cb, u, t, integrator)` SciMLBase protocol.

# Arguments
- `cb`: the symbolic callback descriptor
- `sys`: the parent system providing compilation context
- `tspan`: required for `PeriodicCallback` phase computation; may be `nothing` otherwise
"""
function generate_callback(cb, sys; tspan = nothing, kwargs...)
    is_timed = is_timed_condition(conditions(cb))
    dvs = unknowns(sys)
    ps = parameters(sys; initial_parameters = true)

    trigger = is_timed ? conditions(cb) : compile_condition(cb, sys, dvs, ps; kwargs...)
    affect = compile_affect(cb.affect, cb, sys; default = EMPTY_AFFECT, kwargs...)
    affect_neg = if is_discrete(cb)
        nothing
    else
        (cb.affect === cb.affect_neg) ? affect :
            compile_affect(cb.affect_neg, cb, sys; default = EMPTY_AFFECT, kwargs...)
    end
    init = compile_affect(
        cb.initialize, cb, sys; default = SciMLBase.INITIALIZE_DEFAULT,
        kwargs...
    )
    final = compile_affect(
        cb.finalize, cb, sys; default = SciMLBase.FINALIZE_DEFAULT, kwargs...
    )

    initialize = isnothing(cb.initialize) ? init : InitFinalizeWrapper(init)
    finalize   = isnothing(cb.finalize) ? final : InitFinalizeWrapper(final)

    saved_clock_partitions = if is_split(sys)
        get(get_index_cache(sys).callback_to_clocks, cb, ())
    else
        ()
    end
    if is_discrete(cb)
        if is_timed && conditions(cb) isa AbstractVector
            return PresetTimeCallback(
                trigger, affect; initialize,
                finalize, initializealg = cb.reinitializealg, saved_clock_partitions,
                initialize_save_discretes = cb.initialize_save_discretes
            )
        elseif is_timed && trigger isa SciMLBase.TimeDomain
            trigger_at_init = iszero((tspan[1] - trigger.phase) % trigger.dt)
            return PeriodicCallback(
                affect, trigger.dt; phase = trigger.phase, initial_affect = trigger_at_init,
                initialize, finalize,
                initializealg = cb.reinitializealg, saved_clock_partitions,
                initialize_save_discretes = trigger_at_init
            )
        elseif is_timed
            return PeriodicCallback(
                affect, trigger; initialize, finalize, initializealg = cb.reinitializealg,
                saved_clock_partitions, initialize_save_discretes = cb.initialize_save_discretes
            )
        else
            return DiscreteCallback(
                trigger, affect; initialize,
                finalize, initializealg = cb.reinitializealg, saved_clock_partitions,
                initialize_save_discretes = cb.initialize_save_discretes
            )
        end
    else
        return ContinuousCallback(
            trigger, affect, affect_neg; initialize, finalize,
            rootfind = cb.rootfind, initializealg = cb.reinitializealg, saved_clock_partitions,
            initialize_save_discretes = cb.initialize_save_discretes
        )
    end
end

"""
    compile_affect(aff, cb, sys; default, kwargs...)

Compile a single symbolic affect into a callable with signature `affect!(integrator)`.
Dispatches on the affect type:
- `nothing` → returns `default`
- `AffectSystem` → delegates to [`compile_equational_affect`](@ref)
- `ImperativeAffect` → delegates to [`compile_functional_affect`](@ref)
"""
function compile_affect(
        aff::Union{Nothing, Affect}, cb::AbstractCallback, sys::AbstractSystem;
        default = nothing, kwargs...
    )
    return if isnothing(aff)
        default
    elseif aff isa AffectSystem
        compile_equational_affect(aff, sys; kwargs...)
    elseif aff isa ImperativeAffect
        compile_functional_affect(aff, sys; kwargs...)
    end
end

"""
    wrap_vector_optional_affect(funs, default)

Create the `initialize` or `finalize` handler for a `VectorContinuousCallback`. If all
entries in `funs` are `nothing`, returns `default` (typically `SciMLBase.INITIALIZE_DEFAULT`
or `SciMLBase.FINALIZE_DEFAULT`). Otherwise returns a [`VectorOptionalAffect`](@ref) that
iterates over `funs` and calls each non-`nothing` entry with the integrator.
"""
function wrap_vector_optional_affect(funs, default)
    all(isnothing, funs) && return default
    return VectorOptionalAffect(funs)
end

"""
    add_integrator_header(sys, integrator, out) -> (oop_wrap, iip_wrap)

Return a pair of code-wrapping functions for use as `wrap_code` in `build_function_wrapper`.
Each wrapper transforms a generated expression by destructuring `integrator` into
`(:u, :p, :t)` for out-of-place (OOP) or `(out, :u, :p, :t)` for in-place (IIP), so the
resulting affect function is directly callable with an integrator rather than explicit
`(u, p, t)` arguments.
"""
function add_integrator_header(
        sys::AbstractSystem, integrator = gensym(:MTKIntegrator), out = :u
    )
    return expr -> Func(
            [DestructuredArgs(expr.args, integrator, inds = [:u, :p, :t])], [],
            expr.body
        ),
        expr -> Func(
            [DestructuredArgs(expr.args, integrator, inds = [out, :u, :p, :t])], [],
            expr.body
        )
end

"""
    default_operating_point(affsys::AffectSystem) -> AnyDict

Construct a default initial operating point for the `ImplicitDiscreteProblem` associated
with `affsys`. All unknowns are initialized to `0.0`; scalar numeric parameters to `false`;
fixed-size array parameters to `zeros(size(p))`. Used as the fallback when no `op` keyword
is passed to [`compile_equational_affect`](@ref).
"""
function default_operating_point(affsys::AffectSystem)
    sys = system(affsys)

    op = AnyDict(unknowns(sys) .=> 0.0)
    for p in parameters(sys)
        T = symtype(p)
        if T <: Number
            op[p] = false
        elseif T <: Array{<:Real} && symbolic_has_known_size(p)
            op[p] = zeros(size(p))
        end
    end
    return op
end

"""
    compile_equational_affect(aff, sys; reset_jumps, eval_expression, eval_module, op, kwargs...)

Compile an affect defined by a set of equations into a callable with signature
`affect!(integrator)`. Dispatches to [`compile_explicit_affect`](@ref) when the affect
system has no algebraic equations, or to [`compile_implicit_affect`](@ref) otherwise.

- Explicit affects generate in-place update functions via `build_function_wrapper` and
  return an [`ExplicitAffect`](@ref) callable struct.
- Implicit affects construct an `ImplicitDiscreteProblem` and return an
  [`ImplicitAffect`](@ref) callable struct that solves it at each callback invocation.

# Arguments
- `aff`: an `AffectSystem` or `Vector{Equation}` defining the affect
- `sys`: the parent `AbstractSystem`
- `reset_jumps`: if `true`, call `reset_aggregated_jumps!` after each invocation
- `op`: optional initial operating point for the implicit discrete problem
- `eval_expression`, `eval_module`: forwarded to code generation
"""
Base.@nospecializeinfer function compile_equational_affect(
        @nospecialize(aff::Union{AffectSystem, Vector{Equation}}), sys;
        reset_jumps = false, eval_expression = false, eval_module = @__MODULE__,
        @nospecialize(op = nothing), kwargs...
    )
    if aff isa AbstractVector
        aff = make_affect(aff; iv = get_iv(sys))
    end
    if op === nothing
        op = default_operating_point(aff)
    end
    if isempty(equations(system(aff)))
        return compile_explicit_affect(
            aff, sys; reset_jumps, eval_expression, eval_module, kwargs...
        )
    else
        return compile_implicit_affect(
            aff, sys; reset_jumps, eval_expression, eval_module, op, kwargs...
        )
    end
end

"""
    compile_explicit_affect(aff::AffectSystem, sys; reset_jumps, eval_expression, eval_module, kwargs...)

Compile an explicit (algebraic-equation-free) equational affect into an
[`ExplicitAffect`](@ref) callable struct. Generates in-place update functions for state
variables and discrete parameters via `build_function_wrapper`.

Called from [`compile_equational_affect`](@ref) when `isempty(equations(system(aff)))`.
"""
Base.@nospecializeinfer function compile_explicit_affect(
        @nospecialize(aff::AffectSystem), sys;
        reset_jumps = false, eval_expression = false, eval_module = @__MODULE__, kwargs...
    )
    affsys = system(aff)
    ps_to_update = discretes(aff)
    dvs_to_update = setdiff(unknowns(aff), getfield.(observed(sys), :lhs))

    _affsys = unhack_system(affsys)
    aff_ir_info = get_ir_info(affsys)
    aff_ir = get_irstructure(affsys)
    obseqs = Equation[]
    sizehint!(obseqs, length(observed(_affsys)))
    for (i, eq) in enumerate(observed(_affsys))
        push!(obseqs, eq.lhs ~ aff_ir_info.obs_subber(eq.rhs))
    end

    update_eqs = substitute(
        obseqs, Dict([p => unPre(p) for p in parameters(affsys)])
    )
    rhss = map(x -> x.rhs, update_eqs)
    lhss = map(x -> x.lhs, update_eqs)
    update_ps_set = Set(ps_to_update)
    is_p = map(lhss) do lhs
        lhs in update_ps_set ||
            iscall(lhs) && operation(lhs) === getindex &&
            arguments(lhs)[1] in update_ps_set
    end
    is_u = [lhs in Set(dvs_to_update) for lhs in lhss]
    dvs = unknowns(sys)
    ps = parameters(sys)
    t = get_iv(sys)

    u_idxs = indexin((@view lhss[is_u]), dvs)

    wrap_mtkparameters = has_index_cache(sys) && (get_index_cache(sys) !== nothing)
    p_idxs = if wrap_mtkparameters
        [parameter_index(sys, p) for (i, p) in enumerate(lhss) if is_p[i]]
    else
        indexin((@view lhss[is_p]), ps)
    end
    _ps = reorder_parameters(sys, ps)
    integ = gensym(:MTKIntegrator)

    u_up,
        u_up! = build_function_wrapper(
        sys, (@view rhss[is_u]), dvs, _ps..., t;
        wrap_code = add_integrator_header(sys, integ, :u),
        outputidxs = u_idxs, wrap_mtkparameters, iip_config = (false, true)
    )
    p_up,
        p_up! = build_function_wrapper(
        sys, (@view rhss[is_p]), dvs, _ps..., t;
        wrap_code = add_integrator_header(sys, integ, :p),
        outputidxs = p_idxs, wrap_mtkparameters, iip_config = (false, true)
    )

    u_up! = eval_or_rgf(u_up!; eval_expression, eval_module)
    p_up! = eval_or_rgf(p_up!; eval_expression, eval_module)

    return ExplicitAffect(dvs_to_update, ps_to_update, reset_jumps, u_up!, p_up!)
end

"""
    compile_implicit_affect(aff::AffectSystem, sys; reset_jumps, eval_expression, eval_module, op, kwargs...)

Compile an implicit (has algebraic equations) equational affect into an
[`ImplicitAffect`](@ref) callable struct. Constructs an `ImplicitDiscreteProblem` at
compile time; the struct's `affu_setter!` and `affp_setter!` mutate the problem's internal
arrays in-place at each callback invocation, and `remake` is used transiently to update the
`tspan` before calling `init`.

Called from [`compile_equational_affect`](@ref) when `!isempty(equations(system(aff)))`.
"""
Base.@nospecializeinfer function compile_implicit_affect(
        @nospecialize(aff::AffectSystem), sys;
        reset_jumps = false, eval_expression = false, eval_module = @__MODULE__,
        @nospecialize(op = nothing), kwargs...
    )
    affsys = system(aff)
    ps_to_update = discretes(aff)
    dvs_to_update = setdiff(unknowns(aff), getfield.(observed(sys), :lhs))

    dvs_to_access = unknowns(affsys)
    ps_to_access = [unPre(p) for p in parameters(affsys)]

    affu_getter  = getsym(sys, dvs_to_access)
    affp_getter  = getsym(sys, ps_to_access)
    affu_setter! = setsym(affsys, unknowns(affsys))
    affp_setter! = setsym(affsys, parameters(affsys))
    u_setter!    = setsym(sys, dvs_to_update)
    p_setter!    = setsym(sys, ps_to_update)
    u_getter     = getsym(affsys, dvs_to_update)
    p_getter     = getsym(affsys, ps_to_update)

    affprob = ImplicitDiscreteProblem(
        affsys, op, (0, 0);
        build_initializeprob = false, check_length = false, eval_expression,
        eval_module, check_compatibility = false, kwargs...
    )

    return ImplicitAffect(
        dvs_to_update, ps_to_update, affsys, aff,
        reset_jumps,
        affu_getter, affp_getter, affu_setter!, affp_setter!,
        u_setter!, p_setter!, u_getter, p_getter,
        affprob
    )
end

struct UnsolvableCallbackError
    eqs::Vector{Equation}
end

function Base.showerror(io::IO, err::UnsolvableCallbackError)
    return println(
        io,
        "The callback defined by the following equations:\n\n$(join(err.eqs, "\n"))\n\nis not solvable. Please check that the algebraic equations and affect equations are correct, and that all parameters intended to be changed are passed in as `discrete_parameters`."
    )
end

"""
    merge_cb(x, y) -> Union{Nothing, AbstractCallback, CallbackSet}

Merge two SciMLBase callbacks into a single `CallbackSet`. Handles `nothing` gracefully:
two `nothing`s return `nothing`; one `nothing` returns the other unchanged; two non-`nothing`
values return `CallbackSet(x, y)`.
"""
merge_cb(::Nothing, ::Nothing) = nothing
merge_cb(::Nothing, x) = merge_cb(x, nothing)
merge_cb(x, ::Nothing) = x
merge_cb(x, y) = CallbackSet(x, y)

"""
    process_events(sys; callback, tspan, kwargs...) -> Union{Nothing, AbstractCallback, CallbackSet}

Entry point for callback code generation when building an `ODEProblem` or `SDEProblem`
(called from `process_kwargs` in `problem_utils.jl`). Compiles all continuous and discrete
symbolic events attached to `sys` into concrete SciMLBase callbacks, then merges them with
any user-supplied `callback`.

Returns `nothing` if `sys` has no events and `callback` is `nothing`.

# Arguments
- `sys`: a compiled `AbstractSystem` with attached symbolic events
- `callback`: an optional pre-built SciMLBase callback to merge into the result
- `tspan`: the integration time span; required for `PeriodicCallback` phase computation
- `kwargs...`: forwarded to [`generate_continuous_callbacks`](@ref) and
  [`generate_discrete_callbacks`](@ref)

# See also
[`generate_continuous_callbacks`](@ref), [`generate_discrete_callbacks`](@ref)
"""
function process_events(sys; callback = nothing, tspan = nothing, kwargs...)
    contin_cbs = generate_continuous_callbacks(sys; kwargs...)
    discrete_cbs = generate_discrete_callbacks(sys; tspan, kwargs...)
    cb = merge_cb(contin_cbs, callback)
    return (discrete_cbs === nothing) ? cb : CallbackSet(contin_cbs, discrete_cbs...)
end

"""
    discrete_events(sys::AbstractSystem) :: Vector{SymbolicDiscreteCallback}

Returns a vector of all the `discrete_events` in an abstract system and its component subsystems.
The `SymbolicDiscreteCallback`s in the returned vector are structs with two fields: `condition` and
`affect` which correspond to the first and second elements of a `Pair` used to define an event, i.e.
`condition => affect`.

See also `get_discrete_events`, which only returns the events of the top-level system.
"""
function discrete_events(sys::AbstractSystem)
    cbs = get_discrete_events(sys)
    systems = get_systems(sys)
    cbs = copy(cbs)
    for s in systems
        append!(cbs, map(Base.Fix2(namespace_callback, s), discrete_events(s)))
    end
    return cbs
end

"""
    $(TYPEDSIGNATURES)

Returns whether the system `sys` has the internal field `discrete_events`.

See also [`get_discrete_events`](@ref).
"""
has_discrete_events(sys::AbstractSystem) = isdefined(sys, :discrete_events)
"""
    $(TYPEDSIGNATURES)

Get the internal field `discrete_events` of a system `sys`.
It only includes `discrete_events` local to `sys`; not those of its subsystems,
like `unknowns(sys)`, `parameters(sys)` and `equations(sys)` does.

See also [`has_discrete_events`](@ref).
"""
function get_discrete_events(sys::AbstractSystem)
    has_discrete_events(sys) || return SymbolicDiscreteCallback[]
    return getfield(sys, :discrete_events)
end

"""
    discrete_events_toplevel(sys::AbstractSystem)

Replicates the behaviour of `discrete_events`, but ignores events of subsystems.

Notes:
- Cannot be applied to non-complete systems.
"""
function discrete_events_toplevel(sys::AbstractSystem)
    if has_parent(sys) && (parent = get_parent(sys)) !== nothing
        return discrete_events_toplevel(parent)
    end
    return get_discrete_events(sys)
end

"""
    continuous_events(sys::AbstractSystem)::Vector{SymbolicContinuousCallback}

Returns a vector of all the `continuous_events` in an abstract system and its component subsystems.
The `SymbolicContinuousCallback`s in the returned vector are structs with two fields: `eqs` and
`affect` which correspond to the first and second elements of a `Pair` used to define an event, i.e.
`eqs => affect`.

See also `get_continuous_events`, which only returns the events of the top-level system.
"""
function continuous_events(sys::AbstractSystem)
    cbs = get_continuous_events(sys)
    systems = get_systems(sys)
    cbs = copy(cbs)
    for s in systems
        append!(cbs, map(Base.Fix2(namespace_callback, s), continuous_events(s)))
    end
    return cbs
end

"""
    $(TYPEDSIGNATURES)

Returns whether the system `sys` has the internal field `continuous_events`.

See also [`get_continuous_events`](@ref).
"""
has_continuous_events(sys::AbstractSystem) = isdefined(sys, :continuous_events)
"""
    $(TYPEDSIGNATURES)

Get the internal field `continuous_events` of a system `sys`.
It only includes `continuous_events` local to `sys`; not those of its subsystems,
like `unknowns(sys)`, `parameters(sys)` and `equations(sys)` does.

See also [`has_continuous_events`](@ref).
"""
function get_continuous_events(sys::AbstractSystem)
    has_continuous_events(sys) || return SymbolicContinuousCallback[]
    return getfield(sys, :continuous_events)
end

"""
    continuous_events_toplevel(sys::AbstractSystem)

Replicates the behaviour of `continuous_events`, but ignores events of subsystems.

Notes:
- Cannot be applied to non-complete systems.
"""
function continuous_events_toplevel(sys::AbstractSystem)
    if has_parent(sys) && (parent = get_parent(sys)) !== nothing
        return continuous_events_toplevel(parent)
    end
    return get_continuous_events(sys)
end

"""
Process the symbolic events of a system.
"""
function create_symbolic_events(cont_events, disc_events)
    cont_callbacks = to_cb_vector(cont_events; CB_TYPE = SymbolicContinuousCallback)
    disc_callbacks = to_cb_vector(disc_events; CB_TYPE = SymbolicDiscreteCallback)
    return cont_callbacks, disc_callbacks
end
