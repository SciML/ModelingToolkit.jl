abstract type AbstractCallback end

function has_functional_affect(cb)
    affects(cb) isa ImperativeAffect
end

struct SymbolicAffect
    affect::Vector{Equation}
    alg_eqs::Vector{Equation}
    discrete_parameters::Vector{Any}
end

function SymbolicAffect(affect::Vector{Equation}; alg_eqs = Equation[],
        discrete_parameters = Any[], kwargs...)
    if !(discrete_parameters isa AbstractVector)
        discrete_parameters = Any[discrete_parameters]
    elseif !(discrete_parameters isa Vector{Any})
        discrete_parameters = Vector{Any}(discrete_parameters)
    end
    SymbolicAffect(affect, alg_eqs, discrete_parameters)
end
function SymbolicAffect(affect::SymbolicAffect; kwargs...)
    SymbolicAffect(affect.affect; alg_eqs = affect.alg_eqs,
        discrete_parameters = affect.discrete_parameters, kwargs...)
end
SymbolicAffect(affect; kwargs...) = make_affect(affect; kwargs...)

function Symbolics.fast_substitute(aff::SymbolicAffect, rules)
    substituter = Base.Fix2(fast_substitute, rules)
    SymbolicAffect(map(substituter, aff.affect), map(substituter, aff.alg_eqs),
        map(substituter, aff.discrete_parameters))
end

struct AffectSystem
    """The internal implicit discrete system whose equations are solved to obtain values after the affect."""
    system::AbstractSystem
    """Unknowns of the parent ODESystem whose values are modified or accessed by the affect."""
    unknowns::Vector
    """Parameters of the parent ODESystem whose values are accessed by the affect."""
    parameters::Vector
    """Parameters of the parent ODESystem whose values are modified by the affect."""
    discretes::Vector
end

function Symbolics.fast_substitute(aff::AffectSystem, rules)
    substituter = Base.Fix2(fast_substitute, rules)
    sys = aff.system
    @set! sys.eqs = map(substituter, get_eqs(sys))
    @set! sys.parameter_dependencies = map(substituter, get_parameter_dependencies(sys))
    @set! sys.defaults = Dict([k => substituter(v) for (k, v) in defaults(sys)])
    @set! sys.guesses = Dict([k => substituter(v) for (k, v) in guesses(sys)])
    @set! sys.unknowns = map(substituter, get_unknowns(sys))
    @set! sys.ps = map(substituter, get_ps(sys))
    AffectSystem(sys, map(substituter, aff.unknowns),
        map(substituter, aff.parameters), map(substituter, aff.discretes))
end

function AffectSystem(spec::SymbolicAffect; iv = nothing, alg_eqs = Equation[], kwargs...)
    AffectSystem(spec.affect; alg_eqs = vcat(spec.alg_eqs, alg_eqs), iv,
        discrete_parameters = spec.discrete_parameters, kwargs...)
end

function AffectSystem(affect::Vector{Equation}; discrete_parameters = Any[],
        iv = nothing, alg_eqs::Vector{Equation} = Equation[], warn_no_algebraic = true, kwargs...)
    isempty(affect) && return nothing
    if isnothing(iv)
        iv = t_nounits
        @warn "No independent variable specified. Defaulting to t_nounits."
    end

    discrete_parameters isa AbstractVector || (discrete_parameters = [discrete_parameters])
    discrete_parameters = unwrap.(discrete_parameters)

    for p in discrete_parameters
        occursin(unwrap(iv), unwrap(p)) ||
            error("Non-time dependent parameter $p passed in as a discrete. Must be declared as @parameters $p(t).")
    end

    dvs = OrderedSet()
    params = OrderedSet()
    _varsbuf = Set()
    for eq in affect
        if !haspre(eq) && !(symbolic_type(eq.rhs) === NotSymbolic() ||
             symbolic_type(eq.lhs) === NotSymbolic())
            @warn "Affect equation $eq has no `Pre` operator. As such it will be interpreted as an algebraic equation to be satisfied after the callback. If you intended to use the value of a variable x before the affect, use Pre(x). Errors may be thrown if there is no `Pre` and the algebraic equation is unsatisfiable, such as X ~ X + 1."
        end
        collect_vars!(dvs, params, eq, iv; op = Pre)
        empty!(_varsbuf)
        vars!(_varsbuf, eq; op = Pre)
        filter!(x -> iscall(x) && operation(x) isa Pre, _varsbuf)
        union!(params, _varsbuf)
        diffvs = collect_applied_operators(eq, Differential)
        union!(dvs, diffvs)
    end
    for eq in alg_eqs
        collect_vars!(dvs, params, eq, iv)
    end
    pre_params = filter(haspre ∘ value, params)
    discrete_parameters = gather_array_params(OrderedSet(discrete_parameters))
    sys_params = collect(setdiff(params, union(discrete_parameters, pre_params)))
    discrete_parameters = collect(discrete_parameters)
    discretes = map(tovar, discrete_parameters)
    dvs = collect(dvs)
    _dvs = map(default_toterm, dvs)

    rev_map = Dict(zip(discrete_parameters, discretes))
    subs = merge(rev_map, Dict(zip(dvs, _dvs)))
    affect = Symbolics.fast_substitute(affect, subs)
    alg_eqs = Symbolics.fast_substitute(alg_eqs, subs)

    @named affectsys = System(
        vcat(affect, alg_eqs), iv, collect(union(_dvs, discretes)),
        collect(union(pre_params, sys_params)); is_discrete = true)
    affectsys = mtkcompile(affectsys; fully_determined = nothing)
    # get accessed parameters p from Pre(p) in the callback parameters
    accessed_params = Vector{Any}(filter(isparameter, map(unPre, collect(pre_params))))
    union!(accessed_params, sys_params)

    # add scalarized unknowns to the map.
    _dvs = reduce(vcat, map(scalarize, _dvs), init = Any[])

    AffectSystem(affectsys, collect(_dvs), collect(accessed_params),
        collect(discrete_parameters))
end

system(a::AffectSystem) = a.system
discretes(a::AffectSystem) = a.discretes
unknowns(a::AffectSystem) = a.unknowns
parameters(a::AffectSystem) = a.parameters
all_equations(a::AffectSystem) = vcat(equations(system(a)), observed(system(a)))

function Base.show(iio::IO, aff::AffectSystem)
    println(iio, "Affect system defined by equations:")
    eqs = all_equations(aff)
    show(iio, eqs)
end

function Base.:(==)(a1::AffectSystem, a2::AffectSystem)
    isequal(system(a1), system(a2)) &&
        isequal(discretes(a1), discretes(a2)) &&
        isequal(unknowns(a1), unknowns(a2)) &&
        isequal(parameters(a1), parameters(a2))
end

function Base.hash(a::AffectSystem, s::UInt)
    s = hash(system(a), s)
    s = hash(unknowns(a), s)
    s = hash(parameters(a), s)
    hash(discretes(a), s)
end

function vars!(vars, aff::AffectSystem; op = Differential)
    for var in Iterators.flatten((unknowns(aff), parameters(aff), discretes(aff)))
        vars!(vars, var)
    end
    vars
end

"""
    Pre(x)

The `Pre` operator. Used by the callback system to indicate the value of a parameter or variable
before the callback is triggered.
"""
struct Pre <: Symbolics.Operator end
Pre(x) = Pre()(x)
is_timevarying_operator(::Type{Pre}) = false
SymbolicUtils.promote_symtype(::Type{Pre}, T) = T
SymbolicUtils.isbinop(::Pre) = false
Base.nameof(::Pre) = :Pre
Base.show(io::IO, x::Pre) = print(io, "Pre")
unPre(x::Num) = unPre(unwrap(x))
unPre(x::Symbolics.Arr) = unPre(unwrap(x))
unPre(x::Symbolic) = (iscall(x) && operation(x) isa Pre) ? only(arguments(x)) : x

function (p::Pre)(x)
    iw = Symbolics.iswrapped(x)
    x = unwrap(x)
    # non-symbolic values don't change
    if symbolic_type(x) == NotSymbolic()
        return x
    end
    # differential variables are default-toterm-ed
    if iscall(x) && operation(x) isa Differential
        x = default_toterm(x)
    end
    # don't double wrap
    iscall(x) && operation(x) isa Pre && return x
    result = if symbolic_type(x) == ArraySymbolic()
        # create an array for `Pre(array)`
        Symbolics.array_term(p, x)
    elseif iscall(x) && operation(x) == getindex
        # instead of `Pre(x[1])` create `Pre(x)[1]`
        # which allows parameter indexing to handle this case automatically.
        arr = arguments(x)[1]
        term(getindex, p(arr), arguments(x)[2:end]...)
    else
        term(p, x)
    end
    # the result should be a parameter
    result = toparam(result)
    if iw
        result = wrap(result)
    end
    return result
end
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
                               affect_neg = affect, initialize = nothing, finalize = nothing, rootfind = SciMLBase.LeftRootFind, alg_eqs = Equation[])

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
* A list of equations that should be applied when the callback is triggered (e.g. `x ~ 3, y ~ 7`) which must be of the form `unknown ~ observed value` where each `unknown` appears only once. Equations will be applied in the order that they appear in the vector; parameters and state updates will become immediately visible to following equations.
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
        kwargs...)
    conditions = (conditions isa AbstractVector) ? conditions : [conditions]

    if isnothing(reinitializealg)
        if any(a -> a isa ImperativeAffect,
            [affect, affect_neg, initialize, finalize])
            reinitializealg = SciMLBase.CheckInit()
        else
            reinitializealg = SciMLBase.NoInit()
        end
    end

    SymbolicContinuousCallback(conditions, SymbolicAffect(affect; kwargs...),
        SymbolicAffect(affect_neg; kwargs...),
        SymbolicAffect(initialize; kwargs...), SymbolicAffect(
            finalize; kwargs...),
        rootfind, reinitializealg, zero_crossing_id)
end # Default affect to nothing

function SymbolicContinuousCallback(p::Pair, args...; kwargs...)
    SymbolicContinuousCallback(p[1], p[2], args...; kwargs...)
end
SymbolicContinuousCallback(cb::SymbolicContinuousCallback, args...; kwargs...) = cb
SymbolicContinuousCallback(cb::Nothing, args...; kwargs...) = nothing
function SymbolicContinuousCallback(cb::Tuple, args...; kwargs...)
    if length(cb) == 2
        SymbolicContinuousCallback(cb[1]; kwargs..., cb[2]...)
    else
        error("Malformed tuple specifying callback. Should be a condition => affect pair, followed by a vector of kwargs.")
    end
end

function complete(cb::SymbolicContinuousCallback; kwargs...)
    SymbolicContinuousCallback(cb.conditions, make_affect(cb.affect; kwargs...),
        make_affect(cb.affect_neg; kwargs...), make_affect(cb.initialize; kwargs...),
        make_affect(cb.finalize; kwargs...), cb.rootfind, cb.reinitializealg, cb.zero_crossing_id)
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
    print(iio, ")")
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
    if finalize_affects(cb) != nothing
        println(iio, "Finalization affect:")
        show(iio, mime, finalize_affects(cb))
        print(iio, "\n")
    end
end

function vars!(vars, cb::AbstractCallback; op = Differential)
    if symbolic_type(conditions(cb)) == NotSymbolic
        if conditions(cb) isa AbstractArray
            for eq in conditions(cb)
                vars!(vars, eq; op)
            end
        end
    else
        vars!(vars, conditions(cb); op)
    end
    for aff in (affects(cb), initialize_affects(cb), finalize_affects(cb))
        isnothing(aff) || vars!(vars, aff; op)
    end
    !is_discrete(cb) && vars!(vars, affect_negs(cb); op)
    return vars
end

################################
######## Discrete events #######
################################
"""
    SymbolicDiscreteCallback(conditions::Vector{Equation}, affect = nothing, iv = nothing;
                             initialize = nothing, finalize = nothing, alg_eqs = Equation[])

A callback that triggers at the first timestep that the conditions are satisfied.

The condition can be one of: 
- Δt::Real              - periodic events with period Δt
- ts::Vector{Real}      - events trigger at these preset times given by `ts`
- eqs::Vector{Symbolic} - events trigger when the condition evaluates to true

Arguments: 
- iv: The independent variable of the system. This must be specified if the independent variable appears in one of the equations explicitly, as in x ~ t + 1.
- alg_eqs: Algebraic equations of the system that must be satisfied after the callback occurs.
"""
struct SymbolicDiscreteCallback <: AbstractCallback
    conditions::Union{Number, Vector{<:Number}, Symbolic{Bool}}
    affect::Union{Affect, SymbolicAffect, Nothing}
    initialize::Union{Affect, SymbolicAffect, Nothing}
    finalize::Union{Affect, SymbolicAffect, Nothing}
    reinitializealg::SciMLBase.DAEInitializationAlgorithm
end

function SymbolicDiscreteCallback(
        condition::Union{Symbolic{Bool}, Number, Vector{<:Number}}, affect = nothing;
        initialize = nothing, finalize = nothing,
        reinitializealg = nothing, kwargs...)
    # Manual error check (to prevent events like `[X < 5.0] => [X ~ Pre(X) + 10.0]` from being created).
    (condition isa Vector) && (eltype(condition) <: Num) &&
        error("Vectors of symbolic conditions are not allowed for `SymbolicDiscreteCallback`.")

    c = is_timed_condition(condition) ? condition : value(scalarize(condition))
    if isnothing(reinitializealg)
        if any(a -> a isa ImperativeAffect,
            [affect, initialize, finalize])
            reinitializealg = SciMLBase.CheckInit()
        else
            reinitializealg = SciMLBase.NoInit()
        end
    end
    SymbolicDiscreteCallback(c, SymbolicAffect(affect; kwargs...),
        SymbolicAffect(initialize; kwargs...),
        SymbolicAffect(finalize; kwargs...), reinitializealg)
end # Default affect to nothing

function SymbolicDiscreteCallback(p::Pair, args...; kwargs...)
    SymbolicDiscreteCallback(p[1], p[2], args...; kwargs...)
end
SymbolicDiscreteCallback(cb::SymbolicDiscreteCallback, args...; kwargs...) = cb
SymbolicDiscreteCallback(cb::Nothing, args...; kwargs...) = nothing
function SymbolicDiscreteCallback(cb::Tuple, args...; kwargs...)
    if length(cb) == 2
        SymbolicDiscreteCallback(cb[1]; cb[2]...)
    else
        error("Malformed tuple specifying callback. Should be a condition => affect pair, followed by a vector of kwargs.")
    end
end

function complete(cb::SymbolicDiscreteCallback; kwargs...)
    SymbolicDiscreteCallback(cb.conditions, make_affect(cb.affect; kwargs...),
        make_affect(cb.initialize; kwargs...),
        make_affect(cb.finalize; kwargs...), cb.reinitializealg)
end

function is_timed_condition(condition::T) where {T}
    if T === Num
        false
    elseif T <: Real
        true
    elseif T <: AbstractVector
        eltype(condition) <: Real
    else
        false
    end
end

to_cb_vector(cbs::Vector{<:AbstractCallback}; kwargs...) = cbs
to_cb_vector(cbs::Union{Nothing, Vector{Nothing}}; kwargs...) = AbstractCallback[]
to_cb_vector(cb::AbstractCallback; kwargs...) = [cb]
function to_cb_vector(cbs; CB_TYPE = SymbolicContinuousCallback, kwargs...)
    if cbs isa Pair
        [CB_TYPE(cbs; kwargs...)]
    else
        Vector{CB_TYPE}([CB_TYPE(cb; kwargs...) for cb in cbs])
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
    AffectSystem(affsys,
        renamespace.((s,), unknowns(affect)),
        renamespace.((s,), parameters(affect)),
        renamespace.((s,), discretes(affect)))
end
function namespace_affects(affect::SymbolicAffect, s)
    SymbolicAffect(
        namespace_equation.(affect.affect, (s,)), namespace_equation.(affect.alg_eqs, (s,)),
        renamespace.((s,), affect.discrete_parameters))
end

namespace_affects(af::Nothing, s) = nothing

function namespace_callback(cb::SymbolicContinuousCallback, s)::SymbolicContinuousCallback
    SymbolicContinuousCallback(
        namespace_equation.(equations(cb), (s,)),
        namespace_affects(affects(cb), s),
        affect_neg = namespace_affects(affect_negs(cb), s),
        initialize = namespace_affects(initialize_affects(cb), s),
        finalize = namespace_affects(finalize_affects(cb), s),
        rootfind = cb.rootfind, reinitializealg = cb.reinitializealg,
        zero_crossing_id = cb.zero_crossing_id)
end

function namespace_conditions(condition, s)
    is_timed_condition(condition) ? condition : namespace_expr(condition, s)
end

function namespace_callback(cb::SymbolicDiscreteCallback, s)::SymbolicDiscreteCallback
    SymbolicDiscreteCallback(
        namespace_conditions(conditions(cb), s),
        namespace_affects(affects(cb), s),
        initialize = namespace_affects(initialize_affects(cb), s),
        finalize = namespace_affects(finalize_affects(cb), s), reinitializealg = cb.reinitializealg)
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
    reduce(vcat, conditions(cb) for cb in cbs; init = [])
end
equations(cb::AbstractCallback) = conditions(cb)
equations(cb::Vector{<:AbstractCallback}) = conditions(cb)

affects(cb::AbstractCallback) = cb.affect
function affects(cbs::Vector{<:AbstractCallback})
    reduce(vcat, affects(cb) for cb in cbs; init = [])
end

affect_negs(cb::SymbolicContinuousCallback) = cb.affect_neg
function affect_negs(cbs::Vector{SymbolicContinuousCallback})
    reduce(vcat, affect_negs(cb) for cb in cbs; init = [])
end

initialize_affects(cb::AbstractCallback) = cb.initialize
function initialize_affects(cbs::Vector{<:AbstractCallback})
    reduce(initialize_affects, vcat, cbs; init = [])
end

finalize_affects(cb::AbstractCallback) = cb.finalize
function finalize_affects(cbs::Vector{<:AbstractCallback})
    reduce(finalize_affects, vcat, cbs; init = [])
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
    cc.f(out, u, parameter_values(integ), t)
end

function (cc::CompiledCondition{false})(u, t, integ)
    if DiffEqBase.isinplace(SciMLBase.get_sol(integ).prob)
        tmp, = DiffEqBase.get_tmp_cache(integ)
        cc.f(tmp, u, parameter_values(integ), t)
        tmp[1]
    else
        cc.f(u, parameter_values(integ), t)
    end
end

function (cc::CompiledCondition{true})(u, t, integ)
    cc.f(u, parameter_values(integ), t)
end

"""
    compile_condition(cb::AbstractCallback, sys, dvs, ps; expression, kwargs...)

Returns a function `condition(u,t,integrator)`, condition(out,u,t,integrator)` returning the `condition(cb)`.
"""
function compile_condition(
        cbs::Union{AbstractCallback, Vector{<:AbstractCallback}}, sys, dvs, ps;
        eval_expression = false, eval_module = @__MODULE__, kwargs...)
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
        sys, condit, u, p..., t; kwargs..., cse = false)
    if is_discrete(cbs)
        fs = (fs, nothing)
    end
    fs = GeneratedFunctionWrapper{(2, 3, is_split(sys))}(
        Val{false}, fs...; eval_expression, eval_module)
    return CompiledCondition{is_discrete(cbs)}(fs)
end

is_discrete(cb::AbstractCallback) = cb isa SymbolicDiscreteCallback
is_discrete(cb::Vector{<:AbstractCallback}) = eltype(cb) isa SymbolicDiscreteCallback

function generate_continuous_callbacks(sys::AbstractSystem, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true); kwargs...)
    cbs = continuous_events(sys)
    isempty(cbs) && return nothing
    cb_classes = Dict{Tuple{SciMLBase.RootfindOpt, SciMLBase.DAEInitializationAlgorithm},
        Vector{SymbolicContinuousCallback}}()

    # Sort the callbacks by their rootfinding method
    for cb in cbs
        _cbs = get!(() -> SymbolicContinuousCallback[],
            cb_classes, (cb.rootfind, cb.reinitializealg))
        push!(_cbs, cb)
    end
    sort!(OrderedDict(cb_classes), by = cb -> cb[1])
    compiled_callbacks = [generate_callback(cb, sys; kwargs...)
                          for ((rf, reinit), cb) in cb_classes]
    if length(compiled_callbacks) == 1
        return only(compiled_callbacks)
    else
        return CallbackSet(compiled_callbacks...)
    end
end

function generate_discrete_callbacks(sys::AbstractSystem, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true); kwargs...)
    dbs = discrete_events(sys)
    isempty(dbs) && return nothing
    [generate_callback(db, sys; kwargs...) for db in dbs]
end

EMPTY_AFFECT(args...) = nothing

"""
Codegen a DifferentialEquations callback. A (set of) continuous callback with multiple equations becomes a VectorContinuousCallback.
Continuous callbacks with only one equation will become a ContinuousCallback.
Individual discrete callbacks become DiscreteCallback, PresetTimeCallback, PeriodicCallback depending on the case.
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
        cbs, sys, unknowns(sys), parameters(sys; initial_parameters = true); kwargs...)
    affects = []
    affect_negs = []
    inits = []
    finals = []
    saved_clock_partitions = Vector{Int}[]
    for cb in cbs
        affect = compile_affect(cb.affect, cb, sys; default = EMPTY_AFFECT, kwargs...)
        push!(affects, affect)
        affect_neg = (cb.affect_neg === cb.affect) ? affect :
                     compile_affect(
            cb.affect_neg, cb, sys; default = EMPTY_AFFECT, kwargs...)
        push!(affect_negs, affect_neg)
        push!(inits,
            compile_affect(
                cb.initialize, cb, sys; default = nothing, kwargs...))
        push!(finals, compile_affect(cb.finalize, cb, sys; default = nothing, kwargs...))

        if ic !== nothing
            save_idxs = get(ic.callback_to_clocks, cb, Int[])
            for _ in conditions(cb)
                push!(saved_clock_partitions, save_idxs)
            end
        end
    end

    # Since there may be different number of conditions and affects,
    # we build a map that translates the condition eq. number to the affect number
    eq2affect = reduce(vcat,
        [fill(i, num_eqs[i]) for i in eachindex(affects)])
    eqs = reduce(vcat, eqs)

    affect = let eq2affect = eq2affect, affects = affects
        function (integ, idx)
            affects[eq2affect[idx]](integ)
        end
    end
    affect_neg = let eq2affect = eq2affect, affect_negs = affect_negs
        function (integ, idx)
            f = affect_negs[eq2affect[idx]]
            isnothing(f) && return
            f(integ)
        end
    end
    initialize = wrap_vector_optional_affect(inits, SciMLBase.INITIALIZE_DEFAULT)
    finalize = wrap_vector_optional_affect(finals, SciMLBase.FINALIZE_DEFAULT)

    return VectorContinuousCallback(
        trigger, affect, affect_neg, length(eqs); initialize, finalize,
        rootfind = cbs[1].rootfind, initializealg = cbs[1].reinitializealg,
        saved_clock_partitions)
end

function generate_callback(cb, sys; kwargs...)
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
    init = compile_affect(cb.initialize, cb, sys; default = SciMLBase.INITIALIZE_DEFAULT,
        kwargs...)
    final = compile_affect(
        cb.finalize, cb, sys; default = SciMLBase.FINALIZE_DEFAULT, kwargs...)

    initialize = isnothing(cb.initialize) ? init : ((c, u, t, i) -> init(i))
    finalize = isnothing(cb.finalize) ? final : ((c, u, t, i) -> final(i))

    saved_clock_partitions = if is_split(sys)
        get(get_index_cache(sys).callback_to_clocks, cb, ())
    else
        ()
    end
    if is_discrete(cb)
        if is_timed && conditions(cb) isa AbstractVector
            return PresetTimeCallback(trigger, affect; initialize,
                finalize, initializealg = cb.reinitializealg, saved_clock_partitions)
        elseif is_timed
            return PeriodicCallback(
                affect, trigger; initialize, finalize, initializealg = cb.reinitializealg,
                saved_clock_partitions)
        else
            return DiscreteCallback(trigger, affect; initialize,
                finalize, initializealg = cb.reinitializealg, saved_clock_partitions)
        end
    else
        return ContinuousCallback(trigger, affect, affect_neg; initialize, finalize,
            rootfind = cb.rootfind, initializealg = cb.reinitializealg, saved_clock_partitions)
    end
end

"""
    compile_affect(cb::AbstractCallback, sys::AbstractSystem, dvs, ps; expression, outputidxs, kwargs...)

Returns a function that takes an integrator as argument and modifies the state with the
affect. The generated function has the signature `affect!(integrator)`.

Notes
  - `kwargs` are passed through to `Symbolics.build_function`.
"""
function compile_affect(
        aff::Union{Nothing, Affect}, cb::AbstractCallback, sys::AbstractSystem;
        default = nothing, kwargs...)
    if isnothing(aff)
        default
    elseif aff isa AffectSystem
        compile_equational_affect(aff, sys; kwargs...)
    elseif aff isa ImperativeAffect
        compile_functional_affect(aff, sys; kwargs...)
    end
end

"""
Initialize and finalize for VectorContinuousCallback.
"""
function wrap_vector_optional_affect(funs, default)
    all(isnothing, funs) && return default
    return let funs = funs
        function (cb, u, t, integ)
            for func in funs
                isnothing(func) ? continue : func(integ)
            end
        end
    end
end

function add_integrator_header(
        sys::AbstractSystem, integrator = gensym(:MTKIntegrator), out = :u)
    expr -> Func([DestructuredArgs(expr.args, integrator, inds = [:u, :p, :t])], [],
        expr.body),
    expr -> Func(
        [DestructuredArgs(expr.args, integrator, inds = [out, :u, :p, :t])], [],
        expr.body)
end

function default_operating_point(affsys::AffectSystem)
    sys = system(affsys)

    op = AnyDict(unknowns(sys) .=> 0.0)
    for p in parameters(sys)
        T = symtype(p)
        if T <: Number
            op[p] = false
        elseif T <: Array{<:Real} && is_sized_array_symbolic(p)
            op[p] = zeros(size(p))
        end
    end
    return op
end

"""
Compile an affect defined by a set of equations. Systems with algebraic equations will solve implicit discrete problems to obtain their next state. Systems without will generate functions that perform explicit updates.
"""
function compile_equational_affect(
        aff::Union{AffectSystem, Vector{Equation}}, sys; reset_jumps = false,
        eval_expression = false, eval_module = @__MODULE__, op = nothing, kwargs...)
    if aff isa AbstractVector
        aff = make_affect(
            aff; iv = get_iv(sys), warn_no_algebraic = false)
    end
    if op === nothing
        op = default_operating_point(aff)
    end
    affsys = system(aff)
    ps_to_update = discretes(aff)
    dvs_to_update = setdiff(unknowns(aff), getfield.(observed(sys), :lhs))

    obseqs, eqs = unhack_observed(observed(affsys), equations(affsys))
    if isempty(equations(affsys))
        update_eqs = Symbolics.fast_substitute(
            obseqs, Dict([p => unPre(p) for p in parameters(affsys)]))
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
            [parameter_index(sys, p) for (i, p) in enumerate(lhss)
             if is_p[i]]
        else
            indexin((@view lhss[is_p]), ps)
        end
        _ps = reorder_parameters(sys, ps)
        integ = gensym(:MTKIntegrator)

        u_up,
        u_up! = build_function_wrapper(sys, (@view rhss[is_u]), dvs, _ps..., t;
            wrap_code = add_integrator_header(sys, integ, :u), expression = Val{false},
            outputidxs = u_idxs, wrap_mtkparameters, cse = false, eval_expression,
            eval_module)
        p_up,
        p_up! = build_function_wrapper(sys, (@view rhss[is_p]), dvs, _ps..., t;
            wrap_code = add_integrator_header(sys, integ, :p), expression = Val{false},
            outputidxs = p_idxs, wrap_mtkparameters, cse = false, eval_expression,
            eval_module)

        return let dvs_to_update = dvs_to_update, ps_to_update = ps_to_update,
            reset_jumps = reset_jumps, u_up! = u_up!, p_up! = p_up!

            function explicit_affect!(integ)
                isempty(dvs_to_update) || u_up!(integ)
                isempty(ps_to_update) || p_up!(integ)
                reset_jumps && reset_aggregated_jumps!(integ)
            end
        end
    else
        return let dvs_to_update = dvs_to_update, affsys = affsys,
            ps_to_update = ps_to_update, aff = aff, sys = sys, reset_jumps = reset_jumps

            dvs_to_access = unknowns(affsys)
            ps_to_access = [unPre(p) for p in parameters(affsys)]

            affu_getter = getsym(sys, dvs_to_access)
            affp_getter = getsym(sys, ps_to_access)
            affu_setter! = setsym(affsys, unknowns(affsys))
            affp_setter! = setsym(affsys, parameters(affsys))
            u_setter! = setsym(sys, dvs_to_update)
            p_setter! = setsym(sys, ps_to_update)
            u_getter = getsym(affsys, dvs_to_update)
            p_getter = getsym(affsys, ps_to_update)

            affprob = ImplicitDiscreteProblem(
                affsys, op,
                (0, 0);
                build_initializeprob = false, check_length = false, eval_expression,
                eval_module, check_compatibility = false, kwargs...)

            function implicit_affect!(integ)
                new_u0 = affu_getter(integ)
                affu_setter!(affprob, new_u0)
                new_ps = affp_getter(integ)
                affp_setter!(affprob, new_ps)

                affprob = remake(
                    affprob, tspan = (integ.t, integ.t))
                affsol = init(affprob, IDSolve())
                (check_error(affsol) === ReturnCode.InitialFailure) &&
                    throw(UnsolvableCallbackError(all_equations(aff)))

                u_setter!(integ, u_getter(affsol))
                p_setter!(integ, p_getter(affsol))

                reset_jumps && reset_aggregated_jumps!(integ)
            end
        end
    end
end

struct UnsolvableCallbackError
    eqs::Vector{Equation}
end

function Base.showerror(io::IO, err::UnsolvableCallbackError)
    println(io,
        "The callback defined by the following equations:\n\n$(join(err.eqs, "\n"))\n\nis not solvable. Please check that the algebraic equations and affect equations are correct, and that all parameters intended to be changed are passed in as `discrete_parameters`.")
end

merge_cb(::Nothing, ::Nothing) = nothing
merge_cb(::Nothing, x) = merge_cb(x, nothing)
merge_cb(x, ::Nothing) = x
merge_cb(x, y) = CallbackSet(x, y)

"""
Generate the CallbackSet for a ODESystem or SDESystem.
"""
function process_events(sys; callback = nothing, kwargs...)
    contin_cbs = generate_continuous_callbacks(sys; kwargs...)
    discrete_cbs = generate_discrete_callbacks(sys; kwargs...)
    cb = merge_cb(contin_cbs, callback)
    (discrete_cbs === nothing) ? cb : CallbackSet(contin_cbs, discrete_cbs...)
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
    obs = get_discrete_events(sys)
    systems = get_systems(sys)
    cbs = [obs;
           reduce(vcat,
               (map(cb -> namespace_callback(cb, s), discrete_events(s)) for s in systems),
               init = SymbolicDiscreteCallback[])]
    cbs
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
    getfield(sys, :discrete_events)
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
    obs = get_continuous_events(sys)
    filter(!isempty, obs)

    systems = get_systems(sys)
    cbs = [obs;
           reduce(vcat,
               (map(o -> namespace_callback(o, s), continuous_events(s)) for s in systems),
               init = SymbolicContinuousCallback[])]
    filter(!isempty, cbs)
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
    getfield(sys, :continuous_events)
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
    cont_callbacks, disc_callbacks
end
