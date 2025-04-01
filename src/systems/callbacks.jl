abstract type AbstractCallback end

struct FunctionalAffect
    f::Any
    sts::Vector
    sts_syms::Vector{Symbol}
    pars::Vector
    pars_syms::Vector{Symbol}
    discretes::Vector
    ctx::Any
end

function FunctionalAffect(f, sts, pars, discretes, ctx = nothing)
    # sts & pars contain either pairs: resistor.R => R, or Syms: R
    vs = [x isa Pair ? x.first : x for x in sts]
    vs_syms = Symbol[x isa Pair ? Symbol(x.second) : getname(x) for x in sts]
    length(vs_syms) == length(unique(vs_syms)) || error("Variables are not unique")

    ps = [x isa Pair ? x.first : x for x in pars]
    ps_syms = Symbol[x isa Pair ? Symbol(x.second) : getname(x) for x in pars]
    length(ps_syms) == length(unique(ps_syms)) || error("Parameters are not unique")

    FunctionalAffect(f, vs, vs_syms, ps, ps_syms, discretes, ctx)
end

function FunctionalAffect(; f, sts, pars, discretes, ctx = nothing)
    FunctionalAffect(f, sts, pars, discretes, ctx)
end

func(a::FunctionalAffect) = a.f
context(a::FunctionalAffect) = a.ctx
parameters(a::FunctionalAffect) = a.pars
parameters_syms(a::FunctionalAffect) = a.pars_syms
unknowns(a::FunctionalAffect) = a.sts
unknowns_syms(a::FunctionalAffect) = a.sts_syms
discretes(a::FunctionalAffect) = a.discretes

function Base.:(==)(a1::FunctionalAffect, a2::FunctionalAffect)
    isequal(a1.f, a2.f) && isequal(a1.sts, a2.sts) && isequal(a1.pars, a2.pars) &&
        isequal(a1.sts_syms, a2.sts_syms) && isequal(a1.pars_syms, a2.pars_syms) &&
        isequal(a1.ctx, a2.ctx)
end

function Base.hash(a::FunctionalAffect, s::UInt)
    s = hash(a.f, s)
    s = hash(a.sts, s)
    s = hash(a.sts_syms, s)
    s = hash(a.pars, s)
    s = hash(a.pars_syms, s)
    s = hash(a.discretes, s)
    hash(a.ctx, s)
end

function has_functional_affect(cb)
    (affects(cb) isa FunctionalAffect || affects(cb) isa ImperativeAffect)
end

struct AffectSystem
    system::ImplicitDiscreteSystem
    unknowns::Vector
    parameters::Vector
    discretes::Vector
    """Maps the symbols of unknowns/observed in the ImplicitDiscreteSystem to its corresponding unknown/parameter in the parent system."""
    aff_to_sys::Dict
end

system(a::AffectSystem) = a.system
discretes(a::AffectSystem) = a.discretes
unknowns(a::AffectSystem) = a.unknowns
parameters(a::AffectSystem) = a.parameters
aff_to_sys(a::AffectSystem) = a.aff_to_sys
previous_vals(a::AffectSystem) = parameters(system(a))
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
        isequal(parameters(a1), parameters(a2)) &&
        isequal(aff_to_sys(a1), aff_to_sys(a2))
end

function Base.hash(a::AffectSystem, s::UInt)
    s = hash(system(a), s)
    s = hash(unknowns(a), s)
    s = hash(parameters(a), s)
    s = hash(discretes(a), s)
    hash(aff_to_sys(a), s)
end

function vars!(vars, aff::Union{FunctionalAffect, AffectSystem}; op = Differential)
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
SymbolicUtils.promote_symtype(::Type{Pre}, T) = T
SymbolicUtils.isbinop(::Pre) = false
Base.nameof(::Pre) = :Pre
Base.show(io::IO, x::Pre) = print(io, "Pre")
input_timedomain(::Pre, _ = nothing) = ContinuousClock()
output_timedomain(::Pre, _ = nothing) = ContinuousClock()
unPre(x::Num) = unPre(unwrap(x))
unPre(x::BasicSymbolic) = (iscall(x) && operation(x) isa Pre) ? only(arguments(x)) : x

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

###############################
###### Continuous events ######
###############################
const Affect = Union{AffectSystem, FunctionalAffect, ImperativeAffect}

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
    affect::Union{Affect, Nothing}
    affect_neg::Union{Affect, Nothing}
    initialize::Union{Affect, Nothing}
    finalize::Union{Affect, Nothing}
    rootfind::Union{Nothing, SciMLBase.RootfindOpt}
    reinitializealg::SciMLBase.DAEInitializationAlgorithm

    function SymbolicContinuousCallback(
            conditions::Union{Equation, Vector{Equation}},
            affect = nothing;
            discrete_parameters = Any[],
            affect_neg = affect,
            initialize = nothing,
            finalize = nothing,
            rootfind = SciMLBase.LeftRootFind,
            reinitializealg = nothing,
            iv = nothing,
            alg_eqs = Equation[])
        conditions = (conditions isa AbstractVector) ? conditions : [conditions]

        if isnothing(reinitializealg)
            any(a -> (a isa FunctionalAffect || a isa ImperativeAffect),
                [affect, affect_neg, initialize, finalize]) ?
            reinitializealg = SciMLBase.CheckInit() :
            reinitializealg = SciMLBase.NoInit()
        end

        new(conditions, make_affect(affect; iv, alg_eqs, discrete_parameters),
            make_affect(affect_neg; iv, alg_eqs, discrete_parameters),
            make_affect(initialize; iv, alg_eqs, discrete_parameters), make_affect(
                finalize; iv, alg_eqs, discrete_parameters),
            rootfind, reinitializealg)
    end # Default affect to nothing
end

function SymbolicContinuousCallback(p::Pair, args...; kwargs...)
    SymbolicContinuousCallback(p[1], p[2], args...; kwargs...)
end

function SymbolicContinuousCallback(cb::SymbolicContinuousCallback, args...;
        iv = nothing, alg_eqs = Equation[], kwargs...)
    cb
end

make_affect(affect::Nothing; kwargs...) = nothing
make_affect(affect::Tuple; kwargs...) = FunctionalAffect(affect...)
make_affect(affect::NamedTuple; kwargs...) = FunctionalAffect(; affect...)
make_affect(affect::Affect; kwargs...) = affect

function make_affect(affect::Vector{Equation}; discrete_parameters = Any[],
        iv = nothing, alg_eqs::Vector{Equation} = Equation[])
    isempty(affect) && return nothing
    isempty(alg_eqs) &&
        @warn "No algebraic equations were found for the callback defined by $(join(affect, ", ")). If the system has no algebraic equations, this can be disregarded. Otherwise pass in `alg_eqs` to the SymbolicContinuousCallback constructor."
    if isnothing(iv)
        iv = t_nounits
        @warn "No independent variable specified. Defaulting to t_nounits."
    end

    for p in discrete_parameters
        occursin(unwrap(iv), unwrap(p)) ||
            error("Non-time dependent parameter $p passed in as a discrete. Must be declared as @parameters $p(t).")
    end

    dvs = OrderedSet()
    params = OrderedSet()
    for eq in affect
        if !haspre(eq) && !(symbolic_type(eq.rhs) === NotSymbolic() ||
             symbolic_type(eq.lhs) === NotSymbolic())
            @warn "Affect equation $eq has no `Pre` operator. As such it will be interpreted as an algebraic equation to be satisfied after the callback. If you intended to use the value of a variable x before the affect, use Pre(x)."
        end
        collect_vars!(dvs, params, eq, iv; op = Pre)
        diffvs = collect_applied_operators(eq, Differential)
        union!(dvs, diffvs)
    end
    for eq in alg_eqs
        collect_vars!(dvs, params, eq, iv)
    end

    pre_params = filter(haspre ∘ value, params)
    sys_params = collect(setdiff(params, union(discrete_parameters, pre_params)))
    discretes = map(tovar, discrete_parameters)
    dvs = collect(dvs)
    _dvs = map(default_toterm, dvs)

    aff_map = Dict(zip(discretes, discrete_parameters))
    rev_map = Dict(zip(discrete_parameters, discretes))
    subs = merge(rev_map, Dict(zip(dvs, _dvs)))
    affect = Symbolics.fast_substitute(affect, subs)
    alg_eqs = Symbolics.fast_substitute(alg_eqs, subs)

    @named affectsys = ImplicitDiscreteSystem(
        vcat(affect, alg_eqs), iv, collect(union(_dvs, discretes)),
        collect(union(pre_params, sys_params)))
    affectsys = structural_simplify(affectsys; fully_determined = false)
    # get accessed parameters p from Pre(p) in the callback parameters
    accessed_params = filter(isparameter, map(unPre, collect(pre_params)))
    union!(accessed_params, sys_params)

    # add scalarized unknowns to the map.
    _dvs = reduce(vcat, map(scalarize, _dvs), init = Any[])
    @show _dvs
    for u in _dvs
        aff_map[u] = u
    end

    AffectSystem(affectsys, collect(_dvs), collect(accessed_params),
        collect(discrete_parameters), aff_map)
end

function make_affect(affect; kwargs...)
    error("Malformed affect $(affect). This should be a vector of equations or a tuple specifying a functional affect.")
end

"""
Generate continuous callbacks.
"""
function SymbolicContinuousCallbacks(events; discrete_parameters = Any[],
        alg_eqs::Vector{Equation} = Equation[], iv = nothing)
    callbacks = SymbolicContinuousCallback[]
    isnothing(events) && return callbacks

    events isa AbstractVector || (events = [events])
    isempty(events) && return callbacks

    for event in events
        cond, affs = event isa Pair ? (event[1], event[2]) : (event, nothing)
        push!(callbacks,
            SymbolicContinuousCallback(cond, affs; iv, alg_eqs, discrete_parameters))
    end
    callbacks
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

# TODO: Iterative callbacks
"""
    SymbolicDiscreteCallback(conditions::Vector{Equation}, affect = nothing, iv = nothing;
                             initialize = nothing, finalize = nothing, alg_eqs = Equation[])

A callback that triggers at the first timestep that the conditions are satisfied.

The condition can be one of: 
- Δt::Real              - periodic events with period Δt
- ts::Vector{Real}      - events trigger at these preset times given by `ts`
- eqs::Vector{Equation} - events trigger when the condition evaluates to true

Arguments: 
- iv: The independent variable of the system. This must be specified if the independent variable appaers in one of the equations explicitly, as in x ~ t + 1.
- alg_eqs: Algebraic equations of the system that must be satisfied after the callback occurs.
"""
struct SymbolicDiscreteCallback <: AbstractCallback
    conditions::Any
    affect::Union{Affect, Nothing}
    initialize::Union{Affect, Nothing}
    finalize::Union{Affect, Nothing}
    reinitializealg::SciMLBase.DAEInitializationAlgorithm

    function SymbolicDiscreteCallback(
            condition, affect = nothing;
            initialize = nothing, finalize = nothing, iv = nothing,
            alg_eqs = Equation[], discrete_parameters = Any[], reinitializealg = nothing)
        c = is_timed_condition(condition) ? condition : value(scalarize(condition))

        if isnothing(reinitializealg)
                reinitializealg = SciMLBase.CheckInit() : 
                reinitializealg = SciMLBase.NoInit()
        end
        new(c, make_affect(affect; iv, alg_eqs, discrete_parameters),
            make_affect(initialize; iv, alg_eqs, discrete_parameters),
            make_affect(finalize; iv, alg_eqs, discrete_parameters), reinitializealg)
    end # Default affect to nothing
end

function SymbolicDiscreteCallback(p::Pair, args...; kwargs...)
    SymbolicDiscreteCallback(p[1], p[2], args...; kwargs...)
end
SymbolicDiscreteCallback(cb::SymbolicDiscreteCallback, args...; kwargs...) = cb

"""
Generate discrete callbacks.
"""
function SymbolicDiscreteCallbacks(events; discrete_parameters::Vector = Any[],
        alg_eqs::Vector{Equation} = Equation[], iv = nothing)
    callbacks = SymbolicDiscreteCallback[]

    isnothing(events) && return callbacks
    events isa AbstractVector || (events = [events])
    isempty(events) && return callbacks

    for event in events
        cond, affs = event isa Pair ? (event[1], event[2]) : (event, nothing)
        push!(callbacks,
            SymbolicDiscreteCallback(cond, affs; iv, alg_eqs, discrete_parameters))
    end
    callbacks
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

############################################
########## Namespacing Utilities ###########
############################################

function namespace_affects(affect::FunctionalAffect, s)
    FunctionalAffect(func(affect),
        renamespace.((s,), unknowns(affect)),
        unknowns_syms(affect),
        renamespace.((s,), parameters(affect)),
        parameters_syms(affect),
        renamespace.((s,), discretes(affect)),
        context(affect))
end

function namespace_affects(affect::AffectSystem, s)
    AffectSystem(renamespace(s, system(affect)),
        renamespace.((s,), unknowns(affect)),
        renamespace.((s,), parameters(affect)),
        renamespace.((s,), discretes(affect)),
        Dict([k => renamespace(s, v) for (k, v) in aff_to_sys(affect)]))
end
namespace_affects(af::Nothing, s) = nothing

function namespace_callback(cb::SymbolicContinuousCallback, s)::SymbolicContinuousCallback
    SymbolicContinuousCallback(
        namespace_equation.(equations(cb), (s,)),
        namespace_affects(affects(cb), s),
        affect_neg = namespace_affects(affect_negs(cb), s),
        initialize = namespace_affects(initialize_affects(cb), s),
        finalize = namespace_affects(finalize_affects(cb), s),
        rootfind = cb.rootfind)
end

function namespace_conditions(condition, s)
    is_timed_condition(condition) ? condition : namespace_expr(condition, s)
end

function namespace_callback(cb::SymbolicDiscreteCallback, s)::SymbolicDiscreteCallback
    SymbolicDiscreteCallback(
        namespace_conditions(conditions(cb), s),
        namespace_affects(affects(cb), s),
        initialize = namespace_affects(initialize_affects(cb), s),
        finalize = namespace_affects(finalize_affects(cb), s))
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
    (is_discrete(e1) === is_discrete(e2)) || return false
    (isequal(e1.conditions, e2.conditions) && isequal(e1.affect, e2.affect) &&
     isequal(e1.initialize, e2.initialize) && isequal(e1.finalize, e2.finalize)) &&
        isequal(e1.reinitializealg, e2.reinitializealg) ||
        return false
    is_discrete(e1) ||
        (isequal(e1.affect_neg, e2.affect_neg) && isequal(e1.rootfind, e2.rootfind))
end

Base.isempty(cb::AbstractCallback) = isempty(cb.conditions)

####################################
####### Compilation functions ######
####################################
"""
    compile_condition(cb::AbstractCallback, sys, dvs, ps; expression, kwargs...)

Returns a function `condition(u,t,integrator)`, condition(out,u,t,integrator)` returning the `condition(cb)`.

Notes

  - `expression = Val{true}`, causes the generated function to be returned as an expression.
    If  set to `Val{false}` a `RuntimeGeneratedFunction` will be returned.
  - `kwargs` are passed through to `Symbolics.build_function`.
"""
function compile_condition(
        cbs::Union{AbstractCallback, Vector{<:AbstractCallback}}, sys, dvs, ps;
        expression = Val{false}, eval_expression = false, eval_module = @__MODULE__, kwargs...)
    u = map(x -> time_varying_as_func(value(x), sys), dvs)
    p = map.(x -> time_varying_as_func(value(x), sys), reorder_parameters(sys, ps))
    t = get_iv(sys)
    condit = conditions(cbs)
    cs = collect_constants(condit)
    if !isempty(cs)
        cmap = map(x -> x => getdefault(x), cs)
        condit = substitute(condit, Dict(cmap))
    end

    if !is_discrete(cbs)
        condit = reduce(vcat, flatten_equations(condit))
        condit = condit isa AbstractVector ? [c.lhs - c.rhs for c in condit] :
                 [condit.lhs - condit.rhs]
    end

    fs = build_function_wrapper(sys,
        condit, u, p..., t; expression,
        kwargs...)

    if expression == Val{false}
        fs = eval_or_rgf.(fs; eval_expression, eval_module)
    end
    f_oop, f_iip = is_discrete(cbs) ? (fs, nothing) : fs # no iip function for discrete condition.

    cond = if cbs isa AbstractVector
        (out, u, t, integ) -> f_iip(out, u, parameter_values(integ), t)
    elseif is_discrete(cbs)
        (u, t, integ) -> f_oop(u, parameter_values(integ), t)
    else
        function (u, t, integ)
            if DiffEqBase.isinplace(integ.sol.prob)
                tmp, = DiffEqBase.get_tmp_cache(integ)
                f_iip(tmp, u, parameter_values(integ), t)
                tmp[1]
            else
                f_oop(u, parameter_values(integ), t)
            end
        end
    end
end

"""
Compile user-defined functional affect.
"""
function compile_functional_affect(affect::FunctionalAffect, sys; kwargs...)
    dvs = unknowns(sys)
    ps = parameters(sys)
    dvs_ind = Dict(reverse(en) for en in enumerate(dvs))
    v_inds = map(sym -> dvs_ind[sym], unknowns(affect))

    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        p_inds = [(pind = parameter_index(sys, sym)) === nothing ? sym : pind
                  for sym in parameters(affect)]
    else
        ps_ind = Dict(reverse(en) for en in enumerate(ps))
        p_inds = map(sym -> get(ps_ind, sym, sym), parameters(affect))
    end
    # HACK: filter out eliminated symbols. Not clear this is the right thing to do
    # (MTK should keep these symbols)
    u = filter(x -> !isnothing(x[2]), collect(zip(unknowns_syms(affect), v_inds))) |>
        NamedTuple
    p = filter(x -> !isnothing(x[2]), collect(zip(parameters_syms(affect), p_inds))) |>
        NamedTuple

    let u = u, p = p, user_affect = func(affect), ctx = context(affect)
        (integ) -> begin
            user_affect(integ, u, p, ctx)
        end
    end
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

    trigger = compile_condition(
        cbs, sys, unknowns(sys), parameters(sys; initial_parameters = true); kwargs...)
    affects = []
    affect_negs = []
    inits = []
    finals = []
    for cb in cbs
        affect = compile_affect(cb.affect, cb, sys; default = EMPTY_AFFECT, kwargs...)
        push!(affects, affect)
        affect_neg = (cb.affect_neg === cb.affect) ? affect :
                     compile_affect(
            cb.affect_neg, cb, sys; default = EMPTY_AFFECT, kwargs...)
        push!(affect_negs, affect_neg)
        push!(inits,
            compile_affect(
                cb.initialize, cb, sys; default = nothing, is_init = true, kwargs...))
        push!(finals, compile_affect(cb.finalize, cb, sys; default = nothing, kwargs...))
    end

    # Since there may be different number of conditions and affects,
    # we build a map that translates the condition eq. number to the affect number
    eq2affect = reduce(vcat,
        [fill(i, num_eqs[i]) for i in eachindex(affects)])
    eqs = reduce(vcat, eqs)

    affect = function (integ, idx)
        affects[eq2affect[idx]](integ)
    end
    affect_neg = function (integ, idx)
        f = affect_negs[eq2affect[idx]]
        isnothing(f) && return
        f(integ)
    end
    initialize = wrap_vector_optional_affect(inits, SciMLBase.INITIALIZE_DEFAULT)
    finalize = wrap_vector_optional_affect(finals, SciMLBase.FINALIZE_DEFAULT)

    return VectorContinuousCallback(
        trigger, affect, affect_neg, length(eqs); initialize, finalize,
        rootfind = cbs[1].rootfind, initializealg = cbs[1].reinitializealg)
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
        is_init = true, kwargs...)
    final = compile_affect(
        cb.finalize, cb, sys; default = SciMLBase.FINALIZE_DEFAULT, kwargs...)

    initialize = isnothing(cb.initialize) ? init : ((c, u, t, i) -> init(i))
    finalize = isnothing(cb.finalize) ? final : ((c, u, t, i) -> final(i))

    if is_discrete(cb)
        if is_timed && conditions(cb) isa AbstractVector
            return PresetTimeCallback(trigger, affect; initialize,
                finalize, initializealg = cb.reinitializealg)
        elseif is_timed
            return PeriodicCallback(
                affect, trigger; initialize, finalize, initializealg = cb.reinitializealg)
        else
            return DiscreteCallback(trigger, affect; initialize,
                finalize, initializealg = cb.reinitializealg)
        end
    else
        return ContinuousCallback(trigger, affect, affect_neg; initialize, finalize,
            rootfind = cb.rootfind, initializealg = cb.reinitializealg)
    end
end

"""
    compile_affect(cb::AbstractCallback, sys::AbstractSystem, dvs, ps; expression, outputidxs, kwargs...)

Returns a function that takes an integrator as argument and modifies the state with the
affect. The generated function has the signature `affect!(integrator)`.

Notes

  - `expression = Val{true}`, causes the generated function to be returned as an expression.
    If set to `Val{false}` a `RuntimeGeneratedFunction` will be returned.
  - `outputidxs`, a vector of indices of the output variables which should correspond to
    `unknowns(sys)`. If provided, checks that the LHS of affect equations are variables are
    dropped, i.e. it is assumed these indices are correct and affect equations are
    well-formed.
  - `kwargs` are passed through to `Symbolics.build_function`.
"""
function compile_affect(
        aff::Union{Nothing, Affect}, cb::AbstractCallback, sys::AbstractSystem;
        default = nothing, is_init = false, kwargs...)
    save_idxs = if !(has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing)
        Int[]
    else
        get(ic.callback_to_clocks, cb, Int[])
    end

    if isnothing(aff)
        is_init ? wrap_save_discretes(default, save_idxs) : default
    elseif aff isa AffectSystem
        f = compile_equational_affect(aff, sys; kwargs...)
        wrap_save_discretes(f, save_idxs)
    elseif aff isa FunctionalAffect || aff isa ImperativeAffect
        f = compile_functional_affect(aff, sys; kwargs...)
        wrap_save_discretes(f, save_idxs)
    end
end

function wrap_save_discretes(f, save_idxs)
    let save_idxs = save_idxs
        if f === SciMLBase.INITIALIZE_DEFAULT
            (c, u, t, i) -> begin
                isnothing(f) || f(c, u, t, i)
                for idx in save_idxs
                    SciMLBase.save_discretes!(i, idx)
                end
            end
        else
            (i) -> begin
                isnothing(f) || f(i)
                for idx in save_idxs
                    SciMLBase.save_discretes!(i, idx)
                end
            end
        end
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

"""
Compile an affect defined by a set of equations. Systems with algebraic equations will solve implicit discrete problems to obtain their next state. Systems without will generate functions that perform explicit updates.
"""
function compile_equational_affect(
        aff::Union{AffectSystem, Vector{Equation}}, sys; reset_jumps = false, kwargs...)
    if aff isa AbstractVector
        aff = make_affect(aff; iv = get_iv(sys))
    end
    affsys = system(aff)
    ps_to_update = discretes(aff)
    dvs_to_update = setdiff(unknowns(aff), getfield.(observed(sys), :lhs))
    aff_map = aff_to_sys(aff)
    sys_map = Dict([v => k for (k, v) in aff_map])

    if isempty(equations(affsys))
        update_eqs = Symbolics.fast_substitute(
            observed(affsys), Dict([p => unPre(p) for p in parameters(affsys)]))
        rhss = map(x -> x.rhs, update_eqs)
        lhss = map(x -> aff_map[x.lhs], update_eqs)
        is_p = [lhs ∈ Set(ps_to_update) for lhs in lhss]
        is_u = [lhs ∈ Set(dvs_to_update) for lhs in lhss]
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

        u_up, u_up! = build_function_wrapper(sys, (@view rhss[is_u]), dvs, _ps..., t;
            wrap_code = add_integrator_header(sys, integ, :u),
            expression = Val{false}, outputidxs = u_idxs, wrap_mtkparameters)
        p_up, p_up! = build_function_wrapper(sys, (@view rhss[is_p]), dvs, _ps..., t;
            wrap_code = add_integrator_header(sys, integ, :p),
            expression = Val{false}, outputidxs = p_idxs, wrap_mtkparameters)

        return function explicit_affect!(integ)
            isempty(dvs_to_update) || u_up!(integ)
            isempty(ps_to_update) || p_up!(integ)
            reset_jumps && reset_aggregated_jumps!(integ)
        end
    else
        return let dvs_to_update = dvs_to_update, aff_map = aff_map, sys_map = sys_map,
            affsys = affsys, ps_to_update = ps_to_update, aff = aff

            function implicit_affect!(integ)
                pmap = Pair[]
                for pre_p in parameters(affsys)
                    p = unPre(pre_p)
                    pval = isparameter(p) ? integ.ps[p] : integ[p]
                    push!(pmap, pre_p => pval)
                end
                u0 = Pair[]
                for u in unknowns(affsys)
                    uval = isparameter(aff_map[u]) ? integ.ps[aff_map[u]] : integ[u]
                    push!(u0, u => uval)
                end
                affprob = ImplicitDiscreteProblem(affsys, u0, (integ.t, integ.t), pmap;
                    build_initializeprob = false, check_length = false)
                affsol = init(affprob, IDSolve())
                (check_error(affsol) === ReturnCode.InitialFailure) &&
                    throw(UnsolvableCallbackError(all_equations(aff)))
                for u in dvs_to_update
                    integ[u] = affsol[sys_map[u]]
                end
                for p in ps_to_update
                    integ.ps[p] = affsol[sys_map[p]]
                end
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

has_discrete_events(sys::AbstractSystem) = isdefined(sys, :discrete_events)
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

has_continuous_events(sys::AbstractSystem) = isdefined(sys, :continuous_events)
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
