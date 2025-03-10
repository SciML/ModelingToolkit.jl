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

function vars!(vars, aff::FunctionalAffect; op = Differential)
    for var in Iterators.flatten((unknowns(aff), parameters(aff), discretes(aff)))
        vars!(vars, var)
    end
    return vars
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
        Symbolics.array_term(p, toparam(x))
    elseif iscall(x) && operation(x) == getindex
        # instead of `Pre(x[1])` create `Pre(x)[1]`
        # which allows parameter indexing to handle this case automatically.
        arr = arguments(x)[1]
        term(getindex, p(toparam(arr)), arguments(x)[2:end]...)
    else
        term(p, toparam(x))
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
const Affect = Union{ImplicitDiscreteSystem, FunctionalAffect, ImperativeAffect}

"""
    SymbolicContinuousCallback(eqs::Vector{Equation}, affect, affect_neg, rootfind)

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

DAEs will automatically be reinitialized.

Initial and final affects can also be specified with SCC, which are specified identically to positive and negative edge affects. Initialization affects
will run as soon as the solver starts, while finalization affects will be executed after termination.
"""
struct SymbolicContinuousCallback <: AbstractCallback
    conditions::Vector{Equation}
    affect::Union{Affect, Nothing}
    affect_neg::Union{Affect, Nothing}
    initialize::Union{Affect, Nothing}
    finalize::Union{Affect, Nothing}
    rootfind::Union{Nothing, SciMLBase.RootfindOpt}

    function SymbolicContinuousCallback(
            conditions::Vector{Equation},
            affect = nothing;
            affect_neg = affect,
            initialize = nothing,
            finalize = nothing,
            rootfind = SciMLBase.LeftRootFind)
        new(eqs, initialize, finalize, make_affect(affect),
            make_affect(affect_neg), rootfind)
    end # Default affect to nothing
end

make_affect(affect::Tuple, iv) = FunctionalAffect(affects...)
make_affect(affect::NamedTuple, iv) = FunctionalAffect(; affects...)
make_affect(affect::FunctionalAffect, iv) = affect

function make_affect(affect::Vector{Equation}, iv; warn = true)
    affect = scalarize(affect)
    unknowns = OrderedSet()
    params = OrderedSet()
    for eq in affect
        !haspre(eq) && warn && @warn "Equation $eq has no `Pre` operator. As such it will be interpreted as an algebraic equation to be satisfied after the callback. 
            If you intended to use the value of a variable x before the affect, use Pre(x)."
        collect_vars!(unknowns, params, eq, iv; op = Pre)
    end

    # System parameters should become unknowns.
    cb_params = OrderedSet()
    sys_params = OrderedSet()
    for p in params
        if iscall(p) && (operator(p) isa Pre)
            push!(cb_params, p)
        else
            p = Sym{FnType{Tuple{symtype(iv)}, Real}}(nameof(p))(iv)
            p = wrap(tovar(p))
            push!(sys_params, p)
        end
    end
    
    @mtkbuild affect = ImplicitDiscreteSystem(affect, iv, vcat(unknowns, sys_params), cb_params)
end

make_affect(affect, iv) = error("Malformed affect $(affect). This should be a vector of equations or a tuple specifying a functional affect.")

"""
Generate continuous callbacks.
""" 
function SymbolicContinuousCallbacks(events, algeeqs, iv)
    callbacks = SymbolicContinuousCallback[]
    (isnothing(events) || isempty(events)) && return callbacks

    events isa AbstractVector || (events = [events])
    for (cond, affs) in events
        if affs isa AbstractVector
            affs = vcat(affs, algeeqs)
        end
        affect = make_affect(affs, iv)
        push!(callbacks, SymbolicContinuousCallback(cond, affect, affect, nothing, nothing, SciMLBase.LeftRootFind)) 
    end
    callbacks
end

function Base.show(io::IO, cb::SymbolicContinuousCallback)
    indent = get(io, :indent, 0)
    iio = IOContext(io, :indent => indent + 1)
    print(io, "SymbolicContinuousCallback(")
    print(iio, "Equations:")
    show(iio, equations(cb))
    print(iio, "; ")
    if affects(cb) != NULL_AFFECT
        print(iio, "Affect:")
        show(iio, affects(cb))
        print(iio, ", ")
    end
    if affect_negs(cb) != NULL_AFFECT
        print(iio, "Negative-edge affect:")
        show(iio, affect_negs(cb))
        print(iio, ", ")
    end
    if initialize_affects(cb) != NULL_AFFECT
        print(iio, "Initialization affect:")
        show(iio, initialize_affects(cb))
        print(iio, ", ")
    end
    if finalize_affects(cb) != NULL_AFFECT
        print(iio, "Finalization affect:")
        show(iio, finalize_affects(cb))
    end
    print(iio, ")")
end

function Base.show(io::IO, mime::MIME"text/plain", cb::SymbolicContinuousCallback)
    indent = get(io, :indent, 0)
    iio = IOContext(io, :indent => indent + 1)
    println(io, "SymbolicContinuousCallback:")
    println(iio, "Equations:")
    show(iio, mime, equations(cb))
    print(iio, "\n")
    if affects(cb) != NULL_AFFECT
        println(iio, "Affect:")
        show(iio, mime, affects(cb))
        print(iio, "\n")
    end
    if affect_negs(cb) != NULL_AFFECT
        println(iio, "Negative-edge affect:")
        show(iio, mime, affect_negs(cb))
        print(iio, "\n")
    end
    if initialize_affects(cb) != NULL_AFFECT
        println(iio, "Initialization affect:")
        show(iio, mime, initialize_affects(cb))
        print(iio, "\n")
    end
    if finalize_affects(cb) != NULL_AFFECT
        println(iio, "Finalization affect:")
        show(iio, mime, finalize_affects(cb))
        print(iio, "\n")
    end
end

function vars!(vars, cb::SymbolicContinuousCallback; op = Differential)
    for eq in equations(cb)
        vars!(vars, eq; op)
    end
    for aff in (affects(cb), affect_negs(cb), initialize_affects(cb), finalize_affects(cb))
        if aff isa Vector{Equation}
            for eq in aff
                vars!(vars, eq; op)
            end
        elseif aff !== nothing
            vars!(vars, aff; op)
        end
    end
    return vars
end

################################
######## Discrete events #######
################################

# TODO: Iterative callbacks
"""
    SymbolicDiscreteCallback(conditions::Vector{Equation}, affect)

A callback that triggers at the first timestep that the conditions are satisfied.

The condition can be one of: 
- Δt::Real              - periodic events with period Δt
- ts::Vector{Real}      - events trigger at these preset times given by `ts`
- eqs::Vector{Equation} - events trigger when the condition evaluates to true
"""
struct SymbolicDiscreteCallback <: AbstractCallback 
    conditions::Any
    affect::Affect
    initialize::Union{Affect, Nothing}
    finalize::Union{Affect, Nothing}

    function SymbolicDiscreteCallback(
            condition, affect = nothing;
            initialize = nothing, finalize = nothing)
        c = is_timed_condition(condition) ? condition : value(scalarize(condition))
        new(c, make_affect(affect), make_affect(initialize),
            make_affect(finalize))
    end # Default affect to nothing
end

"""
Generate discrete callbacks.
"""
function SymbolicDiscreteCallbacks(events, algeeqs, iv)
    callbacks = SymbolicDiscreteCallback[]
    (isnothing(events) || isempty(events)) && return callbacks
    events isa AbstractVector || (events = [events])

    for (cond, aff) in events
        if aff isa AbstractVector
            aff = vcat(aff, algeeqs)
        end
        affect = make_affect(aff, iv)
        push!(callbacks, SymbolicDiscreteCallback(cond, affect, nothing, nothing))
    end
    callbacks
end

function is_timed_condition(condition::T) where T
    if T <: Real
        true
    elseif T <: AbstractVector
        eltype(condition) <: Real
    else
        false
    end
end

function Base.show(io::IO, db::SymbolicDiscreteCallback)
    indent = get(io, :indent, 0)
    iio = IOContext(io, :indent => indent + 1)
    println(io, "SymbolicDiscreteCallback:")
    println(iio, "Conditions:")
    print(iio, "; ")
    if affects(db) != nothing 
        print(iio, "Affect:")
        show(iio, affects(db))
        print(iio, ", ")
    end
    if initialize_affects(db) != nothing 
        print(iio, "Initialization affect:")
        show(iio, initialize_affects(db))
        print(iio, ", ")
    end
    if finalize_affects(db) != nothing 
        print(iio, "Finalization affect:")
        show(iio, finalize_affects(db))
    end
    print(iio, ")")
end

function vars!(vars, cb::SymbolicDiscreteCallback; op = Differential)
    if symbolic_type(cb.condition) == NotSymbolic
        if cb.condition isa AbstractArray
            for eq in cb.condition
                vars!(vars, eq; op)
            end
        end
    else
        vars!(vars, cb.condition; op)
    end
    for aff in (cb.affects, cb.initialize, cb.finalize)
        if aff isa Vector{Equation}
            for eq in aff
                vars!(vars, eq; op)
            end
        elseif aff !== nothing
            vars!(vars, aff; op)
        end
    end
    return vars
end

############################################
########## Namespacing Utilities ###########
############################################

function namespace_affect(affect::FunctionalAffect, s)
    FunctionalAffect(func(affect),
        renamespace.((s,), unknowns(affect)),
        unknowns_syms(affect),
        renamespace.((s,), parameters(affect)),
        parameters_syms(affect),
        renamespace.((s,), discretes(affect)),
        context(affect))
end

function namespace_affects(af::Affect, s)
    if af isa ImplicitDiscreteSystem
        af
    elseif af isa FunctionalAffect || af isa ImperativeAffect
        namespace_affect(af, s)
    else
        nothing
    end
end

function namespace_callback(cb::SymbolicContinuousCallback, s)::SymbolicContinuousCallback
    SymbolicContinuousCallback(
        namespace_equation.(equations(cb), (s,)),
        namespace_affects(affects(cb), s),
        namespace_affects(affect_negs(cb), s),
        namespace_affects(initialize_affects(cb), s),
        namespace_affects(finalize_affects(cb), s),
        cb.rootfind)
end

function namespace_condition(condition, s)
    is_timed_condition(condition) ? condition : namespace_expr(condition, s)
end

function namespace_callback(cb::SymbolicDiscreteCallback, s)::SymbolicDiscreteCallback
    SymbolicDiscreteCallback(
        namespace_condition(condition(cb), s), 
        namespace_affects(affects(cb), s),
        namespace_affects(initialize_affects(cb), s),
        namespace_affects(finalize_affects(cb), s))
end

function Base.hash(cb::SymbolicContinuousCallback, s::UInt)
    s = foldr(hash, cb.eqs, init = s)
    s = hash(cb.affect, s)
    s = hash(cb.affect_neg, s)
    s = hash(cb.initialize, s)
    s = hash(cb.finalize, s)
    hash(cb.rootfind, s)
end

function Base.hash(cb::SymbolicDiscreteCallback, s::UInt)
    s = hash(cb.condition, s)
    s = hash(cb.affects, s)
    s = hash(cb.initialize, s)
    hash(cb.finalize, s)
end

###########################
######### Helpers #########
###########################

conditions(cb::AbstractCallback) = cb.conditions
conditions(cbs::Vector{<:AbstractCallback}) = reduce(vcat, conditions(cb) for cb in cbs; init = [])
equations(cb::AbstractCallback) = conditions(cb)
equations(cb::Vector{<:AbstractCallback}) = conditions(cb)

affects(cb::AbstractCallback) = cb.affect
affects(cbs::Vector{<:AbstractCallback}) = reduce(vcat, affects(cb) for cb in cbs; init = [])

affect_negs(cb::SymbolicContinuousCallback) = cb.affect_neg
affect_negs(cbs::Vector{SymbolicContinuousCallback}) = reduce(vcat, affect_negs(cb) for cb in cbs; init = [])

initialize_affects(cb::AbstractCallback) = cb.initialize
initialize_affects(cbs::Vector{<:AbstractCallback}) = reduce(initialize_affects, vcat, cbs; init = [])

finalize_affects(cb::AbstractCallback) = cb.finalize
finalize_affects(cbs::Vector{<:AbstractCallback}) = reduce(finalize_affects, vcat, cbs; init = [])

function Base.:(==)(e1::SymbolicDiscreteCallback, e2::SymbolicDiscreteCallback)
    isequal(e1.condition, e2.condition) && isequal(e1.affects, e2.affects) &&
        isequal(e1.initialize, e2.initialize) && isequal(e1.finalize, e2.finalize)
end

function Base.:(==)(e1::SymbolicContinuousCallback, e2::SymbolicContinuousCallback)
    isequal(e1.eqs, e2.eqs) && isequal(e1.affect, e2.affect) &&
        isequal(e1.initialize, e2.initialize) && isequal(e1.finalize, e2.finalize) &&
        isequal(e1.affect_neg, e2.affect_neg) && isequal(e1.rootfind, e2.rootfind)
end

Base.isempty(cb::AbstractCallback) = isempty(cb.conditions)

####################################
####### Compilation functions ######
####################################
function condition_header(sys::AbstractSystem, integrator = gensym(:MTKIntegrator))
    expr -> Func(
        [expr.args[1], expr.args[2],
            DestructuredArgs(expr.args[3:end], integrator, inds = [:p])],
        [],
        expr.body)
end

"""
    compile_condition(cb::AbstractCallback, sys, dvs, ps; expression, kwargs...)

Returns a function `condition(u,t,integrator)` returning the `condition(cb)`.

Notes

  - `expression = Val{true}`, causes the generated function to be returned as an expression.
    If  set to `Val{false}` a `RuntimeGeneratedFunction` will be returned.
  - `kwargs` are passed through to `Symbolics.build_function`.
"""
function compile_condition(cb::AbstractCallback, sys, dvs, ps;
        expression = Val{true}, eval_expression = false, eval_module = @__MODULE__, kwargs...)
    u = map(value, dvs)
    p = map.(value, reorder_parameters(sys, ps))
    t = get_iv(sys)
    condit = conditions(cb)
    cs = collect_constants(condit)
    if !isempty(cs)
        cmap = map(x -> x => getdefault(x), cs)
        condit = substitute(condit, cmap)
    end
    expr = build_function_wrapper(sys,
        condit, u, t, p...; expression = Val{true},
        p_start = 3, p_end = length(p) + 2,
        wrap_code = condition_header(sys),
        kwargs...)
    if expression == Val{true}
        return expr
    end
    return eval_or_rgf(expr; eval_expression, eval_module)
end

"""
Compile user-defined functional affect.
"""
function compile_functional_affect(affect::FunctionalAffect, cb, sys, dvs, ps; kwargs...)
    dvs_ind = Dict(reverse(en) for en in enumerate(dvs))
    v_inds = map(sym -> dvs_ind[sym], unknowns(affect))

    if has_index_cache(sys) && get_index_cache(sys) !== nothing
        p_inds = [(pind = parameter_index(sys, sym)) === nothing ? sym : pind for sym in parameters(affect)]
        save_idxs = get(ic. callback_to_clocks, cb, Int[])
    else
        ps_ind = Dict(reverse(en) for en in enumerate(ps))
        p_inds = map(sym -> get(ps_ind, sym, sym), parameters(affect))
        save_idxs = Int[]
    end
    # HACK: filter out eliminated symbols. Not clear this is the right thing to do
    # (MTK should keep these symbols)
    u = filter(x -> !isnothing(x[2]), collect(zip(unknowns_syms(affect), v_inds))) |>
        NamedTuple
    p = filter(x -> !isnothing(x[2]), collect(zip(parameters_syms(affect), p_inds))) |>
        NamedTuple

    let u = u, p = p, user_affect = func(affect), ctx = context(affect),
        save_idxs = save_idxs

        function (integ)
            user_affect(integ, u, p, ctx)
            for idx in save_idxs
                SciMLBase.save_discretes!(integ, idx)
            end
        end
    end
end

is_discrete(cb::AbstractCallback) = cb isa SymbolicDiscreteCallback

"""
Codegen a DifferentialEquations callback. A set of continuous callbacks becomes a VectorContinuousCallback.
Individual callbacks become DiscreteCallback, PresetTimeCallback, PeriodicCallback, or ContinuousCallback
depending on the case.
"""
function generate_callback(cbs::Vector{SymbolicContinuousCallback}, sys; kwargs...)
    length(cbs) == 1 && return generate_callback(only(cbs), sys)
    eqs = map(cb -> flatten_equations(cb.eqs), cbs)

    _, f_iip = generate_custom_function(
        sys, [eq.lhs - eq.rhs for eq in eqs], unknowns(sys), parameters(sys); 
        expression = Val{false}, kwargs...)
    trigger = (out, u, t, integ) -> f_iip(out, u, parameter_values(integ), t)

    affects = []
    affect_negs = []
    inits = []
    finals = []
    for cb in cbs
        affect = compile_affect(cb.affect)

        push!(affects, affect)
        push!(affect_negs, compile_affect(cb.affect_neg, default = affect))
        push!(inits, compile_affect(cb.initialize, default = SciMLBase.INITALIZE_DEFAULT))
        push!(finals, compile_affect(cb.finalize, default = SciMLBase.FINALIZE_DEFAULT))
    end

    # Since there may be different number of conditions and affects,
    # we build a map that translates the condition eq. number to the affect number
    num_eqs = length.(eqs)
    eq2affect = reduce(vcat,
        [fill(i, num_eqs[i]) for i in eachindex(affects)])
    @assert length(eq2affect) == length(eqs)
    @assert maximum(eq2affect) == length(affect_functions)

    affect = function (integ, idx)
        affects[eq2affect[idx]](integ)
    end
    affect_neg = function (integ, idx)
        f = affect_negs[eq2affect[idx]]
        isnothing(f) && return
        f(integ)
    end
    initialize = compile_vector_optional_affect(inits, SciMLBase.INITIALIZE_DEFAULT)
    finalize = compile_vector_optional_affect(finals, SciMLBase.FINALIZE_DEFAULT)

    return VectorContinuousCallback(trigger, affect, length(cbs); affect_neg, initialize, finalize, rootfind = callback.rootfind, initializealg = SciMLBase.NoInit)
end

function generate_callback(cb, sys; kwargs...)
    is_timed = is_timed_condition(conditions(cb))
    dvs = unknowns(sys)
    ps = parameters(sys; initial_parameters = true)

    trigger = is_timed ? conditions(cb) : compile_condition(cb, sys, dvs, ps; kwargs...)
    affect = compile_affect(cb.affect) 
    affect_neg = hasfield(cb, :affect_neg) ? compile_affect(cb.affect_neg, default = affect) : nothing
    initialize = compile_affect(cb.initialize, default = SciMLBase.INITIALIZE_DEFAULT)
    finalize = compile_affect(cb.finalize, default = SciMLBase.FINALIZE_DEFAULT)

    if is_discrete(cb)
        if is_timed && condition(cb) isa AbstractVector
            return PresetTimeCallback(trigger, affect; affect_neg, initialize, finalize, initializealg = SciMLBase.NoInit)
        elseif is_timed
            return PeriodicCallback(affect, trigger; initialize, finalize)
        else
            return DiscreteCallback(trigger, affect; affect_neg, initialize, finalize, initializealg = SciMLBase.NoInit)
        end
    else
        return ContinuousCallback(trigger, affect; affect_neg, initialize, finalize, rootfind = cb.rootfind, initializealg = SciMLBase.NoInit)
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
function compile_affect(aff::Affect, cb::AbstractCallback, sys::AbstractSystem; default = nothing)
    save_idxs = if !(has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing)
        Int[]
    else
        get(ic.callback_to_clocks, cb, Int[])
    end

    isnothing(aff) && return default

    ps = parameters(aff; initial_parameters = true)
    dvs = unknowns(aff)

    if aff isa ImplicitDiscreteSystem
        function affect!(integrator) 
            pmap = [] 
            for pre_p in ps
                p = only(arguments(unwrap(pre_p)))
                push!(pmap, pre_p => integrator[p])
            end
            guesses = [u => integrator[u] for u in dvs]
            prob = ImplicitDiscreteProblem(aff, [], (0, 1), pmap; guesses)
            sol = init(prob, SimpleIDSolve())
            for u in dvs
                integrator[u] = sol[u]
            end

            for idx in save_idxs
                SciMLBase.save_discretes!(integ, idx)
            end
        end
    elseif aff isa FunctionalAffect || aff isa ImperativeAffect
        compile_functional_affect(aff, callback, sys, dvs, ps; kwargs...)
    end
end

"""
Initialize and Finalize for VectorContinuousCallback.
"""
function compile_vector_optional_affect(funs, default)
    all(isnothing, funs) && return default
    return let funs = funs
        function (cb, u, t, integ)
           for func in funs
               isnothing(func) ? continue : func(integ)
           end
        end
    end
end

merge_cb(::Nothing, ::Nothing) = nothing
merge_cb(::Nothing, x) = merge_cb(x, nothing)
merge_cb(x, ::Nothing) = x
merge_cb(x, y) = CallbackSet(x, y)

"""
Generate the CallbackSet for a ODESystem or SDESystem.
"""
function process_events(sys; callback = nothing, kwargs...)
    if has_continuous_events(sys) && !isempty(continuous_events(sys))
        cbs = continuous_events(sys)
        contin_cbs = generate_callback(cbs, sys; kwargs...)
    else
        contin_cbs = nothing
    end
    if has_discrete_events(sys) && !isempty(discrete_events(sys))
        dbs = discrete_events(sys)
        discrete_cbs = [generate_callback(db, sys; kwargs...) for db in dbs] 
    else
        discrete_cbs = nothing
    end

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
               (map(o -> namespace_callback(o, s), discrete_events(s)) for s in systems),
               init = SymbolicDiscreteCallback[])]
    cbs
end

has_discrete_events(sys::AbstractSystem) = isdefined(sys, :discrete_events)
function get_discrete_events(sys::AbstractSystem)
    has_discrete_events(sys) || return SymbolicDiscreteCallback[]
    getfield(sys, :discrete_events)
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
               (map(o -> namespace_callback(o, s), continuous_events(s))
               for s in systems),
               init = SymbolicContinuousCallback[])]
    filter(!isempty, cbs)
end

has_continuous_events(sys::AbstractSystem) = isdefined(sys, :continuous_events)
function get_continuous_events(sys::AbstractSystem)
    has_continuous_events(sys) || return SymbolicContinuousCallback[]
    getfield(sys, :continuous_events)
end
