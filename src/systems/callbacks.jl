#################################### system operations #####################################
get_continuous_events(sys::AbstractSystem) = SymbolicContinuousCallback[]
get_continuous_events(sys::AbstractODESystem) = getfield(sys, :continuous_events)
has_continuous_events(sys::AbstractSystem) = isdefined(sys, :continuous_events)

has_discrete_events(sys::AbstractSystem) = isdefined(sys, :discrete_events)
function get_discrete_events(sys::AbstractSystem)
    has_discrete_events(sys) || return SymbolicDiscreteCallback[]
    getfield(sys, :discrete_events)
end

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

func(f::FunctionalAffect) = f.f
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

has_functional_affect(cb) = affects(cb) isa FunctionalAffect

namespace_affect(affect, s) = namespace_equation(affect, s)
function namespace_affect(affect::FunctionalAffect, s)
    FunctionalAffect(func(affect),
        renamespace.((s,), unknowns(affect)),
        unknowns_syms(affect),
        renamespace.((s,), parameters(affect)),
        parameters_syms(affect),
        renamespace.((s,), discretes(affect)),
        context(affect))
end

#################################### continuous events #####################################

const NULL_AFFECT = Equation[]
"""
    SymbolicContinuousCallback(eqs::Vector{Equation}, affect, affect_neg, rootfind)

A [`ContinuousCallback`](@ref SciMLBase.ContinuousCallback) specified symbolically. Takes a vector of equations `eq`
as well as the positive-edge `affect` and negative-edge `affect_neg` that apply when *any* of `eq` are satisfied.
By default `affect_neg = affect`; to only get rising edges specify `affect_neg = nothing`.

Assume without loss of generality that the equation is of the form `c(u,p,t) ~ 0`; we denote the integrator state as `i.u`. 
For compactness, we define `prev_sign = sign(c(u[t-1], p[t-1], t-1))` and `cur_sign = sign(c(u[t], p[t], t))`.
A condition edge will be detected and the callback will be invoked iff `prev_sign * cur_sign <= 0`. 
Inter-sample condition activation is not guaranteed; for example if we use the dirac delta function as `c` to insert a 
sharp discontinuity between integrator steps (which in this example would not normally be identified by adaptivity) then the condition is not
guaranteed to be triggered.

Once detected the integrator will "wind back" through a root-finding process to identify the point when the condition became active; the method used 
is specified by `rootfind` from [`SciMLBase.RootfindOpt`](@ref). Multiple callbacks in the same system with different `rootfind` operations will be resolved 
into separate VectorContinuousCallbacks in the enumeration order of `SciMLBase.RootfindOpt`, which may cause some callbacks to not fire if several become
active at the same instant. See the `SciMLBase` documentation for more information on the semantic rules.

The positive edge `affect` will be triggered iff an edge is detected and if `prev_sign < 0`; similarly, `affect_neg` will be
triggered iff an edge is detected `prev_sign > 0`. 

Affects (i.e. `affect` and `affect_neg`) can be specified as either:
* A list of equations that should be applied when the callback is triggered (e.g. `x ~ 3, y ~ 7`) which must be of the form `unknown ~ observed value` where each `unknown` appears only once. Equations will be applied in the order that they appear in the vector; parameters and state updates will become immediately visible to following equations.
* A tuple `(f!, unknowns, read_parameters, modified_parameters, ctx)`, where:
    + `f!` is a function with signature `(integ, u, p, ctx)` that is called with the integrator, a state *index* vector `u` derived from `unknowns`, a parameter *index* vector `p` derived from `read_parameters`, and the `ctx` that was given at construction time. Note that `ctx` is aliased between instances.
    + `unknowns` is a vector of symbolic unknown variables and optionally their aliases (e.g. if the model was defined with `@variables x(t)` then a valid value for `unknowns` would be `[x]`). A variable can be aliased with a pair `x => :y`. The indices of these `unknowns` will be passed to `f!` in `u` in a named tuple; in the earlier example, if we pass `[x]` as `unknowns` then `f!` can access `x` as `integ.u[u.x]`. If no alias is specified the name of the index will be the symbol version of the variable name.
    + `read_parameters` is a vector of the parameters that are *used* by `f!`. Their indices are passed to `f` in `p` similarly to the indices of `unknowns` passed in `u`.
    + `modified_parameters` is a vector of the parameters that are *modified* by `f!`. Note that a parameter will not appear in `p` if it only appears in `modified_parameters`; it must appear in both `parameters` and `modified_parameters` if it is used in the affect definition.
    + `ctx` is a user-defined context object passed to `f!` when invoked. This value is aliased for each problem.
"""
struct SymbolicContinuousCallback
    eqs::Vector{Equation}
    affect::Union{Vector{Equation}, FunctionalAffect}
    affect_neg::Union{Vector{Equation}, FunctionalAffect, Nothing}
    rootfind::SciMLBase.RootfindOpt
    function SymbolicContinuousCallback(eqs::Vector{Equation}, affect = NULL_AFFECT,
            affect_neg = affect, rootfind = SciMLBase.LeftRootFind)
        new(eqs, make_affect(affect), make_affect(affect_neg), rootfind)
    end # Default affect to nothing
end
make_affect(affect) = affect
make_affect(affect::Tuple) = FunctionalAffect(affect...)
make_affect(affect::NamedTuple) = FunctionalAffect(; affect...)

function Base.:(==)(e1::SymbolicContinuousCallback, e2::SymbolicContinuousCallback)
    isequal(e1.eqs, e2.eqs) && isequal(e1.affect, e2.affect) &&
        isequal(e1.affect_neg, e2.affect_neg) && isequal(e1.rootfind, e2.rootfind)
end
Base.isempty(cb::SymbolicContinuousCallback) = isempty(cb.eqs)
function Base.hash(cb::SymbolicContinuousCallback, s::UInt)
    s = foldr(hash, cb.eqs, init = s)
    s = cb.affect isa AbstractVector ? foldr(hash, cb.affect, init = s) : hash(cb.affect, s)
    s = cb.affect_neg isa AbstractVector ? foldr(hash, cb.affect_neg, init = s) :
        hash(cb.affect_neg, s)
    hash(cb.rootfind, s)
end

to_equation_vector(eq::Equation) = [eq]
to_equation_vector(eqs::Vector{Equation}) = eqs
function to_equation_vector(eqs::Vector{Any})
    isempty(eqs) || error("This should never happen")
    Equation[]
end

function SymbolicContinuousCallback(args...)
    SymbolicContinuousCallback(to_equation_vector.(args)...)
end # wrap eq in vector
SymbolicContinuousCallback(p::Pair) = SymbolicContinuousCallback(p[1], p[2])
SymbolicContinuousCallback(cb::SymbolicContinuousCallback) = cb # passthrough
function SymbolicContinuousCallback(eqs::Equation, affect = NULL_AFFECT;
        affect_neg = affect, rootfind = SciMLBase.LeftRootFind)
    SymbolicContinuousCallback([eqs], affect, affect_neg, rootfind)
end
function SymbolicContinuousCallback(eqs::Vector{Equation}, affect;
        affect_neg = affect, rootfind = SciMLBase.LeftRootFind)
    SymbolicContinuousCallback(eqs, affect, affect_neg, rootfind)
end

SymbolicContinuousCallbacks(cb::SymbolicContinuousCallback) = [cb]
SymbolicContinuousCallbacks(cbs::Vector{<:SymbolicContinuousCallback}) = cbs
SymbolicContinuousCallbacks(cbs::Vector) = SymbolicContinuousCallback.(cbs)
function SymbolicContinuousCallbacks(ve::Vector{Equation})
    SymbolicContinuousCallbacks(SymbolicContinuousCallback(ve))
end
function SymbolicContinuousCallbacks(others)
    SymbolicContinuousCallbacks(SymbolicContinuousCallback(others))
end
SymbolicContinuousCallbacks(::Nothing) = SymbolicContinuousCallback[]

equations(cb::SymbolicContinuousCallback) = cb.eqs
function equations(cbs::Vector{<:SymbolicContinuousCallback})
    mapreduce(equations, vcat, cbs, init = Equation[])
end

affects(cb::SymbolicContinuousCallback) = cb.affect
function affects(cbs::Vector{SymbolicContinuousCallback})
    mapreduce(affects, vcat, cbs, init = Equation[])
end

affect_negs(cb::SymbolicContinuousCallback) = cb.affect_neg
function affect_negs(cbs::Vector{SymbolicContinuousCallback})
    mapreduce(affect_negs, vcat, cbs, init = Equation[])
end

namespace_affects(af::Vector, s) = Equation[namespace_affect(a, s) for a in af]
namespace_affects(af::FunctionalAffect, s) = namespace_affect(af, s)
namespace_affects(::Nothing, s) = nothing

function namespace_callback(cb::SymbolicContinuousCallback, s)::SymbolicContinuousCallback
    SymbolicContinuousCallback(
        namespace_equation.(equations(cb), (s,)),
        namespace_affects(affects(cb), s),
        namespace_affects(affect_negs(cb), s))
end

"""
    continuous_events(sys::AbstractSystem)::Vector{SymbolicContinuousCallback}

Returns a vector of all the `continuous_events` in an abstract system and its component subsystems.
The `SymbolicContinuousCallback`s in the returned vector are structs with two fields: `eqs` and
`affect` which correspond to the first and second elements of a `Pair` used to define an event, i.e.
`eqs => affect`.
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

#################################### discrete events #####################################

struct SymbolicDiscreteCallback
    # condition can be one of:
    #   Δt::Real - Periodic with period Δt
    #   Δts::Vector{Real} - events trigger in this times (Preset)
    #   condition::Vector{Equation} - event triggered when condition is true
    # TODO: Iterative
    condition::Any
    affects::Any

    function SymbolicDiscreteCallback(condition, affects = NULL_AFFECT)
        c = scalarize_condition(condition)
        a = scalarize_affects(affects)
        new(c, a)
    end # Default affect to nothing
end

is_timed_condition(cb) = false
is_timed_condition(::R) where {R <: Real} = true
is_timed_condition(::V) where {V <: AbstractVector} = eltype(V) <: Real
is_timed_condition(::Num) = false
is_timed_condition(cb::SymbolicDiscreteCallback) = is_timed_condition(condition(cb))

function scalarize_condition(condition)
    is_timed_condition(condition) ? condition : value(scalarize(condition))
end
function namespace_condition(condition, s)
    is_timed_condition(condition) ? condition : namespace_expr(condition, s)
end

scalarize_affects(affects) = scalarize(affects)
scalarize_affects(affects::Tuple) = FunctionalAffect(affects...)
scalarize_affects(affects::NamedTuple) = FunctionalAffect(; affects...)
scalarize_affects(affects::FunctionalAffect) = affects

SymbolicDiscreteCallback(p::Pair) = SymbolicDiscreteCallback(p[1], p[2])
SymbolicDiscreteCallback(cb::SymbolicDiscreteCallback) = cb # passthrough

function Base.show(io::IO, db::SymbolicDiscreteCallback)
    println(io, "condition: ", db.condition)
    println(io, "affects:")
    if db.affects isa FunctionalAffect
        # TODO
        println(io, " ", db.affects)
    else
        for affect in db.affects
            println(io, "  ", affect)
        end
    end
end

function Base.:(==)(e1::SymbolicDiscreteCallback, e2::SymbolicDiscreteCallback)
    isequal(e1.condition, e2.condition) && isequal(e1.affects, e2.affects)
end
function Base.hash(cb::SymbolicDiscreteCallback, s::UInt)
    s = hash(cb.condition, s)
    cb.affects isa AbstractVector ? foldr(hash, cb.affects, init = s) : hash(cb.affects, s)
end

condition(cb::SymbolicDiscreteCallback) = cb.condition
function conditions(cbs::Vector{<:SymbolicDiscreteCallback})
    reduce(vcat, condition(cb) for cb in cbs)
end

affects(cb::SymbolicDiscreteCallback) = cb.affects

function affects(cbs::Vector{SymbolicDiscreteCallback})
    reduce(vcat, affects(cb) for cb in cbs; init = [])
end

function namespace_callback(cb::SymbolicDiscreteCallback, s)::SymbolicDiscreteCallback
    af = affects(cb)
    af = af isa AbstractVector ? namespace_affect.(af, Ref(s)) : namespace_affect(af, s)
    SymbolicDiscreteCallback(namespace_condition(condition(cb), s), af)
end

SymbolicDiscreteCallbacks(cb::Pair) = SymbolicDiscreteCallback[SymbolicDiscreteCallback(cb)]
SymbolicDiscreteCallbacks(cbs::Vector) = SymbolicDiscreteCallback.(cbs)
SymbolicDiscreteCallbacks(cb::SymbolicDiscreteCallback) = [cb]
SymbolicDiscreteCallbacks(cbs::Vector{<:SymbolicDiscreteCallback}) = cbs
SymbolicDiscreteCallbacks(::Nothing) = SymbolicDiscreteCallback[]

"""
    discrete_events(sys::AbstractSystem) :: Vector{SymbolicDiscreteCallback}

Returns a vector of all the `discrete_events` in an abstract system and its component subsystems.
The `SymbolicDiscreteCallback`s in the returned vector are structs with two fields: `condition` and
`affect` which correspond to the first and second elements of a `Pair` used to define an event, i.e.
`condition => affect`.
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

################################# compilation functions ####################################

# handles ensuring that affect! functions work with integrator arguments
function add_integrator_header(
        sys::AbstractSystem, integrator = gensym(:MTKIntegrator), out = :u)
    if has_index_cache(sys) && get_index_cache(sys) !== nothing
        function (expr)
            p = gensym(:p)
            Func(
                [
                    DestructuredArgs([expr.args[1], p, expr.args[end]],
                    integrator, inds = [:u, :p, :t])
                ],
                [],
                Let(
                    [DestructuredArgs([arg.name for arg in expr.args[2:(end - 1)]], p),
                        expr.args[2:(end - 1)]...],
                    expr.body,
                    false)
            )
        end,
        function (expr)
            p = gensym(:p)
            Func(
                [
                    DestructuredArgs([expr.args[1], expr.args[2], p, expr.args[end]],
                    integrator, inds = [out, :u, :p, :t])
                ],
                [],
                Let(
                    [DestructuredArgs([arg.name for arg in expr.args[3:(end - 1)]], p),
                        expr.args[3:(end - 1)]...],
                    expr.body,
                    false)
            )
        end
    else
        expr -> Func([DestructuredArgs(expr.args, integrator, inds = [:u, :p, :t])], [],
            expr.body),
        expr -> Func(
            [DestructuredArgs(expr.args, integrator, inds = [out, :u, :p, :t])], [],
            expr.body)
    end
end

function condition_header(sys::AbstractSystem, integrator = gensym(:MTKIntegrator))
    if has_index_cache(sys) && get_index_cache(sys) !== nothing
        function (expr)
            p = gensym(:p)
            res = Func(
                [expr.args[1], expr.args[2],
                    DestructuredArgs([p], integrator, inds = [:p])],
                [],
                Let(
                    [
                        DestructuredArgs([arg.name for arg in expr.args[3:end]], p),
                        expr.args[3:end]...
                    ], expr.body, false
                )
            )
            return res
        end
    else
        expr -> Func(
            [expr.args[1], expr.args[2],
                DestructuredArgs(expr.args[3:end], integrator, inds = [:p])],
            [],
            expr.body)
    end
end

"""
    compile_condition(cb::SymbolicDiscreteCallback, sys, dvs, ps; expression, kwargs...)

Returns a function `condition(u,p,t)` returning the `condition(cb)`.

Notes

  - `expression = Val{true}`, causes the generated function to be returned as an expression.
    If  set to `Val{false}` a `RuntimeGeneratedFunction` will be returned.
  - `kwargs` are passed through to `Symbolics.build_function`.
"""
function compile_condition(cb::SymbolicDiscreteCallback, sys, dvs, ps;
        expression = Val{true}, eval_expression = false, eval_module = @__MODULE__, kwargs...)
    u = map(x -> time_varying_as_func(value(x), sys), dvs)
    p = map.(x -> time_varying_as_func(value(x), sys), reorder_parameters(sys, ps))
    t = get_iv(sys)
    condit = condition(cb)
    cs = collect_constants(condit)
    if !isempty(cs)
        cmap = map(x -> x => getdefault(x), cs)
        condit = substitute(condit, cmap)
    end
    expr = build_function(
        condit, u, t, p...; expression = Val{true}, wrap_code = condition_header(sys),
        kwargs...)
    if expression == Val{true}
        return expr
    end
    return eval_or_rgf(expr; eval_expression, eval_module)
end

function compile_affect(cb::SymbolicContinuousCallback, args...; kwargs...)
    compile_affect(affects(cb), args...; kwargs...)
end

"""
    compile_affect(eqs::Vector{Equation}, sys, dvs, ps; expression, outputidxs, kwargs...)
    compile_affect(cb::SymbolicContinuousCallback, args...; kwargs...)

Returns a function that takes an integrator as argument and modifies the state with the
affect. The generated function has the signature `affect!(integrator)`.

Notes

  - `expression = Val{true}`, causes the generated function to be returned as an expression.
    If  set to `Val{false}` a `RuntimeGeneratedFunction` will be returned.
  - `outputidxs`, a vector of indices of the output variables which should correspond to
    `unknowns(sys)`. If provided, checks that the LHS of affect equations are variables are
    dropped, i.e. it is assumed these indices are correct and affect equations are
    well-formed.
  - `kwargs` are passed through to `Symbolics.build_function`.
"""
function compile_affect(eqs::Vector{Equation}, sys, dvs, ps; outputidxs = nothing,
        expression = Val{true}, checkvars = true, eval_expression = false,
        eval_module = @__MODULE__,
        postprocess_affect_expr! = nothing, kwargs...)
    if isempty(eqs)
        if expression == Val{true}
            return :((args...) -> ())
        else
            return (args...) -> () # We don't do anything in the callback, we're just after the event
        end
    else
        eqs = flatten_equations(eqs)
        rhss = map(x -> x.rhs, eqs)
        outvar = :u
        if outputidxs === nothing
            lhss = map(x -> x.lhs, eqs)
            all(isvariable, lhss) ||
                error("Non-variable symbolic expression found on the left hand side of an affect equation. Such equations must be of the form variable ~ symbolic expression for the new value of the variable.")
            update_vars = collect(Iterators.flatten(map(ModelingToolkit.vars, lhss))) # these are the ones we're changing
            length(update_vars) == length(unique(update_vars)) == length(eqs) ||
                error("affected variables not unique, each unknown can only be affected by one equation for a single `root_eqs => affects` pair.")
            alleq = all(isequal(isparameter(first(update_vars))),
                Iterators.map(isparameter, update_vars))
            if !isparameter(first(lhss)) && alleq
                unknownind = Dict(reverse(en) for en in enumerate(dvs))
                update_inds = map(sym -> unknownind[sym], update_vars)
            elseif isparameter(first(lhss)) && alleq
                if has_index_cache(sys) && get_index_cache(sys) !== nothing
                    ic = get_index_cache(sys)
                    update_inds = map(update_vars) do sym
                        pind = parameter_index(sys, sym)
                        discrete_linear_index(ic, pind)
                    end
                else
                    psind = Dict(reverse(en) for en in enumerate(ps))
                    update_inds = map(sym -> psind[sym], update_vars)
                end
                outvar = :p
            else
                error("Error, building an affect function for a callback that wants to modify both parameters and unknowns. This is not currently allowed in one individual callback.")
            end
        else
            update_inds = outputidxs
        end

        ps = reorder_parameters(sys, ps)
        if checkvars
            u = map(x -> time_varying_as_func(value(x), sys), dvs)
            p = map.(x -> time_varying_as_func(value(x), sys), ps)
        else
            u = dvs
            p = ps
        end
        t = get_iv(sys)
        integ = gensym(:MTKIntegrator)
        pre = get_preprocess_constants(rhss)
        rf_oop, rf_ip = build_function(rhss, u, p..., t; expression = Val{true},
            wrap_code = add_integrator_header(sys, integ, outvar),
            outputidxs = update_inds,
            postprocess_fbody = pre,
            kwargs...)
        # applied user-provided function to the generated expression
        if postprocess_affect_expr! !== nothing
            postprocess_affect_expr!(rf_ip, integ)
        end
        if expression == Val{false}
            return eval_or_rgf(rf_ip; eval_expression, eval_module)
        end
        return rf_ip
    end
end

function generate_rootfinding_callback(sys::AbstractODESystem, dvs = unknowns(sys),
        ps = full_parameters(sys); kwargs...)
    cbs = continuous_events(sys)
    isempty(cbs) && return nothing
    generate_rootfinding_callback(cbs, sys, dvs, ps; kwargs...)
end
"""
Generate a single rootfinding callback; this happens if there is only one equation in `cbs` passed to 
generate_rootfinding_callback and thus we can produce a ContinuousCallback instead of a VectorContinuousCallback.
"""
function generate_single_rootfinding_callback(
        eq, cb, sys::AbstractODESystem, dvs = unknowns(sys),
        ps = full_parameters(sys); kwargs...)
    if !isequal(eq.lhs, 0)
        eq = 0 ~ eq.lhs - eq.rhs
    end

    rf_oop, rf_ip = generate_custom_function(
        sys, [eq.rhs], dvs, ps; expression = Val{false}, kwargs...)
    affect_function = compile_affect_fn(cb, sys, dvs, ps, kwargs)
    cond = function (u, t, integ)
        if DiffEqBase.isinplace(integ.sol.prob)
            tmp, = DiffEqBase.get_tmp_cache(integ)
            rf_ip(tmp, u, parameter_values(integ), t)
            tmp[1]
        else
            rf_oop(u, parameter_values(integ), t)
        end
    end
    return ContinuousCallback(
        cond, affect_function.affect, affect_function.affect_neg, rootfind = cb.rootfind)
end

function generate_vector_rootfinding_callback(
        cbs, sys::AbstractODESystem, dvs = unknowns(sys),
        ps = full_parameters(sys); rootfind = SciMLBase.RightRootFind, kwargs...)
    eqs = map(cb -> flatten_equations(cb.eqs), cbs)
    num_eqs = length.(eqs)
    # fuse equations to create VectorContinuousCallback
    eqs = reduce(vcat, eqs)
    # rewrite all equations as 0 ~ interesting stuff
    eqs = map(eqs) do eq
        isequal(eq.lhs, 0) && return eq
        0 ~ eq.lhs - eq.rhs
    end

    rhss = map(x -> x.rhs, eqs)
    _, rf_ip = generate_custom_function(
        sys, rhss, dvs, ps; expression = Val{false}, kwargs...)

    affect_functions = @NamedTuple{affect::Function, affect_neg::Union{Function, Nothing}}[compile_affect_fn(
                                                                                               cb,
                                                                                               sys,
                                                                                               dvs,
                                                                                               ps,
                                                                                               kwargs)
                                                                                           for cb in cbs]
    cond = function (out, u, t, integ)
        rf_ip(out, u, parameter_values(integ), t)
    end

    # since there may be different number of conditions and affects,
    # we build a map that translates the condition eq. number to the affect number
    eq_ind2affect = reduce(vcat,
        [fill(i, num_eqs[i]) for i in eachindex(affect_functions)])
    @assert length(eq_ind2affect) == length(eqs)
    @assert maximum(eq_ind2affect) == length(affect_functions)

    affect = let affect_functions = affect_functions, eq_ind2affect = eq_ind2affect
        function (integ, eq_ind) # eq_ind refers to the equation index that triggered the event, each event has num_eqs[i] equations
            affect_functions[eq_ind2affect[eq_ind]].affect(integ)
        end
    end
    affect_neg = let affect_functions = affect_functions, eq_ind2affect = eq_ind2affect
        function (integ, eq_ind) # eq_ind refers to the equation index that triggered the event, each event has num_eqs[i] equations
            affect_neg = affect_functions[eq_ind2affect[eq_ind]].affect_neg
            if isnothing(affect_neg)
                return # skip if the neg function doesn't exist - don't want to split this into a separate VCC because that'd break ordering
            end
            affect_neg(integ)
        end
    end
    return VectorContinuousCallback(
        cond, affect, affect_neg, length(eqs), rootfind = rootfind)
end

"""
Compile a single continuous callback affect function(s).
"""
function compile_affect_fn(cb, sys::AbstractODESystem, dvs, ps, kwargs)
    eq_aff = affects(cb)
    eq_neg_aff = affect_negs(cb)
    affect = compile_affect(eq_aff, sys, dvs, ps; expression = Val{false}, kwargs...)
    if eq_neg_aff === eq_aff
        affect_neg = affect
    elseif isnothing(eq_neg_aff)
        affect_neg = nothing
    else
        affect_neg = compile_affect(
            eq_neg_aff, sys, dvs, ps; expression = Val{false}, kwargs...)
    end
    (affect = affect, affect_neg = affect_neg)
end

function generate_rootfinding_callback(cbs, sys::AbstractODESystem, dvs = unknowns(sys),
        ps = full_parameters(sys); kwargs...)
    eqs = map(cb -> flatten_equations(cb.eqs), cbs)
    num_eqs = length.(eqs)
    total_eqs = sum(num_eqs)
    (isempty(eqs) || total_eqs == 0) && return nothing
    if total_eqs == 1
        # find the callback with only one eq
        cb_ind = findfirst(>(0), num_eqs)
        if isnothing(cb_ind)
            error("Inconsistent state in affect compilation; one equation but no callback with equations?")
        end
        cb = cbs[cb_ind]
        return generate_single_rootfinding_callback(cb.eqs[], cb, sys, dvs, ps; kwargs...)
    end

    # group the cbs by what rootfind op they use
    # groupby would be very useful here, but alas
    cb_classes = Dict{
        @NamedTuple{rootfind::SciMLBase.RootfindOpt}, Vector{SymbolicContinuousCallback}}()
    for cb in cbs
        push!(
            get!(() -> SymbolicContinuousCallback[], cb_classes, (rootfind = cb.rootfind,)),
            cb)
    end

    # generate the callbacks out; we sort by the equivalence class to ensure a deterministic preference order
    compiled_callbacks = map(collect(pairs(sort!(
        OrderedDict(cb_classes); by = p -> p.rootfind)))) do (equiv_class, cbs_in_class)
        return generate_vector_rootfinding_callback(
            cbs_in_class, sys, dvs, ps; rootfind = equiv_class.rootfind, kwargs...)
    end
    if length(compiled_callbacks) == 1
        return compiled_callbacks[]
    else
        return CallbackSet(compiled_callbacks...)
    end
end

function compile_user_affect(affect::FunctionalAffect, sys, dvs, ps; kwargs...)
    dvs_ind = Dict(reverse(en) for en in enumerate(dvs))
    v_inds = map(sym -> dvs_ind[sym], unknowns(affect))

    if has_index_cache(sys) && get_index_cache(sys) !== nothing
        p_inds = [parameter_index(sys, sym) for sym in parameters(affect)]
    else
        ps_ind = Dict(reverse(en) for en in enumerate(ps))
        p_inds = map(sym -> ps_ind[sym], parameters(affect))
    end
    # HACK: filter out eliminated symbols. Not clear this is the right thing to do
    # (MTK should keep these symbols)
    u = filter(x -> !isnothing(x[2]), collect(zip(unknowns_syms(affect), v_inds))) |>
        NamedTuple
    p = filter(x -> !isnothing(x[2]), collect(zip(parameters_syms(affect), p_inds))) |>
        NamedTuple

    let u = u, p = p, user_affect = func(affect), ctx = context(affect)
        function (integ)
            user_affect(integ, u, p, ctx)
        end
    end
end

function compile_affect(affect::FunctionalAffect, sys, dvs, ps; kwargs...)
    compile_user_affect(affect, sys, dvs, ps; kwargs...)
end

function generate_timed_callback(cb, sys, dvs, ps; postprocess_affect_expr! = nothing,
        kwargs...)
    cond = condition(cb)
    as = compile_affect(affects(cb), sys, dvs, ps; expression = Val{false},
        postprocess_affect_expr!, kwargs...)
    if cond isa AbstractVector
        # Preset Time
        return PresetTimeCallback(cond, as)
    else
        # Periodic
        return PeriodicCallback(as, cond)
    end
end

function generate_discrete_callback(cb, sys, dvs, ps; postprocess_affect_expr! = nothing,
        kwargs...)
    if is_timed_condition(cb)
        return generate_timed_callback(cb, sys, dvs, ps; postprocess_affect_expr!,
            kwargs...)
    else
        c = compile_condition(cb, sys, dvs, ps; expression = Val{false}, kwargs...)
        as = compile_affect(affects(cb), sys, dvs, ps; expression = Val{false},
            postprocess_affect_expr!, kwargs...)
        return DiscreteCallback(c, as)
    end
end

function generate_discrete_callbacks(sys::AbstractSystem, dvs = unknowns(sys),
        ps = full_parameters(sys); kwargs...)
    has_discrete_events(sys) || return nothing
    symcbs = discrete_events(sys)
    isempty(symcbs) && return nothing

    dbs = map(symcbs) do cb
        generate_discrete_callback(cb, sys, dvs, ps; kwargs...)
    end

    dbs
end

merge_cb(::Nothing, ::Nothing) = nothing
merge_cb(::Nothing, x) = merge_cb(x, nothing)
merge_cb(x, ::Nothing) = x
merge_cb(x, y) = CallbackSet(x, y)

function process_events(sys; callback = nothing, kwargs...)
    if has_continuous_events(sys)
        contin_cb = generate_rootfinding_callback(sys; kwargs...)
    else
        contin_cb = nothing
    end
    if has_discrete_events(sys)
        discrete_cb = generate_discrete_callbacks(sys; kwargs...)
    else
        discrete_cb = nothing
    end

    cb = merge_cb(contin_cb, callback)
    (discrete_cb === nothing) ? cb : CallbackSet(contin_cb, discrete_cb...)
end
