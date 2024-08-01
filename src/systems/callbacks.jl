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
struct SymbolicContinuousCallback
    eqs::Vector{Equation}
    affect::Union{Vector{Equation}, FunctionalAffect}
    function SymbolicContinuousCallback(eqs::Vector{Equation}, affect = NULL_AFFECT)
        new(eqs, make_affect(affect))
    end # Default affect to nothing
end
make_affect(affect) = affect
make_affect(affect::Tuple) = FunctionalAffect(affect...)
make_affect(affect::NamedTuple) = FunctionalAffect(; affect...)

function Base.:(==)(e1::SymbolicContinuousCallback, e2::SymbolicContinuousCallback)
    isequal(e1.eqs, e2.eqs) && isequal(e1.affect, e2.affect)
end
Base.isempty(cb::SymbolicContinuousCallback) = isempty(cb.eqs)
function Base.hash(cb::SymbolicContinuousCallback, s::UInt)
    s = foldr(hash, cb.eqs, init = s)
    cb.affect isa AbstractVector ? foldr(hash, cb.affect, init = s) : hash(cb.affect, s)
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

namespace_affects(af::Vector, s) = Equation[namespace_affect(a, s) for a in af]
namespace_affects(af::FunctionalAffect, s) = namespace_affect(af, s)

function namespace_callback(cb::SymbolicContinuousCallback, s)::SymbolicContinuousCallback
    SymbolicContinuousCallback(namespace_equation.(equations(cb), (s,)),
        namespace_affects(affects(cb), s))
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

#################################### continuous events #####################################

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
        condit, u, t, p...; expression = Val{true},
        wrap_code = condition_header(sys) .∘ wrap_array_vars(sys, condit; dvs, ps),
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
                    update_inds = map(update_vars) do sym
                        return parameter_index(sys, sym)
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

        _ps = ps
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
            wrap_code = add_integrator_header(sys, integ, outvar) .∘
                        wrap_array_vars(sys, rhss; dvs, ps = _ps),
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

function generate_rootfinding_callback(cbs, sys::AbstractODESystem, dvs = unknowns(sys),
        ps = full_parameters(sys); kwargs...)
    eqs = map(cb -> flatten_equations(cb.eqs), cbs)
    num_eqs = length.(eqs)
    (isempty(eqs) || sum(num_eqs) == 0) && return nothing
    # fuse equations to create VectorContinuousCallback
    eqs = reduce(vcat, eqs)
    # rewrite all equations as 0 ~ interesting stuff
    eqs = map(eqs) do eq
        isequal(eq.lhs, 0) && return eq
        0 ~ eq.lhs - eq.rhs
    end

    rhss = map(x -> x.rhs, eqs)
    root_eq_vars = unique(collect(Iterators.flatten(map(ModelingToolkit.vars, rhss))))

    rf_oop, rf_ip = generate_custom_function(sys, rhss, dvs, ps; expression = Val{false},
        kwargs...)

    affect_functions = map(cbs) do cb # Keep affect function separate
        eq_aff = affects(cb)
        affect = compile_affect(eq_aff, sys, dvs, ps; expression = Val{false}, kwargs...)
    end

    if length(eqs) == 1
        cond = function (u, t, integ)
            if DiffEqBase.isinplace(integ.sol.prob)
                tmp, = DiffEqBase.get_tmp_cache(integ)
                rf_ip(tmp, u, parameter_values(integ), t)
                tmp[1]
            else
                rf_oop(u, parameter_values(integ), t)
            end
        end
        ContinuousCallback(cond, affect_functions[])
    else
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
                affect_functions[eq_ind2affect[eq_ind]](integ)
            end
        end
        VectorContinuousCallback(cond, affect, length(eqs))
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
