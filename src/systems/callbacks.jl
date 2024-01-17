#################################### system operations #####################################
get_continuous_events(sys::AbstractSystem) = Equation[]
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
    ctx::Any
end

function FunctionalAffect(f, sts, pars, ctx = nothing)
    # sts & pars contain either pairs: resistor.R => R, or Syms: R
    vs = [x isa Pair ? x.first : x for x in sts]
    vs_syms = Symbol[x isa Pair ? Symbol(x.second) : getname(x) for x in sts]
    length(vs_syms) == length(unique(vs_syms)) || error("Variables are not unique")

    ps = [x isa Pair ? x.first : x for x in pars]
    ps_syms = Symbol[x isa Pair ? Symbol(x.second) : getname(x) for x in pars]
    length(ps_syms) == length(unique(ps_syms)) || error("Parameters are not unique")

    FunctionalAffect(f, vs, vs_syms, ps, ps_syms, ctx)
end

FunctionalAffect(; f, sts, pars, ctx = nothing) = FunctionalAffect(f, sts, pars, ctx)

func(f::FunctionalAffect) = f.f
context(a::FunctionalAffect) = a.ctx
parameters(a::FunctionalAffect) = a.pars
parameters_syms(a::FunctionalAffect) = a.pars_syms
states(a::FunctionalAffect) = a.sts
states_syms(a::FunctionalAffect) = a.sts_syms

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
    hash(a.ctx, s)
end

has_functional_affect(cb) = affects(cb) isa FunctionalAffect

namespace_affect(affect, s) = namespace_equation(affect, s)
function namespace_affect(affect::FunctionalAffect, s)
    FunctionalAffect(func(affect),
        renamespace.((s,), states(affect)),
        states_syms(affect),
        renamespace.((s,), parameters(affect)),
        parameters_syms(affect),
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
SymbolicContinuousCallbacks(::Nothing) = SymbolicContinuousCallbacks(Equation[])

equations(cb::SymbolicContinuousCallback) = cb.eqs
function equations(cbs::Vector{<:SymbolicContinuousCallback})
    reduce(vcat, [equations(cb) for cb in cbs])
end

affects(cb::SymbolicContinuousCallback) = cb.affect
function affects(cbs::Vector{SymbolicContinuousCallback})
    reduce(vcat, [affects(cb) for cb in cbs])
end

namespace_affects(af::Vector, s) = Equation[namespace_affect(a, s) for a in af]
namespace_affects(af::FunctionalAffect, s) = namespace_affect(af, s)

function namespace_callback(cb::SymbolicContinuousCallback, s)::SymbolicContinuousCallback
    SymbolicContinuousCallback(namespace_equation.(equations(cb), (s,)),
        namespace_affects(affects(cb), s))
end

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

"Abstract super type to conveniently inject custom symbolic callbacks."
abstract type AbstractSpecialDiscreteCallback end

struct SymbolicDiscreteCallback{T <: Union{Nothing, AbstractSpecialDiscreteCallback}}
    # condition can be one of:
    #   Δt::Real - Periodic with period Δt
    #   Δts::Vector{Real} - events trigger in this times (Preset)
    #   condition::Vector{Equation} - event triggered when condition is true
    condition::Any
    affects::Any

    # The field `wrapped` allows to inject callbacks of custom symbolic callback types 
    # into the vector provided for the `discrete_events` keyword argument of system constructors.
    # By default, if some `sdcb:: AbstractSpecialDiscreteCallback` is encountered, then 
    # `condition` and `affects` are set to nothing by the outer constructor.
    # Therefore, `generate_discrete_callback` has to use `T` for dispatch and implement its 
    # own logic for compilation of `sdcb` into a `DiscreteCallback`.
    wrapped::T

    function SymbolicDiscreteCallback(condition, affects = NULL_AFFECT,
                                      wrapped::T = nothing) where {T}
        c = scalarize_condition(condition)
        a = scalarize_affects(affects)
        # TODO warn about effects of `wrapped` if its not `nothing`?
        new{T}(c, a, wrapped)
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
function SymbolicDiscreteCallback(sdcb::AbstractSpecialDiscreteCallback)
    SymbolicDiscreteCallback(nothing, nothing, sdcb)
end

function Base.show(io::IO, db::SymbolicDiscreteCallback)
    println(io, "condition: ", db.condition)
    println(io, "affects:")
    if db.affects isa FunctionalAffect || isnothing(db.affects)
        # TODO
        println(io, " ", db.affects)
    else
        for affect in db.affects
            println(io, "  ", affect)
        end
    end
    if !isnothing(db.wrapped)
        # TODO
        println(io, "wrapper for:")
        println(io, db.wrapped)
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
    reduce(vcat, affects(cb) for cb in cbs)
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
function SymbolicDiscreteCallbacks(sbc::AbstractSpecialDiscreteCallback)
    [SymbolicDiscreteCallback(sbc)]
end

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
function add_integrator_header(integrator = gensym(:MTKIntegrator), out = :u)
    expr -> Func([DestructuredArgs(expr.args, integrator, inds = [:u, :p, :t])], [],
        expr.body),
    expr -> Func([DestructuredArgs(expr.args, integrator, inds = [out, :u, :p, :t])], [],
        expr.body)
end

function condition_header(integrator = gensym(:MTKIntegrator))
    expr -> Func([expr.args[1], expr.args[2],
            DestructuredArgs(expr.args[3:end], integrator, inds = [:p])], [], expr.body)
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
        expression = Val{true}, kwargs...)
    u = map(x -> time_varying_as_func(value(x), sys), dvs)
    p = map(x -> time_varying_as_func(value(x), sys), ps)
    t = get_iv(sys)
    condit = condition(cb)
    cs = collect_constants(condit)
    if !isempty(cs)
        cmap = map(x -> x => getdefault(x), cs)
        condit = substitute(condit, cmap)
    end
    build_function(condit, u, t, p; expression, wrap_code = condition_header(),
        kwargs...)
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
    `states(sys)`. If provided, checks that the LHS of affect equations are variables are
    dropped, i.e. it is assumed these indices are correct and affect equations are
    well-formed.
  - `kwargs` are passed through to `Symbolics.build_function`.
"""
function compile_affect(eqs::Vector{Equation}, sys, dvs, ps; outputidxs = nothing,
        expression = Val{true}, checkvars = true,
        postprocess_affect_expr! = nothing, kwargs...)
    if isempty(eqs)
        if expression == Val{true}
            return :((args...) -> ())
        else
            return (args...) -> () # We don't do anything in the callback, we're just after the event
        end
    else
        rhss = map(x -> x.rhs, eqs)
        outvar = :u
        if outputidxs === nothing
            lhss = map(x -> x.lhs, eqs)
            all(isvariable, lhss) ||
                error("Non-variable symbolic expression found on the left hand side of an affect equation. Such equations must be of the form variable ~ symbolic expression for the new value of the variable.")
            update_vars = collect(Iterators.flatten(map(ModelingToolkit.vars, lhss))) # these are the ones we're changing
            length(update_vars) == length(unique(update_vars)) == length(eqs) ||
                error("affected variables not unique, each state can only be affected by one equation for a single `root_eqs => affects` pair.")
            alleq = all(isequal(isparameter(first(update_vars))),
                Iterators.map(isparameter, update_vars))
            if !isparameter(first(lhss)) && alleq
                stateind = Dict(reverse(en) for en in enumerate(dvs))
                update_inds = map(sym -> stateind[sym], update_vars)
            elseif isparameter(first(lhss)) && alleq
                psind = Dict(reverse(en) for en in enumerate(ps))
                update_inds = map(sym -> psind[sym], update_vars)
                outvar = :p
            else
                error("Error, building an affect function for a callback that wants to modify both parameters and states. This is not currently allowed in one individual callback.")
            end
        else
            update_inds = outputidxs
        end

        if checkvars
            u = map(x -> time_varying_as_func(value(x), sys), dvs)
            p = map(x -> time_varying_as_func(value(x), sys), ps)
        else
            u = dvs
            p = ps
        end
        t = get_iv(sys)
        integ = gensym(:MTKIntegrator)
        getexpr = (postprocess_affect_expr! === nothing) ? expression : Val{true}
        pre = get_preprocess_constants(rhss)
        rf_oop, rf_ip = build_function(rhss, u, p, t; expression = getexpr,
            wrap_code = add_integrator_header(integ, outvar),
            outputidxs = update_inds,
            postprocess_fbody = pre,
            kwargs...)
        # applied user-provided function to the generated expression
        if postprocess_affect_expr! !== nothing
            postprocess_affect_expr!(rf_ip, integ)
            (expression == Val{false}) &&
                (return drop_expr(@RuntimeGeneratedFunction(rf_ip)))
        end
        rf_ip
    end
end

function generate_rootfinding_callback(sys::AbstractODESystem, dvs = states(sys),
        ps = parameters(sys); kwargs...)
    cbs = continuous_events(sys)
    isempty(cbs) && return nothing
    generate_rootfinding_callback(cbs, sys, dvs, ps; kwargs...)
end

function generate_rootfinding_callback(cbs, sys::AbstractODESystem, dvs = states(sys),
        ps = parameters(sys); kwargs...)
    eqs = map(cb -> cb.eqs, cbs)
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

    u = map(x -> time_varying_as_func(value(x), sys), dvs)
    p = map(x -> time_varying_as_func(value(x), sys), ps)
    t = get_iv(sys)
    pre = get_preprocess_constants(rhss)
    rf_oop, rf_ip = build_function(rhss, u, p, t; expression = Val{false},
        postprocess_fbody = pre, kwargs...)

    affect_functions = map(cbs) do cb # Keep affect function separate
        eq_aff = affects(cb)
        affect = compile_affect(eq_aff, sys, dvs, ps; expression = Val{false}, kwargs...)
    end

    if length(eqs) == 1
        cond = function (u, t, integ)
            if DiffEqBase.isinplace(integ.sol.prob)
                tmp, = DiffEqBase.get_tmp_cache(integ)
                rf_ip(tmp, u, integ.p, t)
                tmp[1]
            else
                rf_oop(u, integ.p, t)
            end
        end
        ContinuousCallback(cond, affect_functions[])
    else
        cond = function (out, u, t, integ)
            rf_ip(out, u, integ.p, t)
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
    v_inds = map(sym -> dvs_ind[sym], states(affect))

    ps_ind = Dict(reverse(en) for en in enumerate(ps))
    p_inds = map(sym -> ps_ind[sym], parameters(affect))

    # HACK: filter out eliminated symbols. Not clear this is the right thing to do
    # (MTK should keep these symbols)
    u = filter(x -> !isnothing(x[2]), collect(zip(states_syms(affect), v_inds))) |>
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

function generate_discrete_callbacks(sys::AbstractSystem, dvs = states(sys),
        ps = parameters(sys); kwargs...)
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

function process_events(sys; callback = nothing, has_difference = false, kwargs...)
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
    difference_cb = has_difference ? generate_difference_cb(sys; kwargs...) : nothing

    cb = merge_cb(contin_cb, difference_cb)
    cb = merge_cb(cb, callback)
    (discrete_cb === nothing) ? cb : CallbackSet(cb, discrete_cb...)
end

################################# special discrete callbacks ####################################
"""
    SymbolicIterativeCallback(time_choice_tuple, user_affect!_tuple, initial_affect=false)

Define an iterative Callback based on the arguments `time_choice_tuple` and 
`user_affect!_tuple`.

# Arguments
- `time_choice_tuple` is a (optionally named) tuple of the form   
   `(time_choice, sts, pars, ctx)`.
- `user_affect!_tuple` is a tuple of the form `(user_affect!, sts, pars, ctx)`.
- `initial_affect::Bool=false` indicates whether or not to apply the affect at `t0`.

Both `time_choice` and `user_affect!` should have signatures 
`(integrator, sys_states, sys_params, ctx)`.
Within the respective function body, the value of a system state `v` can then be accessed by
`integrator.u[sys_states.v]`. It works the same for parameters. The time is `integrator.t`.
Both functions are compiled to suit the requirements of their counterparts for 
`DiffEqCallbacks.IterativeCallback`.

# Example
Below, we model the growth of some population with size ``N(t)``.
If we start with ``N(0) = 0``, then the population will approximate
``α`` in the limit, i.e., ``N(t) ↗ α`` for ``t → ∞``.
With some calculus, we find that ``N(t) ≈ 50`` at ``t = \\ln(2)``.
At this point in time, we decide to inject ``M = 60`` individuals
to the population. A short while later, we can check if that was 
successful:
```jldoctest; output=false
using ModelingToolkit, OrdinaryDiffEq

@parameters α=100 M=60
@variables t N(t)
Dt = Differential(t)
eqs = [Dt(N) ~ α - N]

t_inject = log(2)
t_check = t_inject + 1e-6
t_ref = Base.RefValue{Union{Nothing, Float64}}(-1.0)

function time_choice(integ, sts, pars, ctx)
    if t_ref[] < 0
        t_ref[] = t_inject
    elseif t_ref[] == t_inject
        t_ref[] = t_check
    else
        t_ref[] = nothing
    end
    return t_ref[]
end

function user_affect!(integ, sts, pars, ctx)
    if integ.t == t_inject
        integ.u[ sts.N ] += integ.p[ pars.M ]
    elseif integ.t == t_check
        @assert integ.u[ sts.N ] > integ.p[ pars.α ]
    end
end

# define discrete symbolic callback
cb = ModelingToolkit.SymbolicIterativeCallback(
    (time_choice, [], [], nothing),
    (user_affect!, [N], [α, M], nothing)
)

# setup system, problem, and solve:
@named osys = ODESystem(
    eqs, t, [N], [α, M]; discrete_events=[cb]
)
u0 = [N => 0.0]
tspan = (0.0, 20.0)
oprob = ODEProblem(osys, u0, tspan)
sol = solve(oprob, Tsit5())

# output

retcode: Success
Interpolation: specialized 4th order "free" interpolation
[...]
```

See also [`DiffEqCallbacks.IterativeCallback`](@ref).
"""
struct SymbolicIterativeCallback <: AbstractSpecialDiscreteCallback
    time_choice::Any
    user_affect!::Any

    initial_affect::Bool

    function SymbolicIterativeCallback(time_choice, user_affect!,
                                       initial_affect::Bool = false)
        _time_choice = scalarize_affects(time_choice)
        _user_affect! = scalarize_affects(user_affect!)
        return new(_time_choice, _user_affect!, initial_affect)
    end
end

function generate_discrete_callback(cb::SymbolicDiscreteCallback{T}, sys, dvs, ps;
                                    postprocess_affect_expr! = nothing,
                                    kwargs...) where {T <: SymbolicIterativeCallback}
    time_choice = compile_affect(cb.wrapped.time_choice, sys, dvs, ps;
                                 expression = Val{false},
                                 postprocess_affect_expr!, kwargs...)
    user_affect! = compile_affect(cb.wrapped.user_affect!, sys, dvs, ps;
                                  expression = Val{false}, postprocess_affect_expr!,
                                  kwargs...)
    return DiffEqCallbacks.IterativeCallback(time_choice, user_affect!;
                                             initial_affect = cb.wrapped.initial_affect)
end

"""
    SymbolicPeriodicCallback(
        Δt, user_affect!_tuple, initial_affect=false, final_affect=false
    )

Define a peridic Callback based on the arguments `Δt::Real` and `user_affect!_tuple`.
The affect is applied at times `tspan[1]`, `tspan[1] + Δt`, `tspan[1] + 2*Δt`, etc.

# Arguments
- `Δt::Real` is the periodicity of the event, i.e., the time between event occurrences.
- `user_affect!_tuple` is a tuple of the form `(user_affect!, sts, pars, ctx)`  
  `user_affect!` should have signatures `(integrator, sys_states, sys_params, ctx)`.  
   Within the respective function body, the value of a system state `v` can then be accessed  
   by `integrator.u[sys_states.v]`. It works the same for parameters.   
   The time is `integrator.t`. This function is compiled to suit the requirements of its  
   counterparts for `DiffEqCallbacks.PeriodicCallback`.
- `initial_affect::Bool=false` Whether or not to apply affect at `tspan[1]`.
- `final_affect::Bool=false` Whether or not to apply affect at `tspan[2]`.

# Example
Below, we model the growth of some population with size ``N(t)``.
If we start with ``N(0) = 0``, then the population will approximate
``α`` in the limit, i.e., ``N(t) ↗ α`` for ``t → ∞``.
But we now imagine some periodic event afflicting our population.
If at that point in time, the population is too big, it is reduced by 
substracting ``m=20`` via `user_affect!`:
```jldoctest; output=false
using ModelingToolkit, OrdinaryDiffEq

@parameters α=100 m=20
@variables t N(t)
Dt = Differential(t)
eqs = [Dt(N) ~ α - N]

del_t = log(2)/2

function user_affect!(integ, sts, pars, ctx)
    N = integ.u[ sts.N ]
    if N > 70
        integ.u[ sts.N ] -= integ.p[ pars.m ]
    end
end

# define discrete symbolic callback
cb = ModelingToolkit.SymbolicPeriodicCallback(
    del_t, (user_affect!, [N], [m], nothing)
)

# setup system, problem, and solve:
@named osys = ODESystem(
    eqs, t, [N], [α, m]; discrete_events=[cb]
)

u0 = [N => 0.0]
tspan = (0.0, 20.0)
oprob = ODEProblem(osys, u0, tspan)
sol = solve(oprob, Tsit5())

# output

retcode: Success
Interpolation: specialized 4th order "free" interpolation
[...]
```

See also [`DiffEqCallbacks.PeriodicCallback`](@ref).
"""
struct SymbolicPeriodicCallback <: AbstractSpecialDiscreteCallback
    del_t::Real
    affect::Any

    initial_affect::Bool
    final_affect::Bool
    function SymbolicPeriodicCallback(del_t, affect, initial_affect = false,
                                      final_affect = false)
        _affect = scalarize_affects(affect)
        return new(del_t, _affect, initial_affect, final_affect)
    end
end

function generate_discrete_callback(cb::SymbolicDiscreteCallback{T}, sys, dvs, ps;
                                    postprocess_affect_expr! = nothing,
                                    kwargs...) where {T <: SymbolicPeriodicCallback}
    spc = cb.wrapped
    aff = compile_affect(spc.affect, sys, dvs, ps;
                         expression = Val{false}, postprocess_affect_expr!, kwargs...)
    return PeriodicCallback(aff, spc.del_t;
                            initial_affect = spc.initial_affect,
                            final_affect = spc.final_affect)
end

"""
    SymbolicPresetTimeCallback(tstops, user_affect!_tuple, filter_tstops=true)

Define a Callback with preset times based on the arguments `tstops::Vector{Real}` and 
`user_affect!_tuple`. 

# Arguments
- `tstops::Vector{Real}` are the time stops for event occurrences.
- `user_affect!_tuple` is a tuple of the form `(user_affect!, sts, pars, ctx)`  
  `user_affect!` should have signatures `(integrator, sys_states, sys_params, ctx)`.  
   Within the respective function body, the value of a system state `v` can then be accessed  
   by `integrator.u[sys_states.v]`. It works the same for parameters.   
   The time is `integrator.t`. This function is compiled to suit the requirements of its  
   counterparts for `DiffEqCallbacks.PeriodicCallback`.
- `filter_tstops::Bool=true` Whether or not to filter out tstops beyond `tspan[2]`.

# Example
Below, we model the growth of some population with size ``N(t)``.
If we start with ``N(0) = 0``, then the population will approximate
``α`` in the limit, i.e., ``N(t) ↗ α`` for ``t → ∞``.
At preset times we want to check if ``N(t) > 70`.
If that is the case, we substract ``m=20``:
```jldoctest; output=false
using ModelingToolkit, OrdinaryDiffEq

@parameters α=100 m=20
@variables t N(t)
Dt = Differential(t)
eqs = [Dt(N) ~ α - N]

tspan = (0.0, 20.0)
tstops = LinRange(tspan[1], tspan[2], 100)

function user_affect!(integ, sts, pars, ctx)
    N = integ.u[ sts.N ]
    if N > 70
        integ.u[ sts.N ] -= integ.p[ pars.m ]
    end
end

# define discrete symbolic callback
cb = ModelingToolkit.SymbolicPresetTimeCallback(
    tstops, (user_affect!, [N], [m], nothing)
)

# setup system, problem, and solve:
@named osys = ODESystem(
    eqs, t, [N], [α, m]; discrete_events=[cb]
)
u0 = [N => 0.0]
oprob = ODEProblem(osys, u0, tspan)
sol = solve(oprob, Tsit5())

# output

retcode: Success
Interpolation: specialized 4th order "free" interpolation
[...]
```

See also [`DiffEqCallbacks.PresetTimeCallback`](@ref).
"""
struct SymbolicPresetTimeCallback <: AbstractSpecialDiscreteCallback
    tstops::Any # allow for iterators too (?)
    user_affect!::Any
    filter_tstops::Bool

    function SymbolicPresetTimeCallback(tstops, user_affect!, filter_tstops = true)
        _affect = scalarize_affects(user_affect!)
        return new(tstops, _affect, filter_tstops)
    end
end

function generate_discrete_callback(cb::SymbolicDiscreteCallback{T}, sys, dvs, ps;
                                    postprocess_affect_expr! = nothing,
                                    kwargs...) where {T <: SymbolicPresetTimeCallback}
    sptc = cb.wrapped
    aff = compile_affect(sptc.user_affect!, sys, dvs, ps;
                         expression = Val{false}, postprocess_affect_expr!, kwargs...)
    return DiffEqCallbacks.PresetTimeCallback(sptc.tstops, aff;
                                              filter_tstops = sptc.filter_tstops)
end
