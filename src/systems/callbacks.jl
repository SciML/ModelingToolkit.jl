abstract type AbstractCallback end
function has_functional_affect(cb)
end
struct SymbolicAffect
    affect::Vector{Equation}
    alg_eqs::Vector{Equation}
    discrete_parameters::Vector{SymbolicT}
end
function SymbolicAffect(affect::Vector{Equation}; alg_eqs = Equation[],
        discrete_parameters = SymbolicT[], kwargs...)
    if symbolic_type(discrete_parameters) !== NotSymbolic()
        for p in discrete_parameters
        end
    end
    SymbolicAffect(affect, alg_eqs, discrete_parameters)
end
SymbolicAffect(affect; kwargs...) = make_affect(affect; kwargs...)
function (s::SymbolicUtils.Substituter)(aff::SymbolicAffect)
end
struct AffectSystem
    system::AbstractSystem
    unknowns::Vector{SymbolicT}
    parameters::Vector{SymbolicT}
    discretes::Vector{SymbolicT}
end
function AffectSystem(spec::SymbolicAffect; iv = nothing, alg_eqs = Equation[], kwargs...)
    AffectSystem(spec.affect; alg_eqs = vcat(spec.alg_eqs, alg_eqs), iv,
        discrete_parameters = spec.discrete_parameters, kwargs...)
end
@noinline function warn_algebraic_equation(eq::Equation)
end
function AffectSystem(affect::Vector{Equation}; discrete_parameters = SymbolicT[],
        iv = nothing, alg_eqs::Vector{Equation} = Equation[], warn_no_algebraic = true, kwargs...)
    isempty(affect) && return nothing
    for p in discrete_parameters
    end
    dvs = OrderedSet{SymbolicT}()
    params = OrderedSet{SymbolicT}()
    for eq in affect
        if !haspre(eq) && !(isconst(eq.lhs) && isconst(eq.rhs))
        end
        collect_vars!(dvs, params, eq, iv)
    end
    pre_params = filter(haspre, params)
    sys_params = SymbolicT[]
    discretes = map(tovar, discrete_parameters)
    dvs = collect(dvs)
    _dvs = map(default_toterm, dvs)
    affectsys = System(
        vcat(affect, alg_eqs), iv, collect(union(_dvs, discretes)),
        collect(union(pre_params, sys_params)); is_discrete = true, name = :affectsys)
    affectsys = mtkcompile(affectsys; fully_determined = nothing)
    accessed_params = Vector{SymbolicT}(filter(isparameter, map(unPre, collect(pre_params))))
    AffectSystem(affectsys, _dvs, accessed_params, discrete_parameters)
end
function Base.show(iio::IO, aff::AffectSystem)
end
function vars!(vars, aff::AffectSystem; op = Differential)
    for var in Iterators.flatten((unknowns(aff), parameters(aff), discretes(aff)))
    end
end
struct Pre <: Symbolics.Operator end
Pre(x) = Pre()(x)
unPre(x::SymbolicT) = (iscall(x) && operation(x) isa Pre) ? only(arguments(x)) : x
function (p::Pre)(x)
    iw = Symbolics.iswrapped(x)
end
haspre(O) = recursive_hasoperator(Pre, O)
const Affect = Union{AffectSystem}
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
        rootfind = SciMLBase.LeftRootFind,
        reinitializealg = nothing,
        zero_crossing_id = gensym(),
        kwargs...)
    if isnothing(reinitializealg)
        reinitializealg = SciMLBase.NoInit()
    end
    SymbolicContinuousCallback(conditions, SymbolicAffect(affect; kwargs...),
        SymbolicAffect(affect_neg; kwargs...),
        SymbolicAffect(initialize; kwargs...), SymbolicAffect(
            finalize; kwargs...),
        rootfind, reinitializealg, zero_crossing_id)
end
function SymbolicContinuousCallback(p::Pair, args...; kwargs...)
    SymbolicContinuousCallback(p[1], p[2], args...; kwargs...)
end
function SymbolicContinuousCallback(cb::Tuple, args...; kwargs...)
    if length(cb) == 2
    end
end
function complete(cb::SymbolicContinuousCallback; kwargs...)
    SymbolicContinuousCallback(cb.conditions, make_affect(cb.affect; kwargs...),
        make_affect(cb.finalize; kwargs...), cb.rootfind, cb.reinitializealg, cb.zero_crossing_id)
end
make_affect(affect::SymbolicAffect; kwargs...) = AffectSystem(affect; kwargs...)
function make_affect(affect; kwargs...)
end
struct SymbolicDiscreteCallback <: AbstractCallback
    conditions::Union{Number, Vector{<:Number}, SymbolicT}
end
function SymbolicDiscreteCallback(
        reinitializealg = nothing, kwargs...)
    if isnothing(reinitializealg)
    end
    SymbolicDiscreteCallback(c, SymbolicAffect(affect; kwargs...),
        SymbolicAffect(finalize; kwargs...), reinitializealg)
end
function SymbolicDiscreteCallback(p::Pair, args...; kwargs...)
end
function SymbolicDiscreteCallback(cb::Tuple, args...; kwargs...)
    if length(cb) == 2
    end
end
function complete(cb::SymbolicDiscreteCallback; kwargs...)
    SymbolicDiscreteCallback(cb.conditions, make_affect(cb.affect; kwargs...),
        make_affect(cb.finalize; kwargs...), cb.reinitializealg)
    if T === Num
    end
end
to_cb_vector(cbs::Vector{<:AbstractCallback}; kwargs...) = cbs
function to_cb_vector(cbs; CB_TYPE = SymbolicContinuousCallback, kwargs...)
    if cbs isa Pair
        [CB_TYPE(cbs; kwargs...)]
    end
end
function namespace_affects(affect::AffectSystem, s)
    affsys = system(affect)
    affsys = System(Equation[], get_iv(affsys); systems = [affsys], name = :affectsys)
    affsys = complete(affsys)
    AffectSystem(affsys,
        renamespace.((s,), affect.discrete_parameters))
end
function namespace_callback(cb::SymbolicContinuousCallback, s)::SymbolicContinuousCallback
    SymbolicContinuousCallback(
        affect_neg = namespace_affects(affect_negs(cb), s),
        expr.body)
end
function default_operating_point(affsys::AffectSystem)
    for p in parameters(sys)
        if T <: Number
        end
    end
end
function compile_equational_affect(
        eval_expression = false, eval_module = @__MODULE__, op = nothing, kwargs...)
    if aff isa AbstractVector
        p_idxs = if wrap_mtkparameters
        end
        p_up! = build_function_wrapper(sys, (@view rhss[is_p]), dvs, _ps..., t;
            eval_module)
        return let dvs_to_update = dvs_to_update, ps_to_update = ps_to_update,
            reset_jumps = reset_jumps, u_up! = u_up!, p_up! = p_up!
            function explicit_affect!(integ)
            end
        end
    end
end
has_discrete_events(sys::AbstractSystem) = isdefined(sys, :discrete_events)
function get_discrete_events(sys::AbstractSystem)
    has_discrete_events(sys) || return SymbolicDiscreteCallback[]
end
function discrete_events_toplevel(sys::AbstractSystem)
    if has_parent(sys) && (parent = get_parent(sys)) !== nothing
    end
end
function continuous_events(sys::AbstractSystem)
    obs = get_continuous_events(sys)
    systems = get_systems(sys)
    cbs = [obs;
           reduce(vcat,
               (map(o -> namespace_callback(o, s), continuous_events(s)) for s in systems),
               init = SymbolicContinuousCallback[])]
end
has_continuous_events(sys::AbstractSystem) = isdefined(sys, :continuous_events)
function get_continuous_events(sys::AbstractSystem)
    getfield(sys, :continuous_events)
end
function continuous_events_toplevel(sys::AbstractSystem)
    if has_parent(sys) && (parent = get_parent(sys)) !== nothing
    end
end
function create_symbolic_events(cont_events, disc_events)
    cont_callbacks = to_cb_vector(cont_events; CB_TYPE = SymbolicContinuousCallback)
    disc_callbacks = to_cb_vector(disc_events; CB_TYPE = SymbolicDiscreteCallback)
    cont_callbacks, disc_callbacks
end
