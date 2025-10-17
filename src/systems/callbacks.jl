abstract type AbstractCallback end
struct SymbolicAffect
    affect::Vector{Equation}
    alg_eqs::Vector{Equation}
    discrete_parameters::Vector{SymbolicT}
end
function SymbolicAffect(affect::Vector{Equation}; alg_eqs = Equation[],
        discrete_parameters = SymbolicT[], kwargs...)
    SymbolicAffect(affect, alg_eqs, discrete_parameters)
end
SymbolicAffect(affect; kwargs...) = make_affect(affect; kwargs...)
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
function AffectSystem(affect::Vector{Equation}; discrete_parameters = SymbolicT[],
        iv = nothing, alg_eqs::Vector{Equation} = Equation[], warn_no_algebraic = true, kwargs...)
    isempty(affect) && return nothing
    dvs = OrderedSet{SymbolicT}()
    params = OrderedSet{SymbolicT}()
    for eq in affect
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
function complete(cb::SymbolicContinuousCallback; kwargs...)
    SymbolicContinuousCallback(cb.conditions, make_affect(cb.affect; kwargs...),
        make_affect(cb.finalize; kwargs...), cb.rootfind, cb.reinitializealg, cb.zero_crossing_id)
end
make_affect(affect::SymbolicAffect; kwargs...) = AffectSystem(affect; kwargs...)
function make_affect(affect; kwargs...)
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
function create_symbolic_events(cont_events)
    cont_callbacks = to_cb_vector(cont_events; CB_TYPE = SymbolicContinuousCallback)
end
