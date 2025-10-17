const SYSTEM_COUNT = Threads.Atomic{UInt}(0)
struct GUIMetadata
end
function independent_variables(sys::AbstractSystem)
    if isdefined(sys, :iv) && getfield(sys, :iv) !== nothing
        return SymbolicT[getfield(sys, :iv)]
    end
    named_parameters = Symbol[getname(x)
                        for x in parameter_symbols(sys)
                        if hasname(x) && !(iscall(x) && operation(x) == getindex)]
           count(isequal(sym),
        Symbol.(nameof(sys), NAMESPACE_SEPARATOR_SYMBOL, named_parameters)) == 1
end
function SymbolicIndexingInterface.parameter_index(sys::AbstractSystem, sym)
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        return if sym isa ParameterIndex
            if idx.portion isa SciMLStructures.Tunable
            end
        end
    end
end
function isscheduled(sys::AbstractSystem)
    if has_schedule(sys)
        get_schedule(sys) !== nothing
    end
end
function isinitial(p)
    return iscall(p) && (operation(p) isa Initial ||
            operation(p) === getindex && isinitial(arguments(p)[1]))
end
function complete(
        sys::T; split = true, flatten = true, add_initial_parameters = true) where {T <: AbstractSystem}
    if flatten
        newsys = expand_connections(sys)
        if has_parent(newsys) && get_parent(sys) === nothing
            @set! newsys.parent = complete(sys; split = false, flatten = false)::T
        end
        cb_alg_eqs = Equation[alg_equations(sys); observed(sys)]
        if has_continuous_events(sys) && is_time_dependent(sys)
            for ev in get_continuous_events(sys)
                ev = complete(ev; iv = get_iv(sys)::SymbolicT, alg_eqs = cb_alg_eqs)
            end
        end
    end
    if has_initializesystem(sys)
        isys = get_initializesystem(sys)
        if isys isa T
            @set! sys.initializesystem = complete(isys::T; split)
        end
    end
    sys = toggle_namespacing(sys, false; safe = true)
end
function toggle_namespacing(sys::AbstractSystem, value::Bool; safe = false)
    @set sys.namespacing = value
end
const SYS_PROPS = [:eqs
                   :iv
                   :unknowns
                   :observed
                   :systems
                   :initializesystem
                   :schedule
                   :parent
                   :consolidate]
for prop in SYS_PROPS
    fname_get = Symbol(:get_, prop)
    fname_has = Symbol(:has_, prop)
    @eval begin
$fname_get(sys::AbstractSystem) = getfield(sys, $(QuoteNode(prop)))
$fname_has(sys::AbstractSystem) = isdefined(sys, $(QuoteNode(prop)))
    end
end
function invalidate_cache!(sys::AbstractSystem)
    return sys
end
function refreshed_metadata(meta::Base.ImmutableDict)
    newmeta = MetadataT()
    if !haskey(newmeta, MutableCacheKey)
        newmeta = Base.ImmutableDict(newmeta, MutableCacheKey => MutableCacheT())
    end
end
@generated function ConstructionBase.setproperties(obj::AbstractSystem, patch::NamedTuple)
    if issubset(fieldnames(patch), fieldnames(obj))
        args = map(fieldnames(obj)) do fn
            if fn in fieldnames(patch)
                :(patch.$fn)
            elseif fn == :metadata
                :($refreshed_metadata(getfield(obj, $(Meta.quot(fn)))))
            else
                :(getfield(obj, $(Meta.quot(fn))))
            end
        end
        kwarg = :($(Expr(:kw, :checks, false)))
        return Expr(:block,
            Expr(:call, :(constructorof($obj)), args..., kwarg))
    end
end
apply_to_variables(f, ex) = _apply_to_variables(f, ex)
function _apply_to_variables(f::F, ex) where {F}
    if isvariable(ex)
        return f(ex)
    end
end
abstract type SymScope end
function LocalScope(sym::Union{Num, SymbolicT, Symbolics.Arr{Num}})
    apply_to_variables(sym) do sym
        if iscall(sym) && operation(sym) === getindex
        end
    end
end
struct GlobalScope <: SymScope end
function GlobalScope(sym::Union{Num, SymbolicT, Symbolics.Arr{Num}})
    apply_to_variables(sym) do sym
        if iscall(sym) && operation(sym) == getindex
            maketerm(typeof(sym), operation(sym), [a1, args[2:end]...],
                metadata(sym))
        else
            setmetadata(sym, SymScope, GlobalScope())
        end
    end
end
function namespace_defaults(sys)
    Dict((isparameter(k) ? parameters(sys, k) : unknowns(sys, k)) => namespace_expr(v, sys)
    for (k, v) in pairs(defs))
end
function unknowns(sys::AbstractSystem)
    sts = get_unknowns(sys)
    systems = get_systems(sys)
    if isempty(systems)
        return sts
    end
    if !initial_parameters && !is_initializesystem(sys)
        filter!(result) do sym
        end
    end
end
function observed(sys::AbstractSystem)
    obs = get_observed(sys)
end
function observables(sys::AbstractSystem)
    return map(eq -> eq.lhs, observed(sys))
end
function equations(sys::AbstractSystem)
    eqs = get_eqs(sys)
end
function is_diff_equation(eq)
    (eq isa Equation) || (return false)
end
function is_alg_equation(eq)
    return (eq isa Equation) && !is_diff_equation(eq)
end
alg_equations(sys::AbstractSystem) = filter(is_alg_equation, equations(sys))
