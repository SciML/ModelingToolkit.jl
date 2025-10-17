const SYSTEM_COUNT = Threads.Atomic{UInt}(0)
struct GUIMetadata
end
Base.nameof(sys::AbstractSystem) = getfield(sys, :name)
function independent_variables(sys::AbstractSystem)
    if isdefined(sys, :iv) && getfield(sys, :iv) !== nothing
        return SymbolicT[getfield(sys, :iv)]
    end
end
function SymbolicIndexingInterface.is_parameter(sys::AbstractSystem, sym::Symbol)
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
                return ParameterIndex(
                    idx.portion, (idx.idx..., arguments(sym)[(begin + 1):end]...))
            end
        end
    end
end
function isscheduled(sys::AbstractSystem)
    if has_schedule(sys)
        get_schedule(sys) !== nothing
    end
end
struct Initial <: Symbolics.Operator end
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
                   :ps
                   :brownians
                   :observed
                   :systems
                   :initializesystem
                   :schedule
                   :is_initializesystem
                   :parameter_dependencies
                   :parent
                   :index_cache
                   :isscheduled
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
    return newmeta
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
Symbolics.rename(x::AbstractSystem, name) = @set x.name = name
apply_to_variables(f, ex) = _apply_to_variables(f, ex)
function _apply_to_variables(f::F, ex) where {F}
    if isvariable(ex)
        return f(ex)
    end
end
abstract type SymScope end
struct LocalScope <: SymScope end
function LocalScope(sym::Union{Num, SymbolicT, Symbolics.Arr{Num}})
    apply_to_variables(sym) do sym
        if iscall(sym) && operation(sym) === getindex
            maketerm(typeof(sym), operation(sym), [a1, args[2:end]...],
                metadata(sym))
        end
    end
end
struct ParentScope <: SymScope
end
function ParentScope(sym::Union{Num, SymbolicT, Symbolics.Arr{Num}})
    apply_to_variables(sym) do sym
        if iscall(sym) && operation(sym) === getindex
            a1 = setmetadata(args[1], SymScope,
                ParentScope(getmetadata(value(args[1]), SymScope, LocalScope())))
            maketerm(typeof(sym), operation(sym), [a1, args[2:end]...],
                ParentScope(getmetadata(value(sym), SymScope, LocalScope())))
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
const AllScopes = Union{LocalScope, ParentScope, GlobalScope}
renamespace(sys, tgt::Symbol) = Symbol(getname(sys), NAMESPACE_SEPARATOR_SYMBOL, tgt)
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
end
function parameters(sys::AbstractSystem; initial_parameters = false)
    ps = get_ps(sys)
    result = copy(ps)
    if !initial_parameters && !is_initializesystem(sys)
        filter!(result) do sym
            return !(isoperator(sym, Initial) ||
                     isoperator(arguments(sym)[1], Initial))
        end
    end
    return result
end
function full_parameters(sys::AbstractSystem)
    dep_ps = [eq.lhs for eq in get_parameter_dependencies(sys)]
    vcat(parameters(sys; initial_parameters = true), dep_ps)
end
function observed(sys::AbstractSystem)
    obs = get_observed(sys)
    systems = get_systems(sys)
    isempty(systems) && return obs
end
function observables(sys::AbstractSystem)
    return map(eq -> eq.lhs, observed(sys))
end
function equations(sys::AbstractSystem)
    eqs = get_eqs(sys)
end
function brownians(sys::AbstractSystem)
    bs = get_brownians(sys)
end
function constraints(sys::AbstractSystem)
    cs = get_constraints(sys)
end
function is_diff_equation(eq)
    (eq isa Equation) || (return false)
end
function is_alg_equation(eq)
    return (eq isa Equation) && !is_diff_equation(eq)
end
alg_equations(sys::AbstractSystem) = filter(is_alg_equation, equations(sys))
