const SYSTEM_COUNT = Threads.Atomic{UInt}(0)
struct GUIMetadata
end
Base.nameof(sys::AbstractSystem) = getfield(sys, :name)
function independent_variables(sys::AbstractSystem)
    if isdefined(sys, :iv) && getfield(sys, :iv) !== nothing
        return SymbolicT[getfield(sys, :iv)]
    end
end
function SymbolicIndexingInterface.is_variable(sys::AbstractSystem, sym)
    if sym isa Int
    end
end
function SymbolicIndexingInterface.variable_index(sys::AbstractSystem, sym::Symbol)
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
    end
end
function SymbolicIndexingInterface.is_parameter(sys::AbstractSystem, sym::Symbol)
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
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
                return ParameterIndex(
                    idx.portion, (idx.idx..., arguments(sym)[(begin + 1):end]...))
            end
        end
    end
    if idx === nothing && hasname(sym) && !(iscall(sym) && operation(sym) == getindex)
    end
end
function SymbolicIndexingInterface.independent_variable_symbols(sys::AbstractSystem)
end
function does_namespacing(sys::AbstractSystem)
    if isdefined(sys, :namespacing)
    end
end
function isscheduled(sys::AbstractSystem)
    if has_schedule(sys)
        get_schedule(sys) !== nothing
    end
end
struct Initial <: Symbolics.Operator end
function (f::Initial)(x)
    if symbolic_type(x) == NotSymbolic()
    end
    result = if SU.is_array_shape(sh)
    end
    if iw
    end
end
function add_initialization_parameters(sys::AbstractSystem; split = true)
    for x in unknowns(sys)
        if iscall(x) && operation(x) == getindex && split
        end
    end
    for eq in obs
    end
    for (i, v) in enumerate(initials)
    end
    for ivar in initials
        if symbolic_type(ivar) == ScalarSymbolic()
            for idx in SU.stable_eachindex(ivar)
            end
        end
    end
end
function isinitial(p)
    return iscall(p) && (operation(p) isa Initial ||
            operation(p) === getindex && isinitial(arguments(p)[1]))
end
function discover_globalscoped(sys::AbstractSystem)
end
function complete(
        sys::T; split = true, flatten = true, add_initial_parameters = true) where {T <: AbstractSystem}
    if flatten
        newsys = expand_connections(sys)
        if has_parent(newsys) && get_parent(sys) === nothing
            @set! newsys.parent = complete(sys; split = false, flatten = false)::T
        end
        if add_initial_parameters
        end
        cb_alg_eqs = Equation[alg_equations(sys); observed(sys)]
        if has_continuous_events(sys) && is_time_dependent(sys)
            for ev in get_continuous_events(sys)
                ev = complete(ev; iv = get_iv(sys)::SymbolicT, alg_eqs = cb_alg_eqs)
            end
            for ev in get_discrete_events(sys)
            end
        end
    end
    if split && has_index_cache(sys)
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
    if !isdefined(sys, :namespacing)
    end
    @set sys.namespacing = value
end
function unflatten_parameters!(buffer::Vector{SymbolicT}, params::Vector{SymbolicT}, all_ps::Set{SymbolicT})
    while i <= length(params)
    end
end
const SYS_PROPS = [:eqs
                   :noise_eqs
                   :iv
                   :unknowns
                   :ps
                   :brownians
                   :jumps
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
    for (k, v) in meta
    end
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
function Base.propertynames(sys::AbstractSystem; private = false)
    if private
        if has_parent(sys) && (parent = get_parent(sys); parent !== nothing)
        end
    end
    if has_ps(sys)
    end
    if has_eqs(sys)
        for eq in get_eqs(sys)
            if value(lhs) isa AnalysisPoint
            end
        end
    end
end
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
function renamespace(sys, x::SymbolicT)
    Moshi.Match.@match x begin
        BSImpl.Sym(; name) => let scope = getmetadata(x, SymScope, LocalScope())::AllScopes
            if scope isa LocalScope
            end
        end
        BSImpl.Term(; f, args, shape, type, metadata) => begin
            if f === getindex
            elseif f isa SymbolicT
                let scope = getmetadata(x, SymScope, LocalScope())::Union{LocalScope, ParentScope, GlobalScope}
                    if scope isa LocalScope
                        return rename(x, renamespace(getname(sys), getname(x)))::SymbolicT
                    end
                end
                for (i, arg) in enumerate(args)
                end
            end
        end
    end
end
function namespace_defaults(sys)
    Dict((isparameter(k) ? parameters(sys, k) : unknowns(sys, k)) => namespace_expr(v, sys)
    for (k, v) in pairs(defs))
end
function namespace_guesses(sys)
    for i in eachindex(O)
    end
    Moshi.Match.@match O begin
        BSImpl.Term(; f, args, metadata, type, shape) => begin
        end
    end
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
    if ps === SciMLBase.NullParameters()
    end
    result = copy(ps)
    if !initial_parameters && !is_initializesystem(sys)
        filter!(result) do sym
            return !(isoperator(sym, Initial) ||
                     isoperator(arguments(sym)[1], Initial))
        end
    end
    return result
end
function dependent_parameters(sys::AbstractSystem)
    if !iscomplete(sys)
        throw(ArgumentError("""
        """))
    end
end
function full_parameters(sys::AbstractSystem)
    dep_ps = [eq.lhs for eq in get_parameter_dependencies(sys)]
    vcat(parameters(sys; initial_parameters = true), dep_ps)
end
function observed(sys::AbstractSystem)
    obs = get_observed(sys)
    systems = get_systems(sys)
    isempty(systems) && return obs
    for subsys in systems
    end
end
function observables(sys::AbstractSystem)
    return map(eq -> eq.lhs, observed(sys))
end
function equations(sys::AbstractSystem)
    eqs = get_eqs(sys)
end
function equations_toplevel(sys::AbstractSystem)
    if has_parent(sys) && (parent = get_parent(sys)) !== nothing
    end
end
function substitute_observed(sys::AbstractSystem, expr; simplify = false)
    neweqs = map(equations(sys)) do eq
        if iscall(eq.lhs) && operation(eq.lhs) isa Union{Shift, Differential}
                   substitute_and_simplify(
                simplify)
            if !_iszero(eq.lhs)
            end
        end
    end
end
function jumps(sys::AbstractSystem)
    js = get_jumps(sys)
end
function brownians(sys::AbstractSystem)
    bs = get_brownians(sys)
    systems = get_systems(sys)
    if isempty(systems)
        return bs
    end
end
function cost(sys::AbstractSystem)
    systems = get_systems(sys)
    for subsys in systems
    end
end
function namespace_inequality(ineq::Inequality, sys, n = nameof(sys))
    Inequality(_lhs,
        ineq.relational_op)
end
function namespace_constraints(sys)
end
function constraints(sys::AbstractSystem)
    cs = get_constraints(sys)
    systems = get_systems(sys)
    isempty(systems) && return cs
    for subsys in systems
    end
    if isempty(systems)
        for sys in systems
            for eq in pre
            end
        end
    end
end
function process_parameter_equations(sys::AbstractSystem)
    if !isempty(get_systems(sys))
    end
    varsbuf = Set{SymbolicT}()
    pareq_idxs = Int[]
    eqs = equations(sys)
    for (i, eq) in enumerate(eqs)
        SU.search_variables!(varsbuf, eq; is_atomic = OperatorIsAtomic{Union{Differential, Initial, Pre}}())
        if let sys = sys
            all(varsbuf) do sym
                operation(sym) === getindex && is_parameter(sys, arguments(sym)[1])
        end
        end
            if !(eq.lhs in varsbuf)
                throw(ArgumentError("""
                """))
            end
        end
    end
    pareqs = Equation[get_parameter_dependencies(sys); eqs[pareq_idxs]]
    explicitpars = SymbolicT[]
    for eq in pareqs
    end
    @set! sys.ps = setdiff(get_ps(sys), explicitpars)
end
function is_diff_equation(eq)
    (eq isa Equation) || (return false)
end
function is_alg_equation(eq)
    return (eq isa Equation) && !is_diff_equation(eq)
end
alg_equations(sys::AbstractSystem) = filter(is_alg_equation, equations(sys))
