const SYSTEM_COUNT = Threads.Atomic{UInt}(0)
get_component_type(x::AbstractSystem) = get_gui_metadata(x).type
struct GUIMetadata
    type::GlobalRef
    layout::Any
end
GUIMetadata(type) = GUIMetadata(type, nothing)
const MTKPARAMETERS_ARG = SSym(:___mtkparameters___; type = Vector{Vector{Any}}, shape = SymbolicUtils.Unknown(1))
""""""
Base.nameof(sys::AbstractSystem) = getfield(sys, :name)
""""""
description(sys::AbstractSystem) = has_description(sys) ? get_description(sys) : ""
""""""
function independent_variables(sys::AbstractSystem)
    if isdefined(sys, :iv) && getfield(sys, :iv) !== nothing
        return SymbolicT[getfield(sys, :iv)]
    elseif isdefined(sys, :ivs)
        return unwrap.(getfield(sys, :ivs))::Vector{SymbolicT}
    else
        return SymbolicT[]
    end
end
function SymbolicIndexingInterface.is_variable(sys::AbstractSystem, sym)
    sym = unwrap(sym)
    if sym isa Int
        return sym in 1:length(variable_symbols(sys))
    end
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        return is_variable(ic, sym) ||
               iscall(sym) && operation(sym) === getindex &&
               is_variable(ic, first(arguments(sym)))
    end
    return any(isequal(sym), variable_symbols(sys)) ||
           hasname(sym) && is_variable(sys, getname(sym))
end
function SymbolicIndexingInterface.is_variable(sys::AbstractSystem, sym::Symbol)
    sym = unwrap(sym)
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        return is_variable(ic, sym)
    end
    return any(isequal(sym), getname.(variable_symbols(sys))) ||
           count(NAMESPACE_SEPARATOR, string(sym)) == 1 &&
           count(isequal(sym),
        Symbol.(nameof(sys), NAMESPACE_SEPARATOR_SYMBOL, getname.(variable_symbols(sys)))) ==
           1
end
function SymbolicIndexingInterface.variable_index(sys::AbstractSystem, sym)
    sym = unwrap(sym)
    if sym isa Int
        return sym
    end
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        return if (idx = variable_index(ic, sym)) !== nothing
            idx
        elseif iscall(sym) && operation(sym) === getindex &&
               (idx = variable_index(ic, first(arguments(sym)))) !== nothing
            idx[arguments(sym)[(begin + 1):end]...]
        else
            nothing
        end
    end
    idx = findfirst(isequal(sym), variable_symbols(sys))
    if idx === nothing && hasname(sym)
        idx = variable_index(sys, getname(sym))
    end
    return idx
end
function SymbolicIndexingInterface.variable_index(sys::AbstractSystem, sym::Symbol)
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        return variable_index(ic, sym)
    end
    idx = findfirst(isequal(sym), getname.(variable_symbols(sys)))
    if idx !== nothing
        return idx
    elseif count(NAMESPACE_SEPARATOR, string(sym)) == 1
        return findfirst(isequal(sym),
            Symbol.(
                nameof(sys), NAMESPACE_SEPARATOR_SYMBOL, getname.(variable_symbols(sys))))
    end
    return nothing
end
function SymbolicIndexingInterface.variable_symbols(sys::AbstractSystem)
    return unknowns(sys)
end
function SymbolicIndexingInterface.is_parameter(sys::AbstractSystem, sym::Union{Num, Symbolics.Arr, Symbolics.CallAndWrap})
    is_parameter(sys, unwrap(sym))
end
function SymbolicIndexingInterface.is_parameter(sys::AbstractSystem, sym::Int)
    !is_split(sys) && sym in 1:length(parameter_symbols(sys))
end
function SymbolicIndexingInterface.is_parameter(sys::AbstractSystem, sym::SymbolicT)
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        return sym isa ParameterIndex || is_parameter(ic, sym) ||
               iscall(sym) && operation(sym) === getindex &&
               is_parameter(ic, first(arguments(sym)))
    end
    return any(isequal(sym), parameter_symbols(sys)) ||
           hasname(sym) && !(iscall(sym) && operation(sym) == getindex) &&
           is_parameter(sys, getname(sym))
end
function SymbolicIndexingInterface.is_parameter(sys::AbstractSystem, sym::Symbol)
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        return is_parameter(ic, sym)
    end
    named_parameters = Symbol[getname(x)
                        for x in parameter_symbols(sys)
                        if hasname(x) && !(iscall(x) && operation(x) == getindex)]
    return any(isequal(sym), named_parameters) ||
           count(NAMESPACE_SEPARATOR, string(sym)) == 1 &&
           count(isequal(sym),
        Symbol.(nameof(sys), NAMESPACE_SEPARATOR_SYMBOL, named_parameters)) == 1
end
SymbolicIndexingInterface.is_parameter(sys::AbstractSystem, sym) = false
function SymbolicIndexingInterface.parameter_index(sys::AbstractSystem, sym)
    sym = unwrap(sym)
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        return if sym isa ParameterIndex
            sym
        elseif (idx = parameter_index(ic, sym)) !== nothing
            idx
        elseif iscall(sym) && operation(sym) === getindex &&
               (idx = parameter_index(ic, first(arguments(sym)))) !== nothing
            if idx.portion isa SciMLStructures.Tunable
                return ParameterIndex(
                    idx.portion, idx.idx[arguments(sym)[(begin + 1):end]...])
            else
                return ParameterIndex(
                    idx.portion, (idx.idx..., arguments(sym)[(begin + 1):end]...))
            end
        else
            nothing
        end
    end
    if sym isa Int
        return sym
    end
    idx = findfirst(isequal(sym), parameter_symbols(sys))
    if idx === nothing && hasname(sym) && !(iscall(sym) && operation(sym) == getindex)
        idx = parameter_index(sys, getname(sym))
    end
    return idx
end
function SymbolicIndexingInterface.parameter_index(sys::AbstractSystem, sym::Symbol)
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        idx = parameter_index(ic, sym)
        if idx === nothing ||
           idx.portion isa SciMLStructures.Discrete && idx.idx[2] == idx.idx[3] == 0
            return nothing
        else
            return idx
        end
    end
    pnames = [getname(x)
              for x in parameter_symbols(sys)
              if hasname(x) && !(iscall(x) && operation(x) == getindex)]
    idx = findfirst(isequal(sym), pnames)
    if idx !== nothing
        return idx
    elseif count(NAMESPACE_SEPARATOR, string(sym)) == 1
        return findfirst(isequal(sym),
            Symbol.(
                nameof(sys), NAMESPACE_SEPARATOR_SYMBOL, pnames))
    end
    return nothing
end
function SymbolicIndexingInterface.is_timeseries_parameter(sys::AbstractSystem, sym)
    is_time_dependent(sys) || return false
    has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing || return false
    is_timeseries_parameter(ic, sym)
end
function SymbolicIndexingInterface.timeseries_parameter_index(sys::AbstractSystem, sym)
    is_time_dependent(sys) || return nothing
    has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing || return nothing
    timeseries_parameter_index(ic, sym)
end
function SymbolicIndexingInterface.parameter_observed(sys::AbstractSystem, sym)
    return build_explicit_observed_function(sys, sym; param_only = true)
end
function SymbolicIndexingInterface.parameter_symbols(sys::AbstractSystem)
    return parameters(sys; initial_parameters = true)
end
function SymbolicIndexingInterface.is_independent_variable(sys::AbstractSystem, sym)
    return any(isequal(sym), independent_variable_symbols(sys))
end
function SymbolicIndexingInterface.is_independent_variable(sys::AbstractSystem, sym::Symbol)
    return any(isequal(sym), getname.(independent_variables(sys)))
end
function SymbolicIndexingInterface.independent_variable_symbols(sys::AbstractSystem)
    return independent_variables(sys)
end
function SymbolicIndexingInterface.is_observed(sys::AbstractSystem, sym)
    return !is_variable(sys, sym) && parameter_index(sys, sym) === nothing &&
           !is_independent_variable(sys, sym) && symbolic_type(sym) != NotSymbolic()
end
SymbolicIndexingInterface.supports_tuple_observed(::AbstractSystem) = true
SymbolicIndexingInterface.is_markovian(sys::AbstractSystem) = !is_dde(sys)
SymbolicIndexingInterface.constant_structure(::AbstractSystem) = true
""""""
iscomplete(sys::AbstractSystem) = isdefined(sys, :complete) && getfield(sys, :complete)
""""""
function does_namespacing(sys::AbstractSystem)
    if isdefined(sys, :namespacing)
        getfield(sys, :namespacing)
    else
        !iscomplete(sys)
    end
end
""""""
function isscheduled(sys::AbstractSystem)
    if has_schedule(sys)
        get_schedule(sys) !== nothing
    elseif has_isscheduled(sys)
        get_isscheduled(sys)
    else
        false
    end
end
""""""
struct Initial <: Symbolics.Operator end
is_timevarying_operator(::Type{Initial}) = false
Initial(x) = Initial()(x)
SymbolicUtils.promote_symtype(::Initial, ::Type{T}) where {T} = T
SymbolicUtils.promote_shape(::Initial, @nospecialize(x::SU.ShapeT)) = x
SymbolicUtils.isbinop(::Initial) = false
Base.nameof(::Initial) = :Initial
Base.show(io::IO, x::Initial) = print(io, "Initial")
function (f::Initial)(x)
    iw = Symbolics.iswrapped(x)
    x = unwrap(x)
    if symbolic_type(x) == NotSymbolic()
        return x
    end
    if iscall(x) && operation(x) isa Union{Differential, Shift}
        x = default_toterm(x)
    end
    iscall(x) && operation(x) isa Initial && return x
    sh = SU.shape(x)
    result = if SU.is_array_shape(sh)
        term(f, x; type = symtype(x), shape = sh)
    elseif iscall(x) && operation(x) === getindex
        arr = arguments(x)[1]
        f(arr)[arguments(x)[2:end]...]
    else
        term(f, x; type = symtype(x), shape = sh)
    end
    result = toparam(result)
    if iw
        result = wrap(result)
    end
    return result
end
supports_initialization(sys::AbstractSystem) = true
function add_initialization_parameters(sys::AbstractSystem; split = true)
    @assert !has_systems(sys) || isempty(get_systems(sys))
    supports_initialization(sys) || return sys
    is_initializesystem(sys) && return sys
    all_initialvars = Set{SymbolicT}()
    eqs = equations(sys)
    obs, eqs = unhack_observed(observed(sys), eqs)
    for x in unknowns(sys)
        if iscall(x) && operation(x) == getindex && split
            push!(all_initialvars, arguments(x)[1])
        else
            push!(all_initialvars, x)
        end
    end
    for eq in obs
        x = eq.lhs
        if iscall(x) && operation(x) == getindex && split
            push!(all_initialvars, arguments(x)[1])
        else
            push!(all_initialvars, x)
        end
    end
    if is_time_dependent(sys) && !is_discrete_system(sys)
        D = Differential(get_iv(sys)::SymbolicT)
        for v in collect(all_initialvars)
            iscall(v) && push!(all_initialvars, D(v))
        end
    end
    for eq in get_parameter_dependencies(sys)
        is_variable_floatingpoint(eq.lhs) || continue
        push!(all_initialvars, eq.lhs)
    end
    initials = collect(all_initialvars)
    for (i, v) in enumerate(initials)
        initials[i] = Initial()(v)
    end
    @set! sys.ps = unique!([filter(!isinitial, get_ps(sys)); initials])
    defs = copy(get_defaults(sys))
    for ivar in initials
        if symbolic_type(ivar) == ScalarSymbolic()
            defs[ivar] = false
        else
            defs[ivar] = collect(ivar)
            for idx in SU.stable_eachindex(ivar)
                scal_ivar = ivar[idx]
                defs[scal_ivar] = false
            end
        end
    end
    @set! sys.defaults = defs
    return sys
end
""""""
function isinitial(p)
    p = unwrap(p)
    return iscall(p) && (operation(p) isa Initial ||
            operation(p) === getindex && isinitial(arguments(p)[1]))
end
""""""
function discover_globalscoped(sys::AbstractSystem)
    newunknowns = OrderedSet{SymbolicT}()
    newparams = OrderedSet{SymbolicT}()
    iv::Union{SymbolicT, Nothing} = has_iv(sys) ? get_iv(sys) : nothing
    collect_scoped_vars!(newunknowns, newparams, sys, iv; depth = -1)
    setdiff!(newunknowns, observables(sys))
    @set! sys.ps = unique!(vcat(get_ps(sys), collect(newparams)))
    @set! sys.unknowns = unique!(vcat(get_unknowns(sys), collect(newunknowns)))
    return sys
end
""""""
function complete(
        sys::T; split = true, flatten = true, add_initial_parameters = true) where {T <: AbstractSystem}
    sys = discover_globalscoped(sys)
    if flatten
        newsys = expand_connections(sys)
        newsys = ModelingToolkit.flatten(newsys)
        if has_parent(newsys) && get_parent(sys) === nothing
            @set! newsys.parent = complete(sys; split = false, flatten = false)::T
        end
        sys = newsys
        sys = process_parameter_equations(sys)::T
        if add_initial_parameters
            sys = add_initialization_parameters(sys; split)::T
        end
        cb_alg_eqs = Equation[alg_equations(sys); observed(sys)]
        if has_continuous_events(sys) && is_time_dependent(sys)
            cevts = SymbolicContinuousCallback[]
            for ev in get_continuous_events(sys)
                ev = complete(ev; iv = get_iv(sys)::SymbolicT, alg_eqs = cb_alg_eqs)
                push!(cevts, ev)
            end
            @set! sys.continuous_events = cevts
        end
        if has_discrete_events(sys) && is_time_dependent(sys)
            devts = SymbolicDiscreteCallback[]
            for ev in get_discrete_events(sys)
                ev = complete(ev; iv = get_iv(sys)::SymbolicT, alg_eqs = cb_alg_eqs)
                push!(devts, ev)
            end
            @set! sys.discrete_events = devts
        end
    end
    if split && has_index_cache(sys)
        @set! sys.index_cache = IndexCache(sys)
        all_ps = parameters(sys; initial_parameters = true)
        all_ps_set = Set{SymbolicT}(all_ps)
        input_vars = inputs(sys)
        if !isempty(all_ps)
            ps_split = Vector{Vector{SymbolicT}}(reorder_parameters(sys, all_ps))
            ordered_ps = SymbolicT[]
            offset = 0
            if !isempty(get_index_cache(sys).tunable_idx)
                unflatten_parameters!(ordered_ps, ps_split[offset + 1], all_ps_set)
                offset += 1
            end
            if !isempty(get_index_cache(sys).initials_idx)
                unflatten_parameters!(ordered_ps, ps_split[offset + 1], all_ps_set)
                offset += 1
            end
            for i in (offset+1):length(ps_split)
                append!(ordered_ps, ps_split[i]::Vector{SymbolicT})
            end
            if isscheduled(sys)
                last_idx = 0
                for p in input_vars
                    idx = findfirst(isequal(p), ordered_ps)::Int
                    @assert last_idx < idx
                    last_idx = idx
                end
            end
            @set! sys.ps = ordered_ps
        end
    elseif has_index_cache(sys)
        @set! sys.index_cache = nothing
    end
    if has_initializesystem(sys)
        isys = get_initializesystem(sys)
        if isys isa T
            @set! sys.initializesystem = complete(isys::T; split)
        end
    end
    sys = toggle_namespacing(sys, false; safe = true)
    isdefined(sys, :complete) ? (@set! sys.complete = true) : sys
end
""""""
function toggle_namespacing(sys::AbstractSystem, value::Bool; safe = false)
    if !isdefined(sys, :namespacing)
        safe && return sys
        throw(ArgumentError("The system must define the `namespacing` flag to toggle namespacing"))
    end
    @set sys.namespacing = value
end
""""""
function unflatten_parameters!(buffer::Vector{SymbolicT}, params::Vector{SymbolicT}, all_ps::Set{SymbolicT})
    i = 1
    while i <= length(params)
        sym = params[i]
        if !iscall(sym) || operation(sym) !== getindex || sym in all_ps
            push!(buffer, sym)
            i += 1
            continue
        end
        arrsym = first(arguments(sym))
        for j in (i+1):(i+length(sym)-1)
            p = params[j]
            if !(iscall(p) && operation(p) === getindex && isequal(arguments(p)[1], arrsym))
                error("This should not be possible. Please open an issue in ModelingToolkit.jl with an MWE and stacktrace.")
            end
        end
        push!(buffer, arrsym)
        i += length(arrsym)
    end
end
const SYS_PROPS = [:eqs
                   :tag
                   :noise_eqs
                   :iv
                   :unknowns
                   :ps
                   :tspan
                   :brownians
                   :jumps
                   :name
                   :description
                   :var_to_name
                   :defaults
                   :guesses
                   :observed
                   :systems
                   :constraints
                   :bcs
                   :domain
                   :ivs
                   :dvs
                   :connector_type
                   :preface
                   :initializesystem
                   :initialization_eqs
                   :schedule
                   :tearing_state
                   :metadata
                   :gui_metadata
                   :is_initializesystem
                   :is_discrete
                   :parameter_dependencies
                   :assertions
                   :ignored_connections
                   :parent
                   :is_dde
                   :tstops
                   :inputs
                   :outputs
                   :index_cache
                   :isscheduled
                   :costs
                   :consolidate]
for prop in SYS_PROPS
    fname_get = Symbol(:get_, prop)
    fname_has = Symbol(:has_, prop)
    @eval begin
        """"""
$fname_get(sys::AbstractSystem) = getfield(sys, $(QuoteNode(prop)))
        """"""
$fname_has(sys::AbstractSystem) = isdefined(sys, $(QuoteNode(prop)))
    end
end
has_equations(::AbstractSystem) = true
""""""
function invalidate_cache!(sys::AbstractSystem)
    has_metadata(sys) || return sys
    empty!(getmetadata(sys, MutableCacheKey, nothing))
    return sys
end
function refreshed_metadata(meta::Base.ImmutableDict)
    newmeta = MetadataT()
    for (k, v) in meta
        if k === MutableCacheKey
            v = MutableCacheT()
        end
        newmeta = Base.ImmutableDict(newmeta, k => v)
    end
    if !haskey(newmeta, MutableCacheKey)
        newmeta = Base.ImmutableDict(newmeta, MutableCacheKey => MutableCacheT())
    end
    return newmeta
end
function Setfield.get(obj::AbstractSystem, ::Setfield.PropertyLens{field}) where {field}
    getfield(obj, field)
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
            Expr(:meta, :inline),
            Expr(:call, :(constructorof($obj)), args..., kwarg))
    else
        error("This should never happen. Trying to set $(typeof(obj)) with $patch.")
    end
end
Symbolics.rename(x::AbstractSystem, name) = @set x.name = name
function Base.propertynames(sys::AbstractSystem; private = false)
    if private
        return fieldnames(typeof(sys))
    else
        if has_parent(sys) && (parent = get_parent(sys); parent !== nothing)
            return propertynames(parent; private)
        end
        names = Symbol[]
        for s in get_systems(sys)
            push!(names, getname(s))
        end
        has_unknowns(sys) && for s in get_unknowns(sys)
            push!(names, getname(s))
        end
        has_ps(sys) && for s in get_ps(sys)
            hasname(s) || continue
            push!(names, getname(s))
        end
        has_observed(sys) && for s in get_observed(sys)
            push!(names, getname(s.lhs))
        end
        has_iv(sys) && push!(names, getname(get_iv(sys)))
        return names
    end
end
""""""
function Base.getproperty(
        sys::AbstractSystem, name::Symbol; namespace = does_namespacing(sys))
    if has_parent(sys) && (parent = get_parent(sys); parent !== nothing)
        return getproperty(parent, name; namespace)
    end
    wrap(getvar(sys, name; namespace = namespace))
end
function getvar(sys::AbstractSystem, name::Symbol; namespace = does_namespacing(sys))
    systems = get_systems(sys)
    if !isempty(systems)
        i = findfirst(x -> nameof(x) == name, systems)
        if i !== nothing
            return namespace ? renamespace(sys, systems[i]) : systems[i]
        end
    end
    if has_var_to_name(sys)
        avs = get_var_to_name(sys)
        v = get(avs, name, nothing)
        v === nothing || return namespace ? renamespace(sys, v) : v
    end
    sts = get_unknowns(sys)
    i = findfirst(x -> getname(x) == name, sts)
    if i !== nothing
        return namespace ? renamespace(sys, sts[i]) : sts[i]
    end
    if has_ps(sys)
        ps = get_ps(sys)
        i = findfirst(x -> getname(x) == name, ps)
        if i !== nothing
            return namespace ? renamespace(sys, ps[i]) : ps[i]
        end
    end
    if has_observed(sys)
        obs = get_observed(sys)
        i = findfirst(x -> getname(x.lhs) == name, obs)
        if i !== nothing
            return namespace ? renamespace(sys, obs[i].lhs) : obs[i].lhs
        end
    end
    if has_iv(sys)
        iv = get_iv(sys)
        if getname(iv) == name
            return iv
        end
    end
    if has_eqs(sys)
        for eq in get_eqs(sys)
            eq isa Equation || continue
            lhs = value(eq.lhs)
            rhs = value(eq.rhs)
            if value(lhs) isa AnalysisPoint
                rhs = rhs::AnalysisPoint
                nameof(rhs) == name || continue
                return namespace ? renamespace(sys, rhs) : rhs
            end
        end
    end
    throw(ArgumentError("System $(nameof(sys)): variable $name does not exist"))
end
""""""
apply_to_variables(f, ex) = _apply_to_variables(f, ex)
apply_to_variables(f, ex::Num) = wrap(_apply_to_variables(f, unwrap(ex)))
apply_to_variables(f, ex::Symbolics.Arr) = wrap(_apply_to_variables(f, unwrap(ex)))
function _apply_to_variables(f::F, ex) where {F}
    if isvariable(ex)
        return f(ex)
    end
    iscall(ex) || return ex
    maketerm(typeof(ex), _apply_to_variables(f, operation(ex)),
        map(Base.Fix1(_apply_to_variables, f), arguments(ex)),
        metadata(ex))
end
""""""
abstract type SymScope end
""""""
struct LocalScope <: SymScope end
""""""
function LocalScope(sym::Union{Num, SymbolicT, Symbolics.Arr{Num}})
    apply_to_variables(sym) do sym
        if iscall(sym) && operation(sym) === getindex
            args = arguments(sym)
            a1 = setmetadata(args[1], SymScope, LocalScope())
            maketerm(typeof(sym), operation(sym), [a1, args[2:end]...],
                metadata(sym))
        else
            setmetadata(sym, SymScope, LocalScope())
        end
    end
end
""""""
struct ParentScope <: SymScope
    parent::SymScope
end
""""""
function ParentScope(sym::Union{Num, SymbolicT, Symbolics.Arr{Num}})
    apply_to_variables(sym) do sym
        if iscall(sym) && operation(sym) === getindex
            args = arguments(sym)
            a1 = setmetadata(args[1], SymScope,
                ParentScope(getmetadata(value(args[1]), SymScope, LocalScope())))
            maketerm(typeof(sym), operation(sym), [a1, args[2:end]...],
                metadata(sym))
        else
            setmetadata(sym, SymScope,
                ParentScope(getmetadata(value(sym), SymScope, LocalScope())))
        end
    end
end
""""""
struct GlobalScope <: SymScope end
""""""
function GlobalScope(sym::Union{Num, SymbolicT, Symbolics.Arr{Num}})
    apply_to_variables(sym) do sym
        if iscall(sym) && operation(sym) == getindex
            args = arguments(sym)
            a1 = setmetadata(args[1], SymScope, GlobalScope())
            maketerm(typeof(sym), operation(sym), [a1, args[2:end]...],
                metadata(sym))
        else
            setmetadata(sym, SymScope, GlobalScope())
        end
    end
end
const AllScopes = Union{LocalScope, ParentScope, GlobalScope}
renamespace(sys, eq::Equation) = namespace_equation(eq, sys)
renamespace(names::AbstractVector, x) = foldr(renamespace, names, init = x)
renamespace(sys, tgt::AbstractSystem) = rename(tgt, renamespace(sys, nameof(tgt)))
renamespace(sys, tgt::Symbol) = Symbol(getname(sys), NAMESPACE_SEPARATOR_SYMBOL, tgt)
renamespace(sys, x::Num) = Num(renamespace(sys, unwrap(x)))
renamespace(sys, x::Arr{T, N}) where {T, N} = Arr{T, N}(renamespace(sys, unwrap(x)))
renamespace(sys, x::CallAndWrap{T}) where {T} = CallAndWrap{T}(renamespace(sys, unwrap(x)))
""""""
function renamespace(sys, x::SymbolicT)
    isequal(x, SU.idxs_for_arrayop(VartypeT)) && return x
    Moshi.Match.@match x begin
        BSImpl.Sym(; name) => let scope = getmetadata(x, SymScope, LocalScope())::AllScopes
            if scope isa LocalScope
                return rename(x, renamespace(getname(sys), name))::SymbolicT
            elseif scope isa ParentScope
                return setmetadata(x, SymScope, scope.parent)::SymbolicT
            elseif scope isa GlobalScope
                return x
            end
            error()
        end
        BSImpl.Term(; f, args, shape, type, metadata) => begin
            if f === getindex
                newargs = copy(parent(args))
                newargs[1] = renamespace(sys, args[1])
                return BSImpl.Term{VartypeT}(getindex, newargs; type, shape, metadata)
            elseif f isa SymbolicT
                let scope = getmetadata(x, SymScope, LocalScope())::Union{LocalScope, ParentScope, GlobalScope}
                    if scope isa LocalScope
                        return rename(x, renamespace(getname(sys), getname(x)))::SymbolicT
                    elseif scope isa ParentScope
                        return setmetadata(x, SymScope, scope.parent)::SymbolicT
                    elseif scope isa GlobalScope
                        return x
                    end
                    error()
                end
            elseif f isa Operator
                newargs = copy(parent(args))
                for (i, arg) in enumerate(args)
                    newargs[i] = renamespace(sys, arg)
                end
                return BSImpl.Term{VartypeT}(f, newargs; type, shape, metadata)
            end
            error()
        end
    end
end
namespace_variables(sys::AbstractSystem) = unknowns(sys, unknowns(sys))
namespace_parameters(sys::AbstractSystem) = parameters(sys, parameters(sys))
function namespace_defaults(sys)
    defs = defaults(sys)
    Dict((isparameter(k) ? parameters(sys, k) : unknowns(sys, k)) => namespace_expr(v, sys)
    for (k, v) in pairs(defs))
end
function namespace_guesses(sys)
    guess = guesses(sys)
    Dict(unknowns(sys, k) => namespace_expr(v, sys) for (k, v) in guess)
end
""""""
function namespace_equations(sys::AbstractSystem, ivs = independent_variables(sys))
    eqs = equations(sys)
    isempty(eqs) && return eqs
    if eqs === get_eqs(sys)
        eqs = copy(eqs)
    end
    for i in eachindex(eqs)
        eqs[i] = namespace_equation(eqs[i], sys; ivs)
    end
    return eqs
end
function namespace_initialization_equations(
        sys::AbstractSystem, ivs = independent_variables(sys))
    eqs = initialization_equations(sys)
    isempty(eqs) && return Equation[]
    map(eq -> namespace_equation(eq, sys; ivs), eqs)
end
function namespace_tstops(sys::AbstractSystem)
    tstops = symbolic_tstops(sys)
    isempty(tstops) && return tstops
    map(tstops) do val
        namespace_expr(val, sys)
    end
end
""""""
function namespace_equation(eq::Equation,
        sys,
        n = nameof(sys);
        ivs = independent_variables(sys))
    _lhs = namespace_expr(eq.lhs, sys, n; ivs)
    _rhs = namespace_expr(eq.rhs, sys, n; ivs)
    (_lhs ~ _rhs)::Equation
end
function namespace_jump(j::ConstantRateJump, sys)
    return ConstantRateJump(namespace_expr(j.rate, sys), namespace_expr(j.affect!, sys))
end
function namespace_jump(j::VariableRateJump, sys)
    return VariableRateJump(namespace_expr(j.rate, sys), namespace_expr(j.affect!, sys))
end
function namespace_jump(j::MassActionJump, sys)
    return MassActionJump(namespace_expr(j.scaled_rates, sys),
        [namespace_expr(k, sys) => namespace_expr(v, sys) for (k, v) in j.reactant_stoch],
        [namespace_expr(k, sys) => namespace_expr(v, sys) for (k, v) in j.net_stoch])
end
function namespace_jumps(sys::AbstractSystem)
    js = jumps(sys)
    isempty(js) && return js
    if js === get_jumps(sys)
        js = copy(js)
    end
    for i in eachindex(js)
        js[i] = namespace_jump(js[i], sys)
    end
    return js
end
function namespace_brownians(sys::AbstractSystem)
    bs = brownians(sys)
    if bs === get_brownians(sys)
        bs = copy(bs)
    end
    for i in eachindex(bs)
        bs[i] = renamespace(sys, bs[i])
    end
    return bs
end
function namespace_assignment(eq::Assignment, sys)
    _lhs = namespace_expr(eq.lhs, sys)
    _rhs = namespace_expr(eq.rhs, sys)
    Assignment(_lhs, _rhs)
end
function is_array_of_symbolics(x)
    symbolic_type(x) == ArraySymbolic() && return true
    symbolic_type(x) == ScalarSymbolic() && return false
    x isa AbstractArray &&
        any(y -> symbolic_type(y) != NotSymbolic() || is_array_of_symbolics(y), x)
end
function namespace_expr(O, sys::AbstractSystem, n::Symbol = nameof(sys); kw...)
    return O
end
function namespace_expr(O::Union{Num, Symbolics.Arr, Symbolics.CallAndWrap}, sys::AbstractSystem, n::Symbol = nameof(sys); kw...)
    typeof(O)(namespace_expr(unwrap(O), sys, n; kw...))
end
function namespace_expr(O::AbstractArray, sys::AbstractSystem, n::Symbol = nameof(sys); ivs = independent_variables(sys))
    is_array_of_symbolics(O) || return O
    O = copy(O)
    for i in eachindex(O)
        O[i] = namespace_expr(O[i], sys, n; ivs)
    end
    return O
end
function namespace_expr(O::SymbolicT, sys::AbstractSystem, n::Symbol = nameof(sys); ivs = independent_variables(sys))
    any(isequal(O), ivs) && return O
    isvar = isvariable(O)
    Moshi.Match.@match O begin
        BSImpl.Const(;) => return O
        BSImpl.Sym(;) => return isvar ? renamespace(n, O) : O
        BSImpl.Term(; f, args, metadata, type, shape) => begin
            newargs = copy(parent(args))
            for i in eachindex(args)
                newargs[i] = namespace_expr(newargs[i], sys, n; ivs)
            end
            if isvar
                rescoped = renamespace(n, O)
                f = Moshi.Data.variant_getfield(rescoped, BSImpl.Term{VartypeT}, :f)
                meta = Moshi.Data.variant_getfield(rescoped, BSImpl.Term{VartypeT}, :metadata)
            elseif f isa SymbolicT
                f = renamespace(n, f)
                meta = metadata
            else
                meta = metadata
            end
            return BSImpl.Term{VartypeT}(f, newargs; type, shape, metadata = meta)
        end
        BSImpl.AddMul(; coeff, dict, variant, type, shape, metadata) => begin
            newdict = copy(dict)
            empty!(newdict)
            for (k, v) in dict
                newdict[namespace_expr(k, sys, n; ivs)] = v
            end
            return BSImpl.AddMul{VartypeT}(coeff, newdict, variant; type, shape, metadata)
        end
        BSImpl.Div(; num, den, type, shape, metadata) => begin
            num = namespace_expr(num, sys, n; ivs)
            den = namespace_expr(den, sys, n; ivs)
            return BSImpl.Div{VartypeT}(num, den, false; type, shape, metadata)
        end
        BSImpl.ArrayOp(; output_idx, expr, term, ranges, reduce, type, shape, metadata) => begin
            if term isa SymbolicT
                term = namespace_expr(term, sys, n; ivs)
            end
            expr = namespace_expr(expr, sys, n; ivs)
            return BSImpl.ArrayOp{VartypeT}(output_idx, expr, reduce, term, ranges; type, shape, metadata)
        end
    end
end
""""""
function unknowns(sys::AbstractSystem)
    sts = get_unknowns(sys)
    systems = get_systems(sys)
    if isempty(systems)
        return sts
    end
    result = copy(sts)
    for subsys in systems
        append!(result, namespace_variables(subsys))
    end
    return result
end
""""""
function unknowns_toplevel(sys::AbstractSystem)
    if has_parent(sys) && (parent = get_parent(sys)) !== nothing
        return unknowns_toplevel(parent)
    end
    return get_unknowns(sys)
end
""""""
function parameters(sys::AbstractSystem; initial_parameters = false)
    ps = get_ps(sys)
    if ps === SciMLBase.NullParameters()
        return SymbolicT[]
    end
    if eltype(ps) <: Pair
        ps = Vector{SymbolicT}(unwrap.(first.(ps)))
    end
    systems = get_systems(sys)
    result = copy(ps)
    for subsys in systems
        append!(result, namespace_parameters(subsys))
    end
    if !initial_parameters && !is_initializesystem(sys)
        filter!(result) do sym
            return !(isoperator(sym, Initial) ||
                     iscall(sym) && operation(sym) === getindex &&
                     isoperator(arguments(sym)[1], Initial))
        end
    end
    return result
end
function dependent_parameters(sys::AbstractSystem)
    if !iscomplete(sys)
        throw(ArgumentError("""
        `dependent_parameters` requires that the system is marked as complete. Call
        `complete` or `mtkcompile` on the system.
        """))
    end
    return map(eq -> eq.lhs, parameter_dependencies(sys))
end
""""""
function parameters_toplevel(sys::AbstractSystem)
    if has_parent(sys) && (parent = get_parent(sys)) !== nothing
        return parameters_toplevel(parent)
    end
    return get_ps(sys)
end
""""""
function parameter_dependencies(sys::AbstractSystem)
    if !iscomplete(sys)
        throw(ArgumentError("""
        `parameter_dependencies` requires that the system is marked as complete. Call \
        `complete` or `mtkcompile` on the system.
        """))
    end
    if !has_parameter_dependencies(sys)
        return Equation[]
    end
    get_parameter_dependencies(sys)
end
""""""
function full_parameters(sys::AbstractSystem)
    dep_ps = [eq.lhs for eq in get_parameter_dependencies(sys)]
    vcat(parameters(sys; initial_parameters = true), dep_ps)
end
""""""
function assertions(sys::AbstractSystem)
    has_assertions(sys) || return Dict{BasicSymbolic, String}()
    asserts = get_assertions(sys)
    systems = get_systems(sys)
    namespaced_asserts = mapreduce(
        merge!, systems; init = Dict{BasicSymbolic, String}()) do subsys
        Dict{BasicSymbolic, String}(namespace_expr(k, subsys) => v
        for (k, v) in assertions(subsys))
    end
    return merge(asserts, namespaced_asserts)
end
""""""
struct IgnoredAnalysisPoint
    """
    The input variable/connector.
    """
    input::Union{BasicSymbolic, AbstractSystem}
    """
    The output variables/connectors.
    """
    outputs::Vector{Union{BasicSymbolic, AbstractSystem}}
end
""""""
function guesses(sys::AbstractSystem)
    guess = get_guesses(sys)
    systems = get_systems(sys)
    isempty(systems) && return guess
    for subsys in systems
        guess = merge(guess, namespace_guesses(subsys))
    end
    return guess
end
parameters(_) = []
""""""
function observed(sys::AbstractSystem)
    obs = get_observed(sys)
    systems = get_systems(sys)
    isempty(systems) && return obs
    obs = copy(obs)
    for subsys in systems
        _obs = observed(subsys)
        for eq in _obs
            push!(obs, namespace_equation(eq, subsys))
        end
    end
    return obs
end
""""""
function observables(sys::AbstractSystem)
    return map(eq -> eq.lhs, observed(sys))
end
""""""
function defaults(sys::AbstractSystem)
    systems = get_systems(sys)
    defs = get_defaults(sys)
    isempty(systems) ? defs : mapfoldr(namespace_defaults, merge, systems; init = defs)
end
function defaults_and_guesses(sys::AbstractSystem)
    merge(guesses(sys), defaults(sys))
end
unknowns(sys::Union{AbstractSystem, Nothing}, v) = namespace_expr(v, sys)
unknowns(sys::AbstractSystem, v::Symbolics.Arr) = namespace_expr(v, sys)
parameters(sys::AbstractSystem, v::Symbolics.Arr) = toparam(unknowns(sys, v))
parameters(sys::Union{AbstractSystem, Nothing}, v) = toparam(unknowns(sys, v))
for f in [:unknowns, :parameters]
    @eval function $f(sys::AbstractSystem, vs::AbstractArray)
        map(v -> $f(sys, v), vs)
    end
end
flatten(sys::AbstractSystem, args...) = sys
""""""
function equations(sys::AbstractSystem)
    eqs = get_eqs(sys)
    systems = get_systems(sys)
    isempty(systems) && return eqs
    eqs = copy(eqs)
    for subsys in systems
        append!(eqs, namespace_equations(subsys))
    end
    return eqs
end
""""""
function equations_toplevel(sys::AbstractSystem)
    iscomplete(sys) && error("Cannot apply `equations_toplevel` to complete systems.")
    if has_parent(sys) && (parent = get_parent(sys)) !== nothing
        return equations_toplevel(parent)
    end
    return get_eqs(sys)
end
""""""
function substitute_and_simplify(expr, dict::AbstractDict, simplify::Bool)
    expr = Symbolics.fixpoint_sub(expr, dict; operator = ModelingToolkit.Initial)
    simplify ? Symbolics.simplify(expr) : expr
end
""""""
function substitute_observed(sys::AbstractSystem, expr; simplify = false)
    empty_substitutions(sys) && return expr
    substitutions = get_substitutions(sys)
    return substitute_and_simplify(expr, substitutions, simplify)
end
""""""
function full_equations(sys::AbstractSystem; simplify = false)
    empty_substitutions(sys) && return equations(sys)
    subs = get_substitutions(sys)
    neweqs = map(equations(sys)) do eq
        if iscall(eq.lhs) && operation(eq.lhs) isa Union{Shift, Differential}
            return substitute_and_simplify(eq.lhs, subs, simplify) ~
                   substitute_and_simplify(
                eq.rhs, subs,
                simplify)
        else
            if !_iszero(eq.lhs)
                eq = 0 ~ eq.rhs - eq.lhs
            end
            rhs = substitute_and_simplify(eq.rhs, subs, simplify)
            return 0 ~ rhs
        end
        eq
    end
    return neweqs
end
""""""
function jumps(sys::AbstractSystem)
    js = get_jumps(sys)
    systems = get_systems(sys)
    isempty(systems) && return js
    js = copy(js)
    for subsys in systems
        append!(js, namespace_jumps(subsys))
    end
    return js
end
""""""
function brownians(sys::AbstractSystem)
    bs = get_brownians(sys)
    systems = get_systems(sys)
    if isempty(systems)
        return bs
    end
    bs = copy(bs)
    for subsys in systems
        append!(bs, namespace_brownians(subsys))
    end
    return bs
end
""""""
function cost(sys::AbstractSystem)
    cs = get_costs(sys)
    consolidate = get_consolidate(sys)
    systems = get_systems(sys)
    if isempty(systems)
        return consolidate(cs, Float64[])::SymbolicT
    end
    subcosts = SymbolicT[]
    for subsys in systems
        push!(subcosts, namespace_expr(cost(subsys), subsys))
    end
    return consolidate(cs, subcosts)::SymbolicT
end
namespace_constraint(eq::Equation, sys) = namespace_equation(eq, sys)
namespace_constraint(ineq::Inequality, sys) = namespace_inequality(ineq, sys)
function namespace_inequality(ineq::Inequality, sys, n = nameof(sys))
    _lhs = namespace_expr(ineq.lhs, sys, n)
    _rhs = namespace_expr(ineq.rhs, sys, n)
    Inequality(_lhs,
        _rhs,
        ineq.relational_op)
end
function namespace_constraints(sys)
    cstrs = constraints(sys)
    isempty(cstrs) && return cstrs
    if cstrs === get_constraints(sys)
        cstrs = copy(cstrs)
    end
    for i in eachindex(cstrs)
        cstrs[i] = namespace_constraint(cstrs[i], sys)
    end
    return cstrs
end
""""""
function constraints(sys::AbstractSystem)
    cs = get_constraints(sys)
    systems = get_systems(sys)
    isempty(systems) && return cs
    cs = copy(cs)
    for subsys in systems
        append!(cs, namespace_constraints(subsys))
    end
    return cs
end
""""""
function initialization_equations(sys::AbstractSystem)
    eqs = get_initialization_eqs(sys)
    systems = get_systems(sys)
    if isempty(systems)
        return eqs
    else
        eqs = Equation[eqs;
                       reduce(vcat,
                           namespace_initialization_equations.(get_systems(sys));
                           init = Equation[])]
        return eqs
    end
end
""""""
function symbolic_tstops(sys::AbstractSystem)
    tstops = get_tstops(sys)
    systems = get_systems(sys)
    isempty(systems) && return tstops
    tstops = [tstops; reduce(vcat, namespace_tstops.(get_systems(sys)); init = [])]
    return tstops
end

function preface(sys::AbstractSystem)
    has_preface(sys) || return nothing
    pre = get_preface(sys)
    systems = get_systems(sys)
    if isempty(systems)
        return pre
    else
        pres = pre === nothing ? [] : pre
        for sys in systems
            pre = get_preface(sys)
            pre === nothing && continue
            for eq in pre
                push!(pres, namespace_assignment(eq, sys))
            end
        end
        return isempty(pres) ? nothing : pres
    end
end
""""""
function process_parameter_equations(sys::AbstractSystem)
    if !isempty(get_systems(sys))
        throw(ArgumentError("Expected flattened system"))
    end
    varsbuf = Set{SymbolicT}()
    pareq_idxs = Int[]
    eqs = equations(sys)
    for (i, eq) in enumerate(eqs)
        empty!(varsbuf)
        SU.search_variables!(varsbuf, eq; is_atomic = OperatorIsAtomic{Union{Differential, Initial, Pre}}())
        isempty(varsbuf) && continue
        if let sys = sys
            all(varsbuf) do sym
            is_parameter(sys, sym) ||
                symbolic_type(sym) == ArraySymbolic() &&
                symbolic_has_known_size(sym) &&
                all(Base.Fix1(is_parameter, sys), collect(sym)) ||
                iscall(sym) &&
                operation(sym) === getindex && is_parameter(sys, arguments(sym)[1])
        end
        end
            if !(eq.lhs in varsbuf)
                throw(ArgumentError("""
                LHS of parameter dependency equation must be a single parameter. Found \
                $(eq.lhs).
                """))
            end
            push!(pareq_idxs, i)
        end
    end
    pareqs = Equation[get_parameter_dependencies(sys); eqs[pareq_idxs]]
    explicitpars = SymbolicT[]
    for eq in pareqs
        push!(explicitpars, eq.lhs)
    end
    pareqs = topsort_equations(pareqs, explicitpars)
    eqs = eqs[setdiff(eachindex(eqs), pareq_idxs)]
    @set! sys.eqs = eqs
    @set! sys.parameter_dependencies = pareqs
    @set! sys.ps = setdiff(get_ps(sys), explicitpars)
    return sys
end
function is_diff_equation(eq)
    (eq isa Equation) || (return false)
    isdefined(eq, :lhs) && recursive_hasoperator(Union{Differential, Shift}, eq.lhs) &&
        (return true)
    isdefined(eq, :rhs) && recursive_hasoperator(Union{Differential, Shift}, eq.rhs) &&
        (return true)
    return false
end
function is_alg_equation(eq)
    return (eq isa Equation) && !is_diff_equation(eq)
end
alg_equations(sys::AbstractSystem) = filter(is_alg_equation, equations(sys))
diff_equations(sys::AbstractSystem) = filter(is_diff_equation, equations(sys))
has_alg_equations(sys::AbstractSystem) = any(is_alg_equation, equations(sys))
has_diff_equations(sys::AbstractSystem) = any(is_diff_equation, equations(sys))
get_alg_eqs(sys::AbstractSystem) = filter(is_alg_equation, get_eqs(sys))
get_diff_eqs(sys::AbstractSystem) = filter(is_diff_equation, get_eqs(sys))
has_alg_eqs(sys::AbstractSystem) = any(is_alg_equation, get_eqs(sys))
has_diff_eqs(sys::AbstractSystem) = any(is_diff_equation, get_eqs(sys))
