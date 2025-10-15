const SYSTEM_COUNT = Threads.Atomic{UInt}(0)
get_component_type(x::AbstractSystem) = get_gui_metadata(x).type
struct GUIMetadata
    type::GlobalRef
    layout::Any
end
GUIMetadata(type) = GUIMetadata(type, nothing)
""""""
function generate_custom_function(sys::AbstractSystem, exprs, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true);
        expression = Val{true}, eval_expression = false, eval_module = @__MODULE__,
        cachesyms::Tuple = (), kwargs...)
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `mtkcompile` on the system.")
    end
    p = (reorder_parameters(sys, unwrap.(ps))..., cachesyms...)
    isscalar = !(exprs isa AbstractArray)
    fnexpr = if is_time_dependent(sys)
        build_function_wrapper(sys, exprs,
            dvs,
            p...,
            get_iv(sys);
            kwargs...,
            expression = Val{true}
        )
    else
        build_function_wrapper(sys, exprs,
            dvs,
            p...;
            kwargs...,
            expression = Val{true}
        )
    end
    if expression == Val{true}
        return fnexpr
    end
    if fnexpr isa Tuple
        return eval_or_rgf.(fnexpr; eval_expression, eval_module)
    else
        return eval_or_rgf(fnexpr; eval_expression, eval_module)
    end
end
function wrap_assignments(isscalar, assignments; let_block = false)
    function wrapper(expr)
        Func(expr.args, [], Let(assignments, expr.body, let_block))
    end
    if isscalar
        wrapper
    else
        wrapper, wrapper
    end
end
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
""""""
function has_observed_with_lhs(sys::AbstractSystem, sym)
    has_observed(sys) || return false
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        return haskey(ic.observed_syms_to_timeseries, sym)
    else
        return any(isequal(sym), observables(sys))
    end
end
""""""
function has_parameter_dependency_with_lhs(sys, sym)
    has_parameter_dependencies(sys) || return false
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        return haskey(ic.dependent_pars_to_timeseries, unwrap(sym))
    else
        return any(isequal(sym), [eq.lhs for eq in get_parameter_dependencies(sys)])
    end
end
function _all_ts_idxs!(ts_idxs, ::NotSymbolic, sys, sym)
    if is_variable(sys, sym) || is_independent_variable(sys, sym)
        push!(ts_idxs, ContinuousTimeseries())
    elseif is_timeseries_parameter(sys, sym)
        push!(ts_idxs, timeseries_parameter_index(sys, sym).timeseries_idx)
    end
end
for traitT in [
    ScalarSymbolic,
    ArraySymbolic
]
    @eval function _all_ts_idxs!(ts_idxs, ::$traitT, sys, sym)
        allsyms = Set{SymbolicT}()
        SU.search_variables!(allsyms, sym; is_atomic = OperatorIsAtomic{Symbolics.Operator}())
        for s in allsyms
            s = unwrap(s)
            if is_variable(sys, s) || is_independent_variable(sys, s)
                push!(ts_idxs, ContinuousTimeseries())
            elseif is_timeseries_parameter(sys, s)
                push!(ts_idxs, timeseries_parameter_index(sys, s).timeseries_idx)
            elseif is_time_dependent(sys) && iscall(s) && issym(operation(s)) &&
                   length(arguments(s)) == 1 && is_variable(sys, operation(s)(get_iv(sys)))
                push!(ts_idxs, ContinuousTimeseries())
            else
                if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
                    if (ts = get(ic.observed_syms_to_timeseries, s, nothing)) !== nothing
                        union!(ts_idxs, ts)
                    elseif (ts = get(ic.dependent_pars_to_timeseries, s, nothing)) !==
                           nothing
                        union!(ts_idxs, ts)
                    end
                else
                    if has_observed_with_lhs(sys, sym)
                        push!(ts_idxs, ContinuousTimeseries())
                    end
                end
            end
        end
    end
end
function _all_ts_idxs!(ts_idxs, ::ScalarSymbolic, sys, sym::Symbol)
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        return _all_ts_idxs!(ts_idxs, sys, ic.symbol_to_variable[sym])
    elseif is_variable(sys, sym) || is_independent_variable(sys, sym) ||
           any(isequal(sym), getname.(observables(sys)))
        push!(ts_idxs, ContinuousTimeseries())
    elseif is_timeseries_parameter(sys, sym)
        push!(ts_idxs, timeseries_parameter_index(sys, s).timeseries_idx)
    end
end
function _all_ts_idxs!(ts_idxs, ::NotSymbolic, sys, sym::AbstractArray)
    for s in sym
        _all_ts_idxs!(ts_idxs, sys, s)
    end
end
_all_ts_idxs!(ts_idxs, sys, sym) = _all_ts_idxs!(ts_idxs, symbolic_type(sym), sys, sym)
function SymbolicIndexingInterface.get_all_timeseries_indexes(sys::AbstractSystem, sym)
    if !is_time_dependent(sys)
        return Set()
    end
    ts_idxs = Set()
    _all_ts_idxs!(ts_idxs, sys, sym)
    return ts_idxs
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
function SymbolicIndexingInterface.observed(
        sys::AbstractSystem, sym; eval_expression = false, eval_module = @__MODULE__,
        checkbounds = true, cse = true)
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        if sym isa Symbol
            _sym = get(ic.symbol_to_variable, sym, nothing)
            if _sym === nothing
                throw(ArgumentError("Symbol $sym does not exist in the system"))
            end
            sym = _sym
        elseif (sym isa Tuple ||
                (sym isa AbstractArray && symbolic_type(sym) isa NotSymbolic)) &&
               any(x -> x isa Symbol, sym)
            sym = map(sym) do s
                if s isa Symbol
                    _s = get(ic.symbol_to_variable, s, nothing)
                    if _s === nothing
                        throw(ArgumentError("Symbol $s does not exist in the system"))
                    end
                    return _s
                end
                return unwrap(s)
            end
        end
    end
    return build_explicit_observed_function(
        sys, sym; eval_expression, eval_module, checkbounds, cse)
end
function SymbolicIndexingInterface.default_values(sys::AbstractSystem)
    return merge(
        Dict(eq.lhs => eq.rhs for eq in observed(sys)),
        defaults(sys)
    )
end
SymbolicIndexingInterface.is_markovian(sys::AbstractSystem) = !is_dde(sys)
SymbolicIndexingInterface.constant_structure(::AbstractSystem) = true
function SymbolicIndexingInterface.all_variable_symbols(sys::AbstractSystem)
    syms = variable_symbols(sys)
    obs = observables(sys)
    return isempty(obs) ? syms : vcat(syms, obs)
end
function SymbolicIndexingInterface.all_symbols(sys::AbstractSystem)
    syms = all_variable_symbols(sys)
    for other in (full_parameters(sys), independent_variable_symbols(sys))
        isempty(other) || (syms = vcat(syms, other))
    end
    return syms
end
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
function Base.setproperty!(sys::AbstractSystem, prop::Symbol, val)
    error("""
    `setproperty!` on systems is invalid. Systems are immutable data structures, and \
    modifications to fields should be made by constructing a new system. This can be done \
    easily using packages such as Setfield.jl.

    If you are looking for the old behavior of updating the default of a variable via \
    `setproperty!`, this should now be done by mutating `ModelingToolkit.get_defaults(sys)`.
    """)
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
""""""
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
function islinear(sys::AbstractSystem)
    rhs = [eq.rhs for eq in full_equations(sys)]
    all(islinear(r, unknowns(sys)) for r in rhs)
end
function isaffine(sys::AbstractSystem)
    rhs = [eq.rhs for eq in full_equations(sys)]
    all(isaffine(r, unknowns(sys)) for r in rhs)
end
struct ObservedFunctionCache{S}
    sys::S
    dict::Dict{Any, Any}
    steady_state::Bool
    eval_expression::Bool
    eval_module::Module
    checkbounds::Bool
    cse::Bool
end
function ObservedFunctionCache(
        sys; expression = Val{false}, steady_state = false, eval_expression = false,
        eval_module = @__MODULE__, checkbounds = true, cse = true)
    if expression == Val{true}
        :($ObservedFunctionCache($sys, Dict(), $steady_state, $eval_expression,
            $eval_module, $checkbounds, $cse))
    else
        ObservedFunctionCache(
            sys, Dict(), steady_state, eval_expression, eval_module, checkbounds, cse)
    end
end
function Base.deepcopy_internal(ofc::ObservedFunctionCache, stackdict::IdDict)
    sys = deepcopy(ofc.sys)
    dict = deepcopy(ofc.dict)
    steady_state = ofc.steady_state
    eval_expression = ofc.eval_expression
    eval_module = ofc.eval_module
    checkbounds = ofc.checkbounds
    cse = ofc.cse
    newofc = ObservedFunctionCache(
        sys, dict, steady_state, eval_expression, eval_module, checkbounds, cse)
    stackdict[ofc] = newofc
    return newofc
end
function (ofc::ObservedFunctionCache)(obsvar, args...)
    obs = get!(ofc.dict, value(obsvar)) do
        SymbolicIndexingInterface.observed(
            ofc.sys, obsvar; eval_expression = ofc.eval_expression,
            eval_module = ofc.eval_module, checkbounds = ofc.checkbounds, cse = ofc.cse)
    end
    if ofc.steady_state
        obs = let fn = obs
            fn1(u, p, t = Inf) = fn(u, p, t)
            fn1
        end
    end
    if args === ()
        return obs
    else
        return obs(args...)
    end
end
function push_vars!(stmt, name, typ, vars)
    isempty(vars) && return
    vars_expr = Expr(:macrocall, typ, nothing)
    for s in vars
        if iscall(s)
            f = nameof(operation(s))
            args = arguments(s)
            ex = :($f($(args...)))
        else
            ex = nameof(s)
        end
        push!(vars_expr.args, ex)
        meta_kvps = Expr[]
        if isinput(s)
            push!(meta_kvps, :(input = true))
        end
        if isoutput(s)
            push!(meta_kvps, :(output = true))
        end
        if !isempty(meta_kvps)
            push!(vars_expr.args, Expr(:vect, meta_kvps...))
        end
    end
    push!(stmt, :($name = $collect($vars_expr)))
    return
end
function round_trip_expr(t, var2name)
    name = get(var2name, t, nothing)
    name !== nothing && return name
    issym(t) && return nameof(t)
    iscall(t) || return t
    f = round_trip_expr(operation(t), var2name)
    args = map(Base.Fix2(round_trip_expr, var2name), arguments(t))
    return :($f($(args...)))
end
function round_trip_eq(eq::Equation, var2name)
    if eq.lhs isa Connection
        syss = get_systems(eq.rhs)
        call = Expr(:call, connect)
        for sys in syss
            strs = split(string(nameof(sys)), NAMESPACE_SEPARATOR)
            s = Symbol(strs[1])
            for st in strs[2:end]
                s = Expr(:., s, Meta.quot(Symbol(st)))
            end
            push!(call.args, s)
        end
        call
    else
        Expr(:call, (~), round_trip_expr(eq.lhs, var2name),
            round_trip_expr(eq.rhs, var2name))
    end
end
function push_eqs!(stmt, eqs, var2name)
    eqs_name = gensym(:eqs)
    eqs_expr = Expr(:vcat)
    eqs_blk = Expr(:(=), eqs_name, eqs_expr)
    for eq in eqs
        push!(eqs_expr.args, round_trip_eq(eq, var2name))
    end
    push!(stmt, eqs_blk)
    return eqs_name
end
function push_defaults!(stmt, defs, var2name; name = :defs)
    defs_name = gensym(name)
    defs_expr = Expr(:call, Dict)
    defs_blk = Expr(:(=), defs_name, defs_expr)
    for d in defs
        n = round_trip_expr(d.first, var2name)
        v = round_trip_expr(d.second, var2name)
        push!(defs_expr.args, :($(=>)($n, $v)))
    end
    push!(stmt, defs_blk)
    return defs_name
end
function toexpr(sys::AbstractSystem)
    sys = flatten(sys)
    expr = Expr(:block)
    stmt = expr.args
    name = Meta.quot(nameof(sys))
    ivs = independent_variables(sys)
    ivname = gensym(:iv)
    for iv in ivs
        ivname = gensym(:iv)
        push!(stmt, :($ivname = (@variables $(getname(iv)))[1]))
    end
    stsname = gensym(:sts)
    sts = unknowns(sys)
    push_vars!(stmt, stsname, Symbol("@variables"), sts)
    psname = gensym(:ps)
    ps = parameters(sys)
    push_vars!(stmt, psname, Symbol("@parameters"), ps)
    obs = observed(sys)
    obsvars = [o.lhs for o in obs]
    obsvarsname = gensym(:obs)
    push_vars!(stmt, obsvarsname, Symbol("@variables"), obsvars)
    var2name = Dict{Any, Symbol}()
    for v in Iterators.flatten((sts, ps, obsvars))
        var2name[v] = getname(v)
    end
    eqs_name = push_eqs!(stmt, full_equations(sys), var2name)
    filtered_defs = filter(
        kvp -> !(iscall(kvp[1]) && operation(kvp[1]) isa Initial), defaults(sys))
    filtered_guesses = filter(
        kvp -> !(iscall(kvp[1]) && operation(kvp[1]) isa Initial), guesses(sys))
    defs_name = push_defaults!(stmt, filtered_defs, var2name)
    guesses_name = push_defaults!(stmt, filtered_guesses, var2name; name = :guesses)
    obs_name = push_eqs!(stmt, obs, var2name)
    iv = get_iv(sys)
    if iv === nothing
        ivname = nothing
    else
        ivname = gensym(:iv)
        push!(stmt, :($ivname = (@variables $(getname(iv)))[1]))
    end
    push!(stmt,
        :($System($eqs_name, $ivname, $stsname, $psname; defaults = $defs_name,
            guesses = $guesses_name, observed = $obs_name,
            name = $name, checks = false)))
    expr = :(let
        $expr
    end)
    Base.remove_linenums!(expr)
end
Base.write(io::IO, sys::AbstractSystem) = write(io, readable_code(toexpr(sys)))
""""""
function n_expanded_connection_equations(sys::AbstractSystem)
    isconnector(sys) && return length(get_unknowns(sys))
    sys = remove_analysis_points(sys)
    sys, (csets, _) = generate_connection_set(sys)
    n_extras = 0
    for cset in csets
        rep = cset[1]
        if rep.type <: Union{InputVar, OutputVar, Equality}
            n_extras += length(cset) - 1
        elseif rep.type == Flow
            n_extras += 1
        elseif rep.type == Stream
            n_extras += count(x -> x.isouter, cset)
        end
    end
    return n_extras
end
Base.show(io::IO, sys::AbstractSystem; kws...) = show(io, MIME"text/plain"(), sys; kws...)
function Base.show(
        io::IO, mime::MIME"text/plain", sys::AbstractSystem; hint = true, bold = true)
    limit = get(io, :limit, false)
    rows = first(displaysize(io)) ÷ 5
    desc = description(sys)
    name = nameof(sys)
    printstyled(io, "Model ", name, ":"; bold)
    !isempty(desc) && print(io, " ", desc)
    subs = get_systems(sys)
    nsubs = length(subs)
    nrows = min(nsubs, limit ? rows : nsubs)
    nrows > 0 && printstyled(io, "\nSubsystems ($(nsubs)):"; bold)
    nrows > 0 && hint && print(io, " see hierarchy($name)")
    for i in 1:nrows
        sub = subs[i]
        local name = String(nameof(sub))
        print(io, "\n  ", name)
        desc = description(sub)
        if !isempty(desc)
            maxlen = displaysize(io)[2] - length(name) - 6
            if limit && length(desc) > maxlen
                desc = chop(desc, tail = length(desc) - maxlen) * "…"
            end
            print(io, ": ", desc)
        end
    end
    limited = nrows < nsubs
    limited && print(io, "\n  ⋮")
    eqs = equations(sys)
    if eqs isa AbstractArray && eltype(eqs) <: Equation
        neqs = count(eq -> !(eq.lhs isa Connection), eqs)
        next = n_expanded_connection_equations(sys)
        ntot = neqs + next
        ntot > 0 && printstyled(io, "\nEquations ($ntot):"; bold)
        neqs > 0 && print(io, "\n  $neqs standard", hint ? ": see equations($name)" : "")
        next > 0 && print(io, "\n  $next connecting",
            hint ? ": see equations(expand_connections($name))" : "")
    end
    for varfunc in [unknowns, parameters]
        vars = varfunc(sys)
        nvars = length(vars)
        nvars == 0 && continue
        header = titlecase(String(nameof(varfunc)))
        printstyled(io, "\n$header ($nvars):"; bold)
        hint && print(io, " see $(nameof(varfunc))($name)")
        nrows = min(nvars, limit ? rows : nvars)
        defs = has_defaults(sys) ? defaults(sys) : nothing
        for i in 1:nrows
            s = vars[i]
            print(io, "\n  ", s)
            if !isnothing(defs)
                val = get(defs, s, nothing)
                if !isnothing(val)
                    print(io, " [defaults to ")
                    show(
                        IOContext(io, :compact => true, :limit => true,
                            :displaysize => (1, displaysize(io)[2])),
                        val)
                    print(io, "]")
                end
                desc = getdescription(s)
            end
            if !isnothing(desc) && desc != ""
                print(io, ": ", desc)
            end
        end
        limited = nrows < nvars
        limited && printstyled(io, "\n  ⋮")
    end
    npdeps = has_parameter_dependencies(sys) ? length(get_parameter_dependencies(sys)) : 0
    npdeps > 0 && printstyled(io, "\nParameter dependencies ($npdeps):"; bold)
    npdeps > 0 && hint && print(io, " see parameter_dependencies($name)")
    nobs = has_observed(sys) ? length(observed(sys)) : 0
    nobs > 0 && printstyled(io, "\nObserved ($nobs):"; bold)
    nobs > 0 && hint && print(io, " see observed($name)")
    return nothing
end
function split_assign(expr)
    if !(expr isa Expr && expr.head === :(=) && expr.args[2].head === :call)
        throw(ArgumentError("expression should be of the form `sys = foo(a, b)`"))
    end
    name, call = expr.args
end
varname_fix!(s) = return
function varname_fix!(expr::Expr)
    for arg in expr.args
        MLStyle.@match arg begin
            ::Symbol => continue
            Expr(:kw, a...) || Expr(:kw, a) => varname_sanitization!(arg)
            Expr(:parameters, a...) => begin
                for _arg in arg.args
                    varname_sanitization!(_arg)
                end
            end
            _ => @debug "skipping variable sanitization of $arg"
        end
    end
end
varname_sanitization!(a) = return
function varname_sanitization!(expr::Expr)
    var_splits = split(string(expr.args[1]), ".")
    if length(var_splits) > 1
        expr.args[1] = Symbol(join(var_splits, "__"))
    end
end
function _named(name, call, runtime = false)
    has_kw = false
    call isa Expr || throw(Meta.ParseError("The rhs must be an Expr. Got $call."))
    if length(call.args) >= 2 && call.args[2] isa Expr
        if call.args[2].head === :kw
            call = Expr(call.head, call.args[1], Expr(:parameters, call.args[2:end]...))
            has_kw = true
        elseif call.args[2].head === :parameters
            has_kw = true
        end
    end
    varname_fix!(call)
    if !has_kw
        param = Expr(:parameters)
        if length(call.args) == 1
            push!(call.args, param)
        else
            insert!(call.args, 2, param)
        end
    end
    is_sys_construction = gensym("###__is_system_construction###")
    kws = call.args[2].args
    for (i, kw) in enumerate(kws)
        if Meta.isexpr(kw, (:(=), :kw))
            kw.args[2] = :($is_sys_construction ? $(kw.args[2]) :
                           $default_to_parentscope($(kw.args[2])))
        elseif kw isa Symbol
            rhs = :($is_sys_construction ? $(kw) : $default_to_parentscope($(kw)))
            kws[i] = Expr(:kw, kw, rhs)
        end
    end
    if !any(kw -> (kw isa Symbol ? kw : kw.args[1]) == :name, kws)
        pushfirst!(kws, Expr(:kw, :name, runtime ? name : Meta.quot(name)))
    end
    op = call.args[1]
    quote
        $is_sys_construction = ($op isa $DataType) && ($op <: $AbstractSystem)
        $call
    end
end
function _named_idxs(name::Symbol, idxs, call; extra_args = "")
    if call.head !== :->
        throw(ArgumentError("Not an anonymous function"))
    end
    if !isa(call.args[1], Symbol)
        throw(ArgumentError("not a single-argument anonymous function"))
    end
    sym, ex = call.args
    ex = Base.Cartesian.poplinenum(ex)
    ex = _named(:(Symbol($(Meta.quot(name)), :_, $sym)), ex, true)
    ex = Base.Cartesian.poplinenum(ex)
    :($name = map($sym -> begin
            $extra_args
            $ex
        end, $idxs))
end
function setname(x, name)
    @set x.name = name
end
function single_named_expr(expr)
    name, call = split_assign(expr)
    if Meta.isexpr(name, :ref)
        name, idxs = name.args
        check_name(name)
        var = gensym(name)
        ex = quote
            $var = $(_named(name, call))
            $name = map(i -> $setname($var, Symbol($(Meta.quot(name)), :_, i)), $idxs)
        end
        ex
    else
        check_name(name)
        :($name = $(_named(name, call)))
    end
end
function named_expr(expr)
    if Meta.isexpr(expr, :block)
        newexpr = Expr(:block)
        names = Expr(:vcat)
        for ex in expr.args
            ex isa LineNumberNode && continue
            push!(newexpr.args, single_named_expr(ex))
            push!(names.args, ex.args[1])
        end
        push!(newexpr.args, names)
        newexpr
    else
        single_named_expr(expr)
    end
end
function check_name(name)
    name isa Symbol ||
        throw(Meta.ParseError("The lhs must be a symbol (a) or a ref (a[1:10]). Got $name."))
end
""""""
macro named(expr)
    esc(named_expr(expr))
end
macro named(name::Symbol, idxs, call)
    esc(_named_idxs(name, idxs, call))
end
function default_to_parentscope(v)
    uv = unwrap(v)
    uv isa SymbolicT || return v
    apply_to_variables(v) do sym
        ParentScope(sym)
    end
end
function _config(expr, namespace)
    cn = Base.Fix2(_config, namespace)
    if Meta.isexpr(expr, :.)
        return :($getproperty($(map(cn, expr.args)...); namespace = $namespace))
    elseif Meta.isexpr(expr, :function)
        def = splitdef(expr)
        def[:args] = map(cn, def[:args])
        def[:body] = cn(def[:body])
        combinedef(def)
    elseif expr isa Expr && !isempty(expr.args)
        return Expr(expr.head, map(cn, expr.args)...)
    elseif Meta.isexpr(expr, :(=))
        return Expr(:(=), map(cn, expr.args)...)
    else
        expr
    end
end
""""""
macro nonamespace(expr)
    esc(_config(expr, false))
end
""""""
macro namespace(expr)
    esc(_config(expr, true))
end
function component_post_processing(expr, isconnector)
    @assert expr isa Expr && (expr.head == :function || (expr.head == :(=) &&
              expr.args[1] isa Expr &&
              expr.args[1].head == :call))
    sig = expr.args[1]
    body = expr.args[2]
    fname = sig.args[1]
    args = sig.args[2:end]
    quote
        $Base.@__doc__ function $fname($(args...))
            res = (() -> $body)()
            if $isdefined(res, :gui_metadata) && $getfield(res, :gui_metadata) === nothing
                name = $(Meta.quot(fname))
                if $isconnector
                    $Setfield.@set!(res.connector_type=$connector_type(res))
                end
                $Setfield.@set!(res.gui_metadata=$GUIMetadata($GlobalRef(
                    @__MODULE__, name)))
            else
                res
            end
        end
    end
end
""""""
macro component(expr)
    esc(component_post_processing(expr, false))
end
""""""
macro mtkcompile(exprs...)
    expr = exprs[1]
    named_expr = ModelingToolkit.named_expr(expr)
    name = named_expr.args[1]
    kwargs = Base.tail(exprs)
    kwargs = map(kwargs) do ex
        @assert ex.head == :(=)
        Expr(:kw, ex.args[1], ex.args[2])
    end
    if isempty(kwargs)
        kwargs = ()
    else
        kwargs = (Expr(:parameters, kwargs...),)
    end
    call_expr = Expr(:call, mtkcompile, kwargs..., name)
    esc(quote
        $named_expr
        $name = $call_expr
    end)
end
""""""
function debug_system(
        sys::AbstractSystem; functions = [log, sqrt, (^), /, inv, asin, acos], kw...)
    if !(functions isa Set)
        functions = Set(functions)
    end
    if has_systems(sys) && !isempty(get_systems(sys))
        error("debug_system(sys) only works on systems with no sub-systems! Consider flattening it with flatten(sys) or mtkcompile(sys) first.")
    end
    if has_eqs(sys)
        eqs = debug_sub.(equations(sys), Ref(functions); kw...)
        @set! sys.eqs = eqs
        @set! sys.ps = unique!([get_ps(sys); ASSERTION_LOG_VARIABLE])
        @set! sys.defaults = merge(get_defaults(sys), Dict(ASSERTION_LOG_VARIABLE => true))
    end
    if has_observed(sys)
        @set! sys.observed = debug_sub.(observed(sys), Ref(functions); kw...)
    end
    if iscomplete(sys)
        sys = complete(sys; split = is_split(sys))
    end
    return sys
end
@latexrecipe function f(sys::AbstractSystem)
    return latexify(equations(sys))
end
function Base.show(io::IO, ::MIME"text/latex", x::AbstractSystem)
    print(io, "\$\$ " * latexify(x) * " \$\$")
end
struct InvalidSystemException <: Exception
    msg::String
end
function Base.showerror(io::IO, e::InvalidSystemException)
    print(io, "InvalidSystemException: ", e.msg)
end
struct ExtraVariablesSystemException <: Exception
    msg::String
end
function Base.showerror(io::IO, e::ExtraVariablesSystemException)
    println(io, "ExtraVariablesSystemException: ", e.msg)
    print(io,
        "Note that the process of determining extra variables is a best-effort heuristic. " *
        "The true extra variables are dependent on the model and may not be in this list.")
end
struct ExtraEquationsSystemException <: Exception
    msg::String
end
function Base.showerror(io::IO, e::ExtraEquationsSystemException)
    println(io, "ExtraEquationsSystemException: ", e.msg)
    print(io,
        "Note that the process of determining extra equations is a best-effort heuristic. " *
        "The true extra equations are dependent on the model and may not be in this list.")
end
struct HybridSystemNotSupportedException <: Exception
    msg::String
end
function Base.showerror(io::IO, e::HybridSystemNotSupportedException)
    print(io, "HybridSystemNotSupportedException: ", e.msg)
end
function AbstractTrees.children(sys::AbstractSystem)
    ModelingToolkit.get_systems(sys)
end
function AbstractTrees.printnode(
        io::IO, sys::AbstractSystem; describe = false, bold = false)
    printstyled(io, nameof(sys); bold)
    describe && !isempty(description(sys)) && print(io, ": ", description(sys))
end
""""""
function hierarchy(sys::AbstractSystem; describe = false, bold = describe, kwargs...)
    print_tree(sys; printnode_kw = (describe = describe, bold = bold), kwargs...)
end
function Base.IteratorEltype(::Type{<:TreeIterator{ModelingToolkit.AbstractSystem}})
    Base.HasEltype()
end
function Base.eltype(::Type{<:TreeIterator{ModelingToolkit.AbstractSystem}})
    ModelingToolkit.AbstractSystem
end
function check_array_equations_unknowns(eqs, dvs)
    if any(eq -> eq isa Equation && Symbolics.isarraysymbolic(eq.lhs), eqs)
        throw(ArgumentError("The system has array equations. Call `mtkcompile` to handle such equations or scalarize them manually."))
    end
    if any(x -> Symbolics.isarraysymbolic(x), dvs)
        throw(ArgumentError("The system has array unknowns. Call `mtkcompile` to handle this or scalarize them manually."))
    end
end
function check_eqs_u0(eqs, dvs, u0; check_length = true, kwargs...)
    if u0 !== nothing
        if check_length
            if !(length(eqs) == length(dvs) == length(u0))
                throw(ArgumentError("Equations ($(length(eqs))), unknowns ($(length(dvs))), and initial conditions ($(length(u0))) are of different lengths."))
            end
        elseif length(dvs) != length(u0)
            throw(ArgumentError("Unknowns ($(length(dvs))) and initial conditions ($(length(u0))) are of different lengths."))
        end
    elseif check_length && (length(eqs) != length(dvs))
        throw(ArgumentError("Equations ($(length(eqs))) and Unknowns ($(length(dvs))) are of different lengths."))
    end
    return nothing
end
""""""
function extend(sys::AbstractSystem, basesys::AbstractSystem;
        name::Symbol = nameof(sys), description = description(sys),
        gui_metadata = get_gui_metadata(sys))
    T = SciMLBase.parameterless_type(basesys)
    ivs = independent_variables(basesys)
    if !(sys isa T)
        if length(ivs) == 0
            sys = convert_system(T, sys)
        elseif length(ivs) == 1
            sys = convert_system(T, sys, ivs[1])
        else
            throw("Extending multivariate systems is not supported")
        end
    end
    eqs = union(get_eqs(basesys), get_eqs(sys))
    sts = union(get_unknowns(basesys), get_unknowns(sys))
    ps = union(get_ps(basesys), get_ps(sys))
    dep_ps = union(get_parameter_dependencies(basesys), get_parameter_dependencies(sys))
    obs = union(get_observed(basesys), get_observed(sys))
    cevs = union(get_continuous_events(basesys), get_continuous_events(sys))
    devs = union(get_discrete_events(basesys), get_discrete_events(sys))
    defs = merge(get_defaults(basesys), get_defaults(sys))
    meta = MetadataT()
    for kvp in get_metadata(basesys)
        kvp[1] == MutableCacheKey && continue
        meta = Base.ImmutableDict(meta, kvp)
    end
    for kvp in get_metadata(sys)
        kvp[1] == MutableCacheKey && continue
        meta = Base.ImmutableDict(meta, kvp)
    end
    syss = union(get_systems(basesys), get_systems(sys))
    args = length(ivs) == 0 ? (eqs, sts, ps) : (eqs, ivs[1], sts, ps)
    kwargs = (observed = obs, continuous_events = cevs,
        discrete_events = devs, defaults = defs, systems = syss, metadata = meta,
        name = name, description = description, gui_metadata = gui_metadata)
    ieqs = union(get_initialization_eqs(basesys), get_initialization_eqs(sys))
    guesses = merge(get_guesses(basesys), get_guesses(sys))
    kwargs = merge(kwargs, (initialization_eqs = ieqs, guesses = guesses))
    if has_assertions(basesys)
        kwargs = merge(
            kwargs, (; assertions = merge(get_assertions(basesys), get_assertions(sys))))
    end
    newsys = T(args...; kwargs...)
    @set! newsys.parameter_dependencies = dep_ps
    return newsys
end
""""""
function extend(sys, basesys::Vector{T}) where {T <: AbstractSystem}
    foldl(extend, basesys, init = sys)
end
""""""
function Base.:(&)(sys::AbstractSystem, basesys::AbstractSystem; kwargs...)
    extend(sys, basesys; kwargs...)
end
""""""
function Base.:(&)(
        sys::AbstractSystem, basesys::Vector{T}; kwargs...) where {T <: AbstractSystem}
    extend(sys, basesys; kwargs...)
end
""""""
function compose(sys::AbstractSystem, systems::AbstractArray; name = nameof(sys))
    nsys = length(systems)
    nsys == 0 && return sys
    @set! sys.name = name
    @set! sys.systems = [get_systems(sys); systems]
    if has_is_dde(sys)
        @set! sys.is_dde = _check_if_dde(equations(sys), get_iv(sys), get_systems(sys))
    end
    newunknowns = OrderedSet{SymbolicT}()
    newparams = OrderedSet{SymbolicT}()
    iv = has_iv(sys) ? get_iv(sys) : nothing
    for ssys in systems
        collect_scoped_vars!(newunknowns, newparams, ssys, iv)
    end
    @set! sys.unknowns = unique!(vcat(get_unknowns(sys), collect(newunknowns)))
    @set! sys.ps = unique!(vcat(get_ps(sys), collect(newparams)))
    return sys
end
""""""
function compose(syss...; name = nameof(first(syss)))
    compose(first(syss), collect(syss[2:end]); name = name)
end
""""""
Base.:(∘)(sys1::AbstractSystem, sys2::AbstractSystem) = compose(sys1, sys2)
function UnPack.unpack(sys::ModelingToolkit.AbstractSystem, ::Val{p}) where {p}
    getproperty(sys, p; namespace = false)
end
""""""
function missing_variable_defaults(
        sys::AbstractSystem, default = 0.0; subset = unknowns(sys))
    varmap = get_defaults(sys)
    varmap = Dict(Symbolics.diff2term(value(k)) => value(varmap[k]) for k in keys(varmap))
    missingvars = setdiff(subset, keys(varmap))
    ds = Pair[]
    n = length(missingvars)
    if default isa Vector
        @assert length(default)==n "`default` size ($(length(default))) should match the number of missing variables: $n"
    end
    for (i, missingvar) in enumerate(missingvars)
        if default isa Vector
            push!(ds, missingvar => default[i])
        else
            push!(ds, missingvar => default)
        end
    end
    return ds
end
keytype(::Type{<:Pair{T, V}}) where {T, V} = T
function Symbolics.substitute(sys::AbstractSystem, rules::Union{Vector{<:Pair}, Dict})
    if has_continuous_domain(sys) && get_continuous_events(sys) !== nothing &&
       !isempty(get_continuous_events(sys)) ||
       has_discrete_events(sys) && get_discrete_events(sys) !== nothing &&
       !isempty(get_discrete_events(sys))
        @warn "`substitute` only supports performing substitutions in equations. This system has events, which will not be updated."
    end
    if keytype(eltype(rules)) <: Symbol
        dict = todict(rules)
        systems = get_systems(sys)
        @set! sys.systems = map(Base.Fix2(substitute, dict), systems)
        something(get(rules, nameof(sys), nothing), sys)
    elseif sys isa System
        rules = todict(map(r -> Symbolics.unwrap(r[1]) => Symbolics.unwrap(r[2]),
            collect(rules)))
        newsys = @set sys.eqs = substitute(get_eqs(sys), rules)
        @set! newsys.unknowns = map(get_unknowns(sys)) do var
            get(rules, var, var)
        end
        @set! newsys.ps = map(get_ps(sys)) do var
            get(rules, var, var)
        end
        @set! newsys.parameter_dependencies = substitute(
            get_parameter_dependencies(sys), rules)
        @set! newsys.defaults = Dict(substitute(k, rules) => substitute(v, rules)
        for (k, v) in get_defaults(sys))
        @set! newsys.guesses = Dict(substitute(k, rules) => substitute(v, rules)
        for (k, v) in get_guesses(sys))
        @set! newsys.noise_eqs = substitute(get_noise_eqs(sys), rules)
        @set! newsys.costs = Vector{Union{Real, BasicSymbolic}}(substitute(
            get_costs(sys), rules))
        @set! newsys.observed = substitute(get_observed(sys), rules)
        @set! newsys.initialization_eqs = substitute(
            get_initialization_eqs(sys), rules)
        @set! newsys.constraints = substitute(get_constraints(sys), rules)
        @set! newsys.systems = map(s -> substitute(s, rules), get_systems(sys))
    else
        error("substituting symbols is not supported for $(typeof(sys))")
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
""""""
function dump_parameters(sys::AbstractSystem)
    defs = defaults(sys)
    pdeps = get_parameter_dependencies(sys)
    metas = map(dump_variable_metadata.(parameters(sys))) do meta
        if haskey(defs, meta.var)
            meta = merge(meta, (; default = defs[meta.var]))
        end
        meta
    end
    pdep_metas = map(pdeps) do eq
        sym = eq.lhs
        val = eq.rhs
        meta = dump_variable_metadata(sym)
        defs[eq.lhs] = eq.rhs
        meta = merge(meta,
            (; dependency = val,
                default = symbolic_evaluate(val, defs)))
        return meta
    end
    return vcat(metas, pdep_metas)
end
""""""
function dump_unknowns(sys::AbstractSystem)
    defs = add_toterms(defaults(sys))
    gs = add_toterms(guesses(sys))
    map(dump_variable_metadata.(unknowns(sys))) do meta
        if haskey(defs, meta.var)
            meta = merge(meta, (; default = defs[meta.var]))
        end
        if haskey(gs, meta.var)
            meta = merge(meta, (; guess = gs[meta.var]))
        end
        meta
    end
end
""""""
function parse_variable(sys::AbstractSystem, str::AbstractString)
    iv = has_iv(sys) ? string(getname(get_iv(sys))) : nothing
    str = strip(str)
    derivative_level = 0
    while ((cond1 = startswith(str, "D(")) || startswith(str, "Differential(")) &&
        endswith(str, ")")
        if cond1
            derivative_level += 1
            str = _string_view_inner(str, 2, 1)
            continue
        end
        _tmpstr = _string_view_inner(str, 13, 1)
        if !startswith(_tmpstr, "$iv)(")
            throw(ArgumentError("Expected differential with respect to independent variable $iv in $str"))
        end
        derivative_level += 1
        str = _string_view_inner(_tmpstr, length(iv) + 2, 0)
    end
    arr_idxs = nothing
    if endswith(str, ']')
        open_idx = only(findfirst('[', str))
        idxs_range = nextind(str, open_idx):prevind(str, lastindex(str))
        idxs_str = view(str, idxs_range)
        str = view(str, firstindex(str):prevind(str, open_idx))
        arr_idxs = map(Base.Fix1(parse, Int), eachsplit(idxs_str, ","))
    end
    if iv !== nothing && endswith(str, "($iv)")
        str = _string_view_inner(str, 0, 2 + length(iv))
    end
    dummyderivative_level = 0
    if iv !== nothing && (dd_idx = findfirst('ˍ', str)) !== nothing
        t_idx = findnext(iv, str, dd_idx)
        while t_idx !== nothing
            dummyderivative_level += 1
            t_idx = findnext(iv, str, nextind(str, last(t_idx)))
        end
        str = view(str, firstindex(str):prevind(str, dd_idx))
    end
    if iv !== nothing && endswith(str, "($iv)")
        str = _string_view_inner(str, 0, 2 + length(iv))
    end
    cur = sys
    for ident in eachsplit(str, ('.', NAMESPACE_SEPARATOR))
        ident = Symbol(ident)
        hasproperty(cur, ident) ||
            throw(ArgumentError("System $(nameof(cur)) does not have a subsystem/variable named $(ident)"))
        cur = getproperty(cur, ident)
    end
    if arr_idxs !== nothing
        cur = cur[arr_idxs...]
    end
    for i in 1:(derivative_level + dummyderivative_level)
        cur = Differential(get_iv(sys))(cur)
    end
    return cur
end
function _string_view_inner(str, startoffset, endoffset)
    view(str,
        nextind(str, firstindex(str), startoffset):prevind(str, lastindex(str), endoffset))
end
""""""
function is_diff_equation(eq)
    (eq isa Equation) || (return false)
    isdefined(eq, :lhs) && recursive_hasoperator(Union{Differential, Shift}, eq.lhs) &&
        (return true)
    isdefined(eq, :rhs) && recursive_hasoperator(Union{Differential, Shift}, eq.rhs) &&
        (return true)
    return false
end
""""""
function is_alg_equation(eq)
    return (eq isa Equation) && !is_diff_equation(eq)
end
""""""
alg_equations(sys::AbstractSystem) = filter(is_alg_equation, equations(sys))
""""""
diff_equations(sys::AbstractSystem) = filter(is_diff_equation, equations(sys))
""""""
has_alg_equations(sys::AbstractSystem) = any(is_alg_equation, equations(sys))
""""""
has_diff_equations(sys::AbstractSystem) = any(is_diff_equation, equations(sys))
""""""
get_alg_eqs(sys::AbstractSystem) = filter(is_alg_equation, get_eqs(sys))
""""""
get_diff_eqs(sys::AbstractSystem) = filter(is_diff_equation, get_eqs(sys))
""""""
has_alg_eqs(sys::AbstractSystem) = any(is_alg_equation, get_eqs(sys))
""""""
has_diff_eqs(sys::AbstractSystem) = any(is_diff_equation, get_eqs(sys))
""""""
function validate_replacement_rule(
        rule::Pair{T, T}; namespace = []) where {T <: AbstractSystem}
    lhs, rhs = rule
    iscomplete(lhs) && throw(ArgumentError("LHS of replacement rule cannot be completed."))
    iscomplete(rhs) && throw(ArgumentError("RHS of replacement rule cannot be completed."))
    rhs_h = namespace_hierarchy(nameof(rhs))
    if length(rhs_h) != 1
        throw(ArgumentError("RHS of replacement rule must not be namespaced."))
    end
    rhs_h[1] == namespace_hierarchy(nameof(lhs))[end] ||
        throw(ArgumentError("LHS and RHS must have the same name."))
    if !isequal(get_iv(lhs), get_iv(rhs))
        throw(ArgumentError("LHS and RHS of replacement rule must have the same independent variable."))
    end
    lhs_u = get_unknowns(lhs)
    rhs_u = Dict(get_unknowns(rhs) .=> nothing)
    for u in lhs_u
        if !haskey(rhs_u, u)
            if isempty(namespace)
                throw(ArgumentError("RHS of replacement rule does not contain unknown $u."))
            else
                throw(ArgumentError("Subsystem $(join([namespace; nameof(lhs)], NAMESPACE_SEPARATOR)) of RHS does not contain unknown $u."))
            end
        end
        ru = getkey(rhs_u, u, nothing)
        name = join([namespace; nameof(lhs); (hasname(u) ? getname(u) : Symbol(u))],
            NAMESPACE_SEPARATOR)
        l_connect = something(getconnect(u), Equality)
        r_connect = something(getconnect(ru), Equality)
        if l_connect != r_connect
            throw(ArgumentError("Variable $(name) should have connection metadata $(l_connect),"))
        end
        l_input = isinput(u)
        r_input = isinput(ru)
        if l_input != r_input
            throw(ArgumentError("Variable $name has differing causality. Marked as `input = $l_input` in LHS and `input = $r_input` in RHS."))
        end
        l_output = isoutput(u)
        r_output = isoutput(ru)
        if l_output != r_output
            throw(ArgumentError("Variable $name has differing causality. Marked as `output = $l_output` in LHS and `output = $r_output` in RHS."))
        end
    end
    lhs_p = get_ps(lhs)
    rhs_p = Set(get_ps(rhs))
    for p in lhs_p
        if !(p in rhs_p)
            if isempty(namespace)
                throw(ArgumentError("RHS of replacement rule does not contain parameter $p"))
            else
                throw(ArgumentError("Subsystem $(join([namespace; nameof(lhs)], NAMESPACE_SEPARATOR)) of RHS does not contain parameter $p."))
            end
        end
    end
    lhs_s = get_systems(lhs)
    rhs_s = Dict(nameof(s) => s for s in get_systems(rhs))
    for s in lhs_s
        if haskey(rhs_s, nameof(s))
            rs = rhs_s[nameof(s)]
            if isconnector(s)
                name = join([namespace; nameof(lhs); nameof(s)], NAMESPACE_SEPARATOR)
                if !isconnector(rs)
                    throw(ArgumentError("Subsystem $name of RHS is not a connector."))
                end
                if (lct = get_connector_type(s)) !== (rct = get_connector_type(rs))
                    throw(ArgumentError("Subsystem $name of RHS has connection type $rct but LHS has $lct."))
                end
            end
            validate_replacement_rule(s => rs; namespace = [namespace; nameof(rhs)])
            continue
        end
        name1 = join([namespace; nameof(lhs)], NAMESPACE_SEPARATOR)
        throw(ArgumentError("$name1 of replacement rule does not contain subsystem $(nameof(s))."))
    end
end
""""""
function recursive_getproperty(
        root::AbstractSystem, hierarchy::Vector{Symbol}; skip_namespace_first = true)
    cur = root
    for (i, name) in enumerate(hierarchy)
        cur = getproperty(cur, name; namespace = i > 1 || !skip_namespace_first)
    end
    return unwrap(cur)
end
""""""
function recreate_connections(sys::AbstractSystem)
    eqs = map(get_eqs(sys)) do eq
        eq.lhs isa Union{Connection, AnalysisPoint} || return eq
        if eq.lhs isa Connection
            oldargs = get_systems(eq.rhs)
        else
            ap::AnalysisPoint = eq.rhs
            oldargs = [ap.input; ap.outputs]
        end
        newargs = map(get_systems(eq.rhs)) do arg
            rewrap_nameof = arg isa SymbolicWithNameof
            if rewrap_nameof
                arg = arg.var
            end
            name = arg isa AbstractSystem ? nameof(arg) : getname(arg)
            hierarchy = namespace_hierarchy(name)
            newarg = recursive_getproperty(sys, hierarchy)
            if rewrap_nameof
                newarg = SymbolicWithNameof(newarg)
            end
            return newarg
        end
        if eq.lhs isa Connection
            return eq.lhs ~ Connection(newargs)
        else
            return eq.lhs ~ AnalysisPoint(newargs[1], eq.rhs.name, newargs[2:end])
        end
    end
    @set! sys.eqs = eqs
    @set! sys.systems = map(recreate_connections, get_systems(sys))
    return sys
end
""""""
function substitute_component(sys::T, rule::Pair{T, T}) where {T <: AbstractSystem}
    iscomplete(sys) &&
        throw(ArgumentError("Cannot replace subsystems of completed systems"))
    validate_replacement_rule(rule)
    lhs, rhs = rule
    hierarchy = namespace_hierarchy(nameof(lhs))
    newsys, _ = modify_nested_subsystem(sys, hierarchy) do inner
        return rhs, ()
    end
    return recreate_connections(newsys)
end
