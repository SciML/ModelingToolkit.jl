const SYSTEM_COUNT = Threads.Atomic{UInt}(0)

get_component_type(x::AbstractSystem) = get_gui_metadata(x).type
struct GUIMetadata
    type::GlobalRef
    layout::Any
end

GUIMetadata(type) = GUIMetadata(type, nothing)

"""
```julia
generate_custom_function(sys::AbstractSystem, exprs, dvs = unknowns(sys),
                         ps = parameters(sys); kwargs...)
```

Generate a function to evaluate `exprs`. `exprs` is a symbolic expression or
array of symbolic expression involving symbolic variables in `sys`. The symbolic variables
may be subsetted using `dvs` and `ps`. All `kwargs` are passed to the internal
[`build_function`](@ref) call. The returned function can be called as `f(u, p, t)` or
`f(du, u, p, t)` for time-dependent systems and `f(u, p)` or `f(du, u, p)` for
time-independent systems. If `split=true` (the default) was passed to [`complete`](@ref),
[`mtkcompile`](@ref) or [`@mtkcompile`](@ref), `p` is expected to be an `MTKParameters`
object.
"""
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

const MTKPARAMETERS_ARG = Sym{Vector{Vector}}(:___mtkparameters___)

"""
    $(TYPEDSIGNATURES)

Obtain the name of `sys`.
"""
Base.nameof(sys::AbstractSystem) = getfield(sys, :name)
"""
    $(TYPEDSIGNATURES)

Obtain the description associated with `sys` if present, and an empty
string otherwise.
"""
description(sys::AbstractSystem) = has_description(sys) ? get_description(sys) : ""

"""
$(TYPEDSIGNATURES)

Get the independent variable(s) of the system `sys`.

See also [`@independent_variables`](@ref) and [`ModelingToolkit.get_iv`](@ref).
"""
function independent_variables(sys::AbstractSystem)
    if isdefined(sys, :iv) && getfield(sys, :iv) !== nothing
        return [getfield(sys, :iv)]
    elseif isdefined(sys, :ivs)
        return getfield(sys, :ivs)
    else
        return []
    end
end

#Treat the result as a vector of symbols always
function SymbolicIndexingInterface.is_variable(sys::AbstractSystem, sym)
    sym = unwrap(sym)
    if sym isa Int    # [x, 1] coerces 1 to a Num
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

function SymbolicIndexingInterface.is_parameter(sys::AbstractSystem, sym)
    sym = unwrap(sym)
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        return sym isa ParameterIndex || is_parameter(ic, sym) ||
               iscall(sym) &&
               operation(sym) === getindex &&
               is_parameter(ic, first(arguments(sym)))
    end
    if unwrap(sym) isa Int
        return unwrap(sym) in 1:length(parameter_symbols(sys))
    end
    return any(isequal(sym), parameter_symbols(sys)) ||
           hasname(sym) && !(iscall(sym) && operation(sym) == getindex) &&
           is_parameter(sys, getname(sym))
end

function SymbolicIndexingInterface.is_parameter(sys::AbstractSystem, sym::Symbol)
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        return is_parameter(ic, sym)
    end

    named_parameters = [getname(x)
                        for x in parameter_symbols(sys)
                        if hasname(x) && !(iscall(x) && operation(x) == getindex)]
    return any(isequal(sym), named_parameters) ||
           count(NAMESPACE_SEPARATOR, string(sym)) == 1 &&
           count(isequal(sym),
        Symbol.(nameof(sys), NAMESPACE_SEPARATOR_SYMBOL, named_parameters)) == 1
end

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

"""
    $(TYPEDSIGNATURES)

Check if the system `sys` contains an observed equation with LHS `sym`.
"""
function has_observed_with_lhs(sys::AbstractSystem, sym)
    has_observed(sys) || return false
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        return haskey(ic.observed_syms_to_timeseries, sym)
    else
        return any(isequal(sym), observables(sys))
    end
end

"""
    $(TYPEDSIGNATURES)

Check if the system `sys` contains a parameter dependency equation with LHS `sym`.
"""
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
# Need this to avoid ambiguity with the array case
for traitT in [
    ScalarSymbolic,
    ArraySymbolic
]
    @eval function _all_ts_idxs!(ts_idxs, ::$traitT, sys, sym)
        allsyms = vars(sym; op = Symbolics.Operator)
        for s in allsyms
            s = unwrap(s)
            if is_variable(sys, s) || is_independent_variable(sys, s)
                push!(ts_idxs, ContinuousTimeseries())
            elseif is_timeseries_parameter(sys, s)
                push!(ts_idxs, timeseries_parameter_index(sys, s).timeseries_idx)
            elseif is_time_dependent(sys) && iscall(s) && issym(operation(s)) &&
                   length(arguments(s)) == 1 && is_variable(sys, operation(s)(get_iv(sys)))
                # DDEs case, to detect x(t - k)
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
                    # for split=false systems
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

"""
    $(TYPEDSIGNATURES)

Check whether a system is marked as `complete`.
"""
iscomplete(sys::AbstractSystem) = isdefined(sys, :complete) && getfield(sys, :complete)
"""
    $(TYPEDSIGNATURES)

Check whether a system performs namespacing.
"""
function does_namespacing(sys::AbstractSystem)
    if isdefined(sys, :namespacing)
        getfield(sys, :namespacing)
    else
        !iscomplete(sys)
    end
end

"""
$(TYPEDSIGNATURES)

If a system is scheduled, then changing its equations, variables, and
parameters is no longer legal.
"""
function isscheduled(sys::AbstractSystem)
    if has_schedule(sys)
        get_schedule(sys) !== nothing
    elseif has_isscheduled(sys)
        get_isscheduled(sys)
    else
        false
    end
end

"""
    Initial(x)

The `Initial` operator. Used by initialization to store constant constraints on variables
of a system. See the documentation section on initialization for more information.
"""
struct Initial <: Symbolics.Operator end
is_timevarying_operator(::Type{Initial}) = false
Initial(x) = Initial()(x)
SymbolicUtils.promote_symtype(::Type{Initial}, T) = T
SymbolicUtils.isbinop(::Initial) = false
Base.nameof(::Initial) = :Initial
Base.show(io::IO, x::Initial) = print(io, "Initial")

function (f::Initial)(x)
    # wrap output if wrapped input
    iw = Symbolics.iswrapped(x)
    x = unwrap(x)
    # non-symbolic values don't change
    if symbolic_type(x) == NotSymbolic()
        return x
    end
    # differential variables are default-toterm-ed
    if iscall(x) && operation(x) isa Union{Differential, Shift}
        x = default_toterm(x)
    end
    # don't double wrap
    iscall(x) && operation(x) isa Initial && return x
    result = if symbolic_type(x) == ArraySymbolic()
        # create an array for `Initial(array)`
        Symbolics.array_term(f, x)
    elseif iscall(x) && operation(x) == getindex
        # instead of `Initial(x[1])` create `Initial(x)[1]`
        # which allows parameter indexing to handle this case automatically.
        arr = arguments(x)[1]
        term(getindex, f(arr), arguments(x)[2:end]...)
    else
        term(f, x)
    end
    # the result should be a parameter
    result = toparam(result)
    if iw
        result = wrap(result)
    end
    return result
end

# This is required so `fast_substitute` works
function SymbolicUtils.maketerm(::Type{<:BasicSymbolic}, ::Initial, args, meta)
    val = Initial()(args...)
    if symbolic_type(val) == NotSymbolic()
        return val
    end
    return metadata(val, meta)
end

supports_initialization(sys::AbstractSystem) = true

function add_initialization_parameters(sys::AbstractSystem; split = true)
    @assert !has_systems(sys) || isempty(get_systems(sys))
    supports_initialization(sys) || return sys
    is_initializesystem(sys) && return sys

    all_initialvars = Set{BasicSymbolic}()
    # time-independent systems don't initialize unknowns
    # but may initialize parameters using guesses for unknowns
    eqs = equations(sys)
    if !(eqs isa Vector{Equation})
        eqs = Equation[x for x in eqs if x isa Equation]
    end
    obs, eqs = unhack_observed(observed(sys), eqs)
    for x in Iterators.flatten((unknowns(sys), Iterators.map(eq -> eq.lhs, obs)))
        x = unwrap(x)
        if iscall(x) && operation(x) == getindex && split
            push!(all_initialvars, arguments(x)[1])
        else
            push!(all_initialvars, x)
        end
    end

    # add derivatives of all variables for steady-state initial conditions
    if is_time_dependent(sys) && !is_discrete_system(sys)
        D = Differential(get_iv(sys))
        union!(all_initialvars, [D(v) for v in all_initialvars if iscall(v)])
    end
    for eq in get_parameter_dependencies(sys)
        is_variable_floatingpoint(eq.lhs) || continue
        push!(all_initialvars, eq.lhs)
    end
    all_initialvars = collect(all_initialvars)
    initials = map(Initial(), all_initialvars)
    @set! sys.ps = unique!([get_ps(sys); initials])
    defs = copy(get_defaults(sys))
    for ivar in initials
        if symbolic_type(ivar) == ScalarSymbolic()
            defs[ivar] = false
        else
            defs[ivar] = collect(ivar)
            for scal_ivar in defs[ivar]
                defs[scal_ivar] = false
            end
        end
    end
    @set! sys.defaults = defs
    return sys
end

"""
Returns true if the parameter `p` is of the form `Initial(x)`.
"""
function isinitial(p)
    p = unwrap(p)
    return iscall(p) && (operation(p) isa Initial ||
            operation(p) === getindex && isinitial(arguments(p)[1]))
end

"""
    $(TYPEDSIGNATURES)

Find [`GlobalScope`](@ref)d variables in `sys` and add them to the unknowns/parameters.
"""
function discover_globalscoped(sys::AbstractSystem)
    newunknowns = OrderedSet()
    newparams = OrderedSet()
    iv = has_iv(sys) ? get_iv(sys) : nothing
    collect_scoped_vars!(newunknowns, newparams, sys, iv; depth = -1)
    setdiff!(newunknowns, observables(sys))
    @set! sys.ps = unique!(vcat(get_ps(sys), collect(newparams)))
    @set! sys.unknowns = unique!(vcat(get_unknowns(sys), collect(newunknowns)))
    return sys
end

"""
$(TYPEDSIGNATURES)

Mark a system as completed. A completed system is a system which is done being
defined/modified and is ready for structural analysis or other transformations.
This allows for analyses and optimizations to be performed which require knowing
the global structure of the system.

One property to note is that if a system is complete, the system will no longer
namespace its subsystems or variables, i.e. `isequal(complete(sys).v.i, v.i)`.

This namespacing functionality can also be toggled independently of `complete`
using [`toggle_namespacing`](@ref).
"""
function complete(
        sys::AbstractSystem; split = true, flatten = true, add_initial_parameters = true)
    sys = discover_globalscoped(sys)

    if flatten
        eqs = equations(sys)
        if eqs isa AbstractArray && eltype(eqs) <: Equation
            newsys = expand_connections(sys)
        else
            newsys = sys
        end
        newsys = ModelingToolkit.flatten(newsys)
        if has_parent(newsys) && get_parent(sys) === nothing
            @set! newsys.parent = complete(sys; split = false, flatten = false)
        end
        sys = newsys
        sys = process_parameter_equations(sys)
        if add_initial_parameters
            sys = add_initialization_parameters(sys; split)
        end
        if has_continuous_events(sys) && is_time_dependent(sys)
            @set! sys.continuous_events = complete.(
                get_continuous_events(sys); iv = get_iv(sys),
                alg_eqs = [alg_equations(sys); observed(sys)])
        end
        if has_discrete_events(sys) && is_time_dependent(sys)
            @set! sys.discrete_events = complete.(
                get_discrete_events(sys); iv = get_iv(sys),
                alg_eqs = [alg_equations(sys); observed(sys)])
        end
    end
    if split && has_index_cache(sys)
        @set! sys.index_cache = IndexCache(sys)
        # Ideally we'd do `get_ps` but if `flatten = false`
        # we don't get all of them. So we call `parameters`.
        all_ps = parameters(sys; initial_parameters = true)
        # inputs have to be maintained in a specific order
        input_vars = inputs(sys)
        if !isempty(all_ps)
            # reorder parameters by portions
            ps_split = reorder_parameters(sys, all_ps)
            # if there are tunables, they will all be in `ps_split[1]`
            # and the arrays will have been scalarized
            ordered_ps = eltype(all_ps)[]
            # if there are no tunables, vcat them
            if !isempty(get_index_cache(sys).tunable_idx)
                unflatten_parameters!(ordered_ps, ps_split[1], all_ps)
                ps_split = Base.tail(ps_split)
            end
            # unflatten initial parameters
            if !isempty(get_index_cache(sys).initials_idx)
                unflatten_parameters!(ordered_ps, ps_split[1], all_ps)
                ps_split = Base.tail(ps_split)
            end
            ordered_ps = vcat(
                ordered_ps, reduce(vcat, ps_split; init = eltype(ordered_ps)[]))
            if isscheduled(sys)
                # ensure inputs are sorted
                input_idxs = findfirst.(isequal.(input_vars), (ordered_ps,))
                @assert all(!isnothing, input_idxs)
                @assert issorted(input_idxs)
            end
            @set! sys.ps = ordered_ps
        end
    elseif has_index_cache(sys)
        @set! sys.index_cache = nothing
    end
    if isdefined(sys, :initializesystem) && get_initializesystem(sys) !== nothing
        @set! sys.initializesystem = complete(get_initializesystem(sys); split)
    end
    sys = toggle_namespacing(sys, false; safe = true)
    isdefined(sys, :complete) ? (@set! sys.complete = true) : sys
end

"""
    $(TYPEDSIGNATURES)

Return a new `sys` with namespacing enabled or disabled, depending on `value`. The
keyword argument `safe` denotes whether systems that do not support such a toggle
should error or be ignored.
"""
function toggle_namespacing(sys::AbstractSystem, value::Bool; safe = false)
    if !isdefined(sys, :namespacing)
        safe && return sys
        throw(ArgumentError("The system must define the `namespacing` flag to toggle namespacing"))
    end
    @set sys.namespacing = value
end

"""
    $(TYPEDSIGNATURES)

Given a flattened array of parameters `params` and a collection of all (unscalarized)
parameters in the system `all_ps`, unscalarize the elements in `params` and append
to `buffer` in the same order as they are present in `params`. Effectively, if
`params = [p[1], p[2], p[3], q]` then this is equivalent to `push!(buffer, p, q)`.
"""
function unflatten_parameters!(buffer, params, all_ps)
    i = 1
    # go through all the tunables
    while i <= length(params)
        sym = params[i]
        # if the sym is not a scalarized array symbolic OR it was already scalarized,
        # just push it as-is
        if !iscall(sym) || operation(sym) != getindex ||
           any(isequal(sym), all_ps)
            push!(buffer, sym)
            i += 1
            continue
        end
        # the next `length(sym)` symbols should be scalarized versions of the same
        # array symbolic
        if !allequal(first(arguments(x))
        for x in view(params, i:(i + length(sym) - 1)))
            error("This should not be possible. Please open an issue in ModelingToolkit.jl with an MWE and stacktrace.")
        end
        arrsym = first(arguments(sym))
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
        """
        $(TYPEDSIGNATURES)

        Get the internal field `$($(QuoteNode(prop)))` of a system `sys`.
        It only includes `$($(QuoteNode(prop)))` local to `sys`; not those of its subsystems,
        like `unknowns(sys)`, `parameters(sys)` and `equations(sys)` does.

        See also [`has_$($(QuoteNode(prop)))`](@ref).
        """
        $fname_get(sys::AbstractSystem) = getfield(sys, $(QuoteNode(prop)))

        """
        $(TYPEDSIGNATURES)

        Returns whether the system `sys` has the internal field `$($(QuoteNode(prop)))`.

        See also [`get_$($(QuoteNode(prop)))`](@ref).
        """
        $fname_has(sys::AbstractSystem) = isdefined(sys, $(QuoteNode(prop)))
    end
end

has_equations(::AbstractSystem) = true

"""
    $(TYPEDSIGNATURES)

Invalidate cached jacobians, etc.
"""
function invalidate_cache!(sys::AbstractSystem)
    has_metadata(sys) || return sys
    empty!(getmetadata(sys, MutableCacheKey, nothing))
    return sys
end

# `::MetadataT` but that is defined later
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
        kwarg = :($(Expr(:kw, :checks, false))) # Inputs should already be checked
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

"""
    Base.getproperty(sys::AbstractSystem, name::Symbol)

Access the subsystem, variable or analysis point of `sys` named `name`. To check if `sys`
will namespace the returned value, use `ModelingToolkit.does_namespacing(sys)`.

See also: [`ModelingToolkit.does_namespacing`](@ref).
"""
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
            if eq.lhs isa AnalysisPoint && nameof(eq.rhs) == name
                return namespace ? renamespace(sys, eq.rhs) : eq.rhs
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

"""
    $(TYPEDSIGNATURES)

Apply function `f` to each variable in expression `ex`. `f` should be a function that takes
a variable and returns the replacement to use. A "variable" in this context refers to a
symbolic quantity created directly from a variable creation macro such as
[`Symbolics.@variables`](@ref), [`@independent_variables`](@ref), [`@parameters`](@ref),
[`@constants`](@ref) or [`@brownians`](@ref).
"""
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

"""
Variable metadata key which contains information about scoping/namespacing of the
variable in a hierarchical system.
"""
abstract type SymScope end

"""
    $(TYPEDEF)

The default scope of a variable. It belongs to the system whose equations it is involved
in and is namespaced by every level of the hierarchy.
"""
struct LocalScope <: SymScope end

"""
    $(TYPEDSIGNATURES)

Apply `LocalScope` to `sym`.
"""
function LocalScope(sym::Union{Num, Symbolic, Symbolics.Arr{Num}})
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

"""
    $(TYPEDEF)

Denotes that the variable does not belong to the system whose equations it is involved
in. It is not namespaced by this system. In the immediate parent of this system, the
scope of this variable is given by `parent`.

# Fields

$(TYPEDFIELDS)
"""
struct ParentScope <: SymScope
    parent::SymScope
end
"""
    $(TYPEDSIGNATURES)

Apply `ParentScope` to `sym`, with `parent` being `LocalScope`.
"""
function ParentScope(sym::Union{Num, Symbolic, Symbolics.Arr{Num}})
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

"""
    $(TYPEDEF)

Denotes that a variable belongs to the root system in the hierarchy, regardless of which
equations of subsystems in the hierarchy it is involved in. Variables with this scope
are never namespaced and only added to the unknowns/parameters of a system when calling
`complete` or `mtkcompile`.
"""
struct GlobalScope <: SymScope end

"""
    $(TYPEDSIGNATURES)

Apply `GlobalScope` to `sym`.
"""
function GlobalScope(sym::Union{Num, Symbolic, Symbolics.Arr{Num}})
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

renamespace(sys, eq::Equation) = namespace_equation(eq, sys)

renamespace(names::AbstractVector, x) = foldr(renamespace, names, init = x)

"""
    $(TYPEDSIGNATURES)

Namespace `x` with the name of `sys`.
"""
function renamespace(sys, x)
    sys === nothing && return x
    x = unwrap(x)
    if x isa Symbolic
        T = typeof(x)
        if iscall(x) && operation(x) isa Operator
            return maketerm(typeof(x), operation(x),
                Any[renamespace(sys, only(arguments(x)))],
                metadata(x))::T
        end
        if iscall(x) && operation(x) === getindex
            args = arguments(x)
            return maketerm(
                typeof(x), operation(x), vcat(renamespace(sys, args[1]), args[2:end]),
                metadata(x))::T
        end
        let scope = getmetadata(x, SymScope, LocalScope())
            if scope isa LocalScope
                rename(x, renamespace(getname(sys), getname(x)))::T
            elseif scope isa ParentScope
                setmetadata(x, SymScope, scope.parent)::T
            else # GlobalScope
                x::T
            end
        end
    elseif x isa AbstractSystem
        rename(x, renamespace(sys, nameof(x)))
    else
        Symbol(getname(sys), NAMESPACE_SEPARATOR_SYMBOL, x)
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

"""
    $(TYPEDSIGNATURES)

Return `equations(sys)`, namespaced by the name of `sys`.
"""
function namespace_equations(sys::AbstractSystem, ivs = independent_variables(sys))
    eqs = equations(sys)
    isempty(eqs) && return Equation[]
    map(eq -> namespace_equation(eq, sys; ivs), eqs)
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

"""
    $(TYPEDSIGNATURES)

Namespace the given equation with the name of the given system `sys`.
"""
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
    return [namespace_jump(j, sys) for j in get_jumps(sys)]
end

function namespace_brownians(sys::AbstractSystem)
    return [renamespace(sys, b) for b in brownians(sys)]
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

function namespace_expr(
        O, sys, n = (sys === nothing ? nothing : nameof(sys));
        ivs = sys === nothing ? nothing : independent_variables(sys))
    sys === nothing && return O
    O = unwrap(O)
    # Exceptions for arrays of symbolic and Ref of a symbolic, the latter
    # of which shows up in broadcasts
    if symbolic_type(O) == NotSymbolic() && !(O isa AbstractArray) && !(O isa Ref)
        return O
    end
    if any(isequal(O), ivs)
        return O
    elseif iscall(O)
        T = typeof(O)
        renamed = let sys = sys, n = n, T = T
            map(a -> namespace_expr(a, sys, n; ivs)::Any, arguments(O))
        end
        if isvariable(O)
            # Use renamespace so the scope is correct, and make sure to use the
            # metadata from the rescoped variable
            rescoped = renamespace(n, O)
            maketerm(typeof(rescoped), operation(rescoped), renamed,
                metadata(rescoped))
        elseif Symbolics.isarraysymbolic(O)
            # promote_symtype doesn't work for array symbolics
            maketerm(typeof(O), operation(O), renamed, metadata(O))
        else
            maketerm(typeof(O), operation(O), renamed, metadata(O))
        end
    elseif isvariable(O)
        renamespace(n, O)
    elseif O isa AbstractArray && is_array_of_symbolics(O)
        let sys = sys, n = n
            map(o -> namespace_expr(o, sys, n; ivs), O)
        end
    else
        O
    end
end

_nonum(@nospecialize x) = x isa Num ? x.val : x

"""
$(TYPEDSIGNATURES)

Get the unknown variables of the system `sys` and its subsystems.
These must be explicitly solved for, unlike `observables(sys)`.

See also [`ModelingToolkit.get_unknowns`](@ref).
"""
function unknowns(sys::AbstractSystem)
    sts = get_unknowns(sys)
    systems = get_systems(sys)
    nonunique_unknowns = if isempty(systems)
        sts
    else
        system_unknowns = reduce(vcat, namespace_variables.(systems))
        isempty(sts) ? system_unknowns : [sts; system_unknowns]
    end
    isempty(nonunique_unknowns) && return nonunique_unknowns
    # `Vector{Any}` is incompatible with the `SymbolicIndexingInterface`, which uses
    # `elsymtype = symbolic_type(eltype(_arg))`
    # which inappropriately returns `NotSymbolic()`
    if nonunique_unknowns isa Vector{Any}
        nonunique_unknowns = _nonum.(nonunique_unknowns)
    end
    @assert typeof(nonunique_unknowns) !== Vector{Any}
    unique(nonunique_unknowns)
end

"""
    unknowns_toplevel(sys::AbstractSystem)

Replicates the behaviour of `unknowns`, but ignores unknowns of subsystems.
"""
function unknowns_toplevel(sys::AbstractSystem)
    if has_parent(sys) && (parent = get_parent(sys)) !== nothing
        return unknowns_toplevel(parent)
    end
    return get_unknowns(sys)
end

"""
$(TYPEDSIGNATURES)

Get the parameters of the system `sys` and its subsystems.

See also [`@parameters`](@ref) and [`ModelingToolkit.get_ps`](@ref).
"""
function parameters(sys::AbstractSystem; initial_parameters = false)
    ps = get_ps(sys)
    if ps == SciMLBase.NullParameters()
        return []
    end
    if eltype(ps) <: Pair
        ps = first.(ps)
    end
    systems = get_systems(sys)
    result = unique(isempty(systems) ? ps :
                    [ps; reduce(vcat, namespace_parameters.(systems))])
    if !initial_parameters && !is_initializesystem(sys)
        filter!(result) do sym
            return !(isoperator(sym, Initial) ||
                     iscall(sym) && operation(sym) == getindex &&
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

"""
    parameters_toplevel(sys::AbstractSystem)

Replicates the behaviour of `parameters`, but ignores parameters of subsystems.
"""
function parameters_toplevel(sys::AbstractSystem)
    if has_parent(sys) && (parent = get_parent(sys)) !== nothing
        return parameters_toplevel(parent)
    end
    return get_ps(sys)
end

"""
    $(TYPEDSIGNATURES)

Get the parameter dependencies of the system `sys` and its subsystems. Requires that the
system is `complete`d.
"""
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

"""
    $(TYPEDSIGNATURES)

Return all of the parameters of the system, including hidden initial parameters and ones
eliminated via `parameter_dependencies`.
"""
function full_parameters(sys::AbstractSystem)
    dep_ps = [eq.lhs for eq in get_parameter_dependencies(sys)]
    vcat(parameters(sys; initial_parameters = true), dep_ps)
end

"""
$(TYPEDSIGNATURES)

Get the assertions for a system `sys` and its subsystems.
"""
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

"""
    $(TYPEDEF)

Information about an `AnalysisPoint` for which the corresponding connection must be
ignored during `expand_connections`, since the analysis point has been transformed.

# Fields

$(TYPEDFIELDS)
"""
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

"""
$(TYPEDSIGNATURES)

Get the guesses for variables in the initialization system of the system `sys` and its subsystems.

See also [`initialization_equations`](@ref) and [`ModelingToolkit.get_guesses`](@ref).
"""
function guesses(sys::AbstractSystem)
    guess = get_guesses(sys)
    systems = get_systems(sys)
    isempty(systems) && return guess
    for subsys in systems
        guess = merge(guess, namespace_guesses(subsys))
    end
    return guess
end

# required in `src/connectors.jl:437`
parameters(_) = []

"""
$(TYPEDSIGNATURES)

Get the observed equations of the system `sys` and its subsystems.
These can be expressed in terms of `unknowns(sys)`, and do not have to be explicitly solved for.

See also [`observables`](@ref) and [`ModelingToolkit.get_observed()`](@ref).
"""
function observed(sys::AbstractSystem)
    obs = get_observed(sys)
    systems = get_systems(sys)
    [obs;
     reduce(vcat,
         (map(o -> namespace_equation(o, s), observed(s)) for s in systems),
         init = Equation[])]
end

"""
$(TYPEDSIGNATURES)

Get the observed variables of the system `sys` and its subsystems.
These can be expressed in terms of `unknowns(sys)`, and do not have to be explicitly solved for.
It is equivalent to all left hand sides of `observed(sys)`.

See also [`observed`](@ref).
"""
function observables(sys::AbstractSystem)
    return map(eq -> eq.lhs, observed(sys))
end

"""
$(TYPEDSIGNATURES)

Get the default values of the system sys and its subsystems.
If they are not explicitly provided, variables and parameters are initialized to these values.

See also [`initialization_equations`](@ref) and [`ModelingToolkit.get_defaults`](@ref).
"""
function defaults(sys::AbstractSystem)
    systems = get_systems(sys)
    defs = get_defaults(sys)
    # `mapfoldr` is really important!!! We should prefer the base model for
    # defaults, because people write:
    #
    # `compose(System(...; defaults=defs), ...)`
    #
    # Thus, right associativity is required and crucial for correctness.
    isempty(systems) ? defs : mapfoldr(namespace_defaults, merge, systems; init = defs)
end

function defaults_and_guesses(sys::AbstractSystem)
    merge(guesses(sys), defaults(sys))
end

unknowns(sys::Union{AbstractSystem, Nothing}, v) = namespace_expr(v, sys)
for vType in [Symbolics.Arr, Symbolics.Symbolic{<:AbstractArray}]
    @eval unknowns(sys::AbstractSystem, v::$vType) = namespace_expr(v, sys)
    @eval parameters(sys::AbstractSystem, v::$vType) = toparam(unknowns(sys, v))
end
parameters(sys::Union{AbstractSystem, Nothing}, v) = toparam(unknowns(sys, v))
for f in [:unknowns, :parameters]
    @eval function $f(sys::AbstractSystem, vs::AbstractArray)
        map(v -> $f(sys, v), vs)
    end
end

flatten(sys::AbstractSystem, args...) = sys

"""
$(TYPEDSIGNATURES)

Get the flattened equations of the system `sys` and its subsystems.
It may include some abbreviations and aliases of observables.
It is often the most useful way to inspect the equations of a system.

See also [`full_equations`](@ref) and [`ModelingToolkit.get_eqs`](@ref).
"""
function equations(sys::AbstractSystem)
    eqs = get_eqs(sys)
    systems = get_systems(sys)
    if isempty(systems)
        return eqs
    else
        eqs = Equation[eqs;
                       reduce(vcat,
                           namespace_equations.(get_systems(sys));
                           init = Equation[])]
        return eqs
    end
end

"""
    equations_toplevel(sys::AbstractSystem)

Replicates the behaviour of `equations`, but ignores equations of subsystems.

Notes:
- Cannot be applied to non-complete systems.
"""
function equations_toplevel(sys::AbstractSystem)
    iscomplete(sys) && error("Cannot apply `equations_toplevel` to complete systems.")
    if has_parent(sys) && (parent = get_parent(sys)) !== nothing
        return equations_toplevel(parent)
    end
    return get_eqs(sys)
end

"""
    $(TYPEDSIGNATURES)

Recursively substitute `dict` into `expr`. Use `Symbolics.simplify` on the expression
if `simplify == true`.
"""
function substitute_and_simplify(expr, dict::AbstractDict, simplify::Bool)
    expr = Symbolics.fixpoint_sub(expr, dict; operator = ModelingToolkit.Initial)
    simplify ? Symbolics.simplify(expr) : expr
end

"""
    $(TYPEDSIGNATURES)

Recursively substitute the observed equations of `sys` into `expr`. If `simplify`, call
`Symbolics.simplify` on the resultant expression.
"""
function substitute_observed(sys::AbstractSystem, expr; simplify = false)
    empty_substitutions(sys) && return expr
    substitutions = get_substitutions(sys)
    return substitute_and_simplify(expr, substitutions, simplify)
end

"""
$(TYPEDSIGNATURES)

Like `equations(sys)`, but also substitutes the observed equations eliminated from the
equations during `mtkcompile`. These equations matches generated numerical code.

See also [`equations`](@ref) and [`ModelingToolkit.get_eqs`](@ref).
"""
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

"""
    $(TYPEDSIGNATURES)

Get the flattened jumps of the system. In other words, obtain all of the jumps in `sys` and
all the subsystems of `sys` (appropriately namespaced).
"""
function jumps(sys::AbstractSystem)
    js = get_jumps(sys)
    systems = get_systems(sys)
    if isempty(systems)
        return js
    end
    return [js; reduce(vcat, namespace_jumps.(systems); init = [])]
end

"""
    $(TYPEDSIGNATURES)

Get all of the brownian variables involved in the system `sys` and all subsystems,
appropriately namespaced.
"""
function brownians(sys::AbstractSystem)
    bs = get_brownians(sys)
    systems = get_systems(sys)
    if isempty(systems)
        return bs
    end
    return [bs; reduce(vcat, namespace_brownians.(systems); init = [])]
end

"""
    $(TYPEDSIGNATURES)

Recursively consolidate the cost vector of `sys` and all subsystems of `sys`, returning the
final scalar cost function.
"""
function cost(sys::AbstractSystem)
    cs = get_costs(sys)
    consolidate = get_consolidate(sys)
    systems = get_systems(sys)
    if isempty(systems)
        return consolidate(cs, Float64[])
    end
    subcosts = [namespace_expr(cost(subsys), subsys) for subsys in systems]
    return consolidate(cs, subcosts)
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
    isempty(cstrs) && return Vector{Union{Equation, Inequality}}(undef, 0)
    map(cstr -> namespace_constraint(cstr, sys), cstrs)
end

"""
    $(TYPEDSIGNATURES)

Get all constraints in the system `sys` and all of its subsystems, appropriately namespaced.
"""
function constraints(sys::AbstractSystem)
    cs = get_constraints(sys)
    systems = get_systems(sys)
    isempty(systems) ? cs : [cs; reduce(vcat, namespace_constraints.(systems))]
end

"""
$(TYPEDSIGNATURES)

Get the initialization equations of the system `sys` and its subsystems.

See also [`guesses`](@ref), [`defaults`](@ref) and [`ModelingToolkit.get_initialization_eqs`](@ref).
"""
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

"""
    $(TYPEDSIGNATURES)

Get the tstops present in `sys` and its subsystems, appropriately namespaced.
"""
function symbolic_tstops(sys::AbstractSystem)
    tstops = get_tstops(sys)
    systems = get_systems(sys)
    isempty(systems) && return tstops
    tstops = [tstops; reduce(vcat, namespace_tstops.(get_systems(sys)); init = [])]
    return tstops
end

"""
    $(TYPEDSIGNATURES)

Obtain the preface associated with `sys` and all of its subsystems, appropriately
namespaced.
"""
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

###
### System utils
###
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

# This is hit because ensemble problems do a deepcopy
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

###
### System I/O
###
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
    Base.remove_linenums!(expr) # keeping the line numbers is never helpful
end

Base.write(io::IO, sys::AbstractSystem) = write(io, readable_code(toexpr(sys)))

"""
    n_expanded_connection_equations(sys::AbstractSystem)

Returns the number of equations that the connections in `sys` expands to.
Equivalent to `length(equations(expand_connections(sys))) - length(filter(eq -> !(eq.lhs isa Connection), equations(sys)))`.
"""
function n_expanded_connection_equations(sys::AbstractSystem)
    # TODO: what about inputs?
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
    limit = get(io, :limit, false) # if output should be limited,
    rows = first(displaysize(io)) ÷ 5 # then allocate ≈1/5 of display height to each list

    # Print name and description
    desc = description(sys)
    name = nameof(sys)
    printstyled(io, "Model ", name, ":"; bold)
    !isempty(desc) && print(io, " ", desc)

    # Print subsystems
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
            maxlen = displaysize(io)[2] - length(name) - 6 # remaining length of line
            if limit && length(desc) > maxlen
                desc = chop(desc, tail = length(desc) - maxlen) * "…" # too long
            end
            print(io, ": ", desc)
        end
    end
    limited = nrows < nsubs
    limited && print(io, "\n  ⋮") # too many to print

    # Print equations
    eqs = equations(sys)
    if eqs isa AbstractArray && eltype(eqs) <: Equation
        neqs = count(eq -> !(eq.lhs isa Connection), eqs)
        next = n_expanded_connection_equations(sys)
        ntot = neqs + next
        ntot > 0 && printstyled(io, "\nEquations ($ntot):"; bold)
        neqs > 0 && print(io, "\n  $neqs standard", hint ? ": see equations($name)" : "")
        next > 0 && print(io, "\n  $next connecting",
            hint ? ": see equations(expand_connections($name))" : "")
        #Base.print_matrix(io, eqs) # usually too long and not useful to print all equations
    end

    # Print variables
    for varfunc in [unknowns, parameters]
        vars = varfunc(sys)
        nvars = length(vars)
        nvars == 0 && continue # skip
        header = titlecase(String(nameof(varfunc))) # e.g. "Unknowns"
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
        limited && printstyled(io, "\n  ⋮") # too many variables to print
    end

    # Print parameter dependencies
    npdeps = has_parameter_dependencies(sys) ? length(get_parameter_dependencies(sys)) : 0
    npdeps > 0 && printstyled(io, "\nParameter dependencies ($npdeps):"; bold)
    npdeps > 0 && hint && print(io, " see parameter_dependencies($name)")

    # Print observed
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
        # canonicalize to use `:parameters`
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

    if !any(kw -> (kw isa Symbol ? kw : kw.args[1]) == :name, kws) # don't overwrite `name` kwarg
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

"""
    @named y = foo(x)
    @named y[1:10] = foo(x)
    @named begin
        y[1:10] = foo(x)
        z = foo(x)
    end # returns `[y; z]`
    @named y 1:10 i -> foo(x*i)  # This is not recommended

Pass the LHS name to the model. When it's calling anything that's not an
AbstractSystem, it wraps all keyword arguments in `default_to_parentscope` so
that namespacing works intuitively when passing a symbolic default into a
component.

Examples:

```julia-repl
julia> using ModelingToolkit

julia> foo(i; name) = (; i, name)
foo (generic function with 1 method)

julia> x = 41
41

julia> @named y = foo(x)
(i = 41, name = :y)

julia> @named y[1:3] = foo(x)
3-element Vector{NamedTuple{(:i, :name), Tuple{Int64, Symbol}}}:
 (i = 41, name = :y_1)
 (i = 41, name = :y_2)
 (i = 41, name = :y_3)
```
"""
macro named(expr)
    esc(named_expr(expr))
end

macro named(name::Symbol, idxs, call)
    esc(_named_idxs(name, idxs, call))
end

function default_to_parentscope(v)
    uv = unwrap(v)
    uv isa Symbolic || return v
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

"""
$(SIGNATURES)

Rewrite `@nonamespace a.b.c` to
`getvar(getvar(a, :b; namespace = false), :c; namespace = false)`.

This is the default behavior of `getvar`. This should be used when inheriting unknowns from a model.
"""
macro nonamespace(expr)
    esc(_config(expr, false))
end

"""
$(SIGNATURES)

Rewrite `@namespace a.b.c` to
`getvar(getvar(a, :b; namespace = true), :c; namespace = true)`.
"""
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
            # we need to create a closure to escape explicit return in `body`.
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

"""
    $(TYPEDSIGNATURES)

Mark a system constructor function as building a component. For example,

```julia
@component function AddOne(; name)
    @variables in(t) out(t)
    eqs = [out ~ in + 1]
    return System(eqs, t, [in, out], []; name)
end
```

ModelingToolkit systems are either components or connectors. Components define dynamics of
the model. Connectors are used to connect components together. See the
[Model building reference](@ref model_building_api) section of the documentation for more
information.

See also: [`@connector`](@ref).
"""
macro component(expr)
    esc(component_post_processing(expr, false))
end

"""
    $(TYPEDSIGNATURES)

Macro shorthand for building and compiling a system in one step.

```julia
@mtkcompile sys = Constructor(args...; kwargs....)
```

Is shorthand for

```julia
@named sys = Constructor(args...; kwargs...)
sys = mtkcompile(sys)
```
"""
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

"""
    debug_system(sys::AbstractSystem; functions = [log, sqrt, (^), /, inv, asin, acos], error_nonfinite = true)

Wrap `functions` in `sys` so any error thrown in them shows helpful symbolic-numeric
information about its input. If `error_nonfinite`, functions that output nonfinite
values (like `Inf` or `NaN`) also display errors, even though the raw function itself
does not throw an exception (like `1/0`). For example:

```julia-repl
julia> sys = debug_system(complete(sys))

julia> prob = ODEProblem(sys, [0.0, 2.0], (0.0, 1.0))

julia> prob.f(prob.u0, prob.p, 0.0)
ERROR: Function /(1, sin(P(t))) output non-finite value Inf with input
  1 => 1
  sin(P(t)) => 0.0
```

Additionally, all assertions in the system are optionally logged when they fail.
A new parameter is also added to the system which controls whether the message associated
with each assertion will be logged when the assertion fails. This parameter defaults to
`true` and can be toggled by symbolic indexing with
`ModelingToolkit.ASSERTION_LOG_VARIABLE`. For example,
`prob.ps[ModelingToolkit.ASSERTION_LOG_VARIABLE] = false` will disable logging.
"""
function debug_system(
        sys::AbstractSystem; functions = [log, sqrt, (^), /, inv, asin, acos], kw...)
    if !(functions isa Set)
        functions = Set(functions) # more efficient "in" lookup
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
"""
    hierarchy(sys::AbstractSystem; describe = false, bold = describe, kwargs...)

Print a tree of a system's hierarchy of subsystems.

# Keyword arguments

- `describe`: Whether to also print the description of each subsystem, if present.
- `bold`: Whether to print the name of the system in **bold** font.
"""
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

###
### Inheritance & composition
###
"""
$(TYPEDSIGNATURES)

Extend `basesys` with `sys`. This can be thought of as the `merge` operation on systems.
Values in `sys` take priority over duplicates in `basesys` (for example, defaults).

By default, the resulting system inherits `sys`'s name and description.

The `&` operator can also be used for this purpose. `sys & basesys` is equivalent to
`extend(sys, basesys)`.

See also [`compose`](@ref).
"""
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

    # collect fields common to all system types
    eqs = union(get_eqs(basesys), get_eqs(sys))
    sts = union(get_unknowns(basesys), get_unknowns(sys))
    ps = union(get_ps(basesys), get_ps(sys))
    dep_ps = union(get_parameter_dependencies(basesys), get_parameter_dependencies(sys))
    obs = union(get_observed(basesys), get_observed(sys))
    cevs = union(get_continuous_events(basesys), get_continuous_events(sys))
    devs = union(get_discrete_events(basesys), get_discrete_events(sys))
    defs = merge(get_defaults(basesys), get_defaults(sys)) # prefer `sys`
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

    # collect fields specific to some system types
    ieqs = union(get_initialization_eqs(basesys), get_initialization_eqs(sys))
    guesses = merge(get_guesses(basesys), get_guesses(sys)) # prefer `sys`
    kwargs = merge(kwargs, (initialization_eqs = ieqs, guesses = guesses))

    if has_assertions(basesys)
        kwargs = merge(
            kwargs, (; assertions = merge(get_assertions(basesys), get_assertions(sys))))
    end

    newsys = T(args...; kwargs...)
    @set! newsys.parameter_dependencies = dep_ps

    return newsys
end

"""
    $(TYPEDSIGNATURES)

Extend `sys` with all systems in `basesys` in order.
"""
function extend(sys, basesys::Vector{T}) where {T <: AbstractSystem}
    foldl(extend, basesys, init = sys)
end

"""
    $(TYPEDSIGNATURES)

Syntactic sugar for `extend(sys, basesys)`.

See also: [`extend`](@ref).
"""
function Base.:(&)(sys::AbstractSystem, basesys::AbstractSystem; kwargs...)
    extend(sys, basesys; kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Syntactic sugar for `extend(sys, basesys)`.
"""
function Base.:(&)(
        sys::AbstractSystem, basesys::Vector{T}; kwargs...) where {T <: AbstractSystem}
    extend(sys, basesys; kwargs...)
end

"""
$(SIGNATURES)

Compose multiple systems together. This adds all of `systems` as subsystems of `sys`.
The resulting system inherits the name of `sys` by default.

The `∘` operator can also be used for this purpose. `sys ∘ basesys` is equivalent to
`compose(sys, basesys)`.

See also [`extend`](@ref).
"""
function compose(sys::AbstractSystem, systems::AbstractArray; name = nameof(sys))
    nsys = length(systems)
    nsys == 0 && return sys
    @set! sys.name = name
    @set! sys.systems = [get_systems(sys); systems]
    if has_is_dde(sys)
        @set! sys.is_dde = _check_if_dde(equations(sys), get_iv(sys), get_systems(sys))
    end
    newunknowns = OrderedSet()
    newparams = OrderedSet()
    iv = has_iv(sys) ? get_iv(sys) : nothing
    for ssys in systems
        collect_scoped_vars!(newunknowns, newparams, ssys, iv)
    end
    @set! sys.unknowns = unique!(vcat(get_unknowns(sys), collect(newunknowns)))
    @set! sys.ps = unique!(vcat(get_ps(sys), collect(newparams)))

    newinputs = map(systems) do sys
        map(x -> namespace_expr(x, sys), inputs(sys))
    end
    newoutputs = map(systems) do sys
        map(x -> namespace_expr(x, sys), outputs(sys))
    end

    newinputs = reduce(vcat, newinputs)
    newoutputs = reduce(vcat, newoutputs)

    @set! sys.inputs = OrderedSet(vcat(inputs(sys), newinputs))
    @set! sys.outputs = OrderedSet(vcat(outputs(sys), newoutputs))
    return sys
end
"""
    $(TYPEDSIGNATURES)

Syntactic sugar for adding all systems in `syss` as the subsystems of `first(syss)`.
"""
function compose(syss...; name = nameof(first(syss)))
    compose(first(syss), collect(syss[2:end]); name = name)
end

"""
    $(TYPEDSIGNATURES)

Syntactic sugar for `compose(sys1, sys2)`.

See also: [`compose`](@ref).
"""
Base.:(∘)(sys1::AbstractSystem, sys2::AbstractSystem) = compose(sys1, sys2)

function UnPack.unpack(sys::ModelingToolkit.AbstractSystem, ::Val{p}) where {p}
    getproperty(sys, p; namespace = false)
end

"""
    missing_variable_defaults(sys::AbstractSystem, default = 0.0; subset = unknowns(sys))

Returns a `Vector{Pair}` of variables set to `default` which are missing from `get_defaults(sys)`.  The `default` argument can be a single value or vector to set the missing defaults respectively.
"""
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
        # post-walk to avoid infinite recursion
        @set! sys.systems = map(Base.Fix2(substitute, dict), systems)
        something(get(rules, nameof(sys), nothing), sys)
    elseif sys isa System
        rules = todict(map(r -> Symbolics.unwrap(r[1]) => Symbolics.unwrap(r[2]),
            collect(rules)))
        newsys = @set sys.eqs = fast_substitute(get_eqs(sys), rules)
        @set! newsys.unknowns = map(get_unknowns(sys)) do var
            get(rules, var, var)
        end
        @set! newsys.ps = map(get_ps(sys)) do var
            get(rules, var, var)
        end
        @set! newsys.parameter_dependencies = fast_substitute(
            get_parameter_dependencies(sys), rules)
        @set! newsys.defaults = Dict(fast_substitute(k, rules) => fast_substitute(v, rules)
        for (k, v) in get_defaults(sys))
        @set! newsys.guesses = Dict(fast_substitute(k, rules) => fast_substitute(v, rules)
        for (k, v) in get_guesses(sys))
        @set! newsys.noise_eqs = fast_substitute(get_noise_eqs(sys), rules)
        @set! newsys.costs = Vector{Union{Real, BasicSymbolic}}(fast_substitute(
            get_costs(sys), rules))
        @set! newsys.observed = fast_substitute(get_observed(sys), rules)
        @set! newsys.initialization_eqs = fast_substitute(
            get_initialization_eqs(sys), rules)
        @set! newsys.constraints = fast_substitute(get_constraints(sys), rules)
        @set! newsys.systems = map(s -> substitute(s, rules), get_systems(sys))
    else
        error("substituting symbols is not supported for $(typeof(sys))")
    end
end

"""
    $(TYPEDSIGNATURES)

Find equations of `sys` involving only parameters and separate them out into the
`parameter_dependencies` field. Relative ordering of equations is maintained.
Parameter-only equations are validated to be explicit and sorted topologically. All such
explicitly determined parameters are removed from the parameters of `sys`. Return the new
system.
"""
function process_parameter_equations(sys::AbstractSystem)
    if !isempty(get_systems(sys))
        throw(ArgumentError("Expected flattened system"))
    end
    varsbuf = Set()
    pareq_idxs = Int[]
    eqs = equations(sys)
    for (i, eq) in enumerate(eqs)
        empty!(varsbuf)
        vars!(varsbuf, eq; op = Union{Differential, Initial, Pre})
        # singular equations
        isempty(varsbuf) && continue
        if all(varsbuf) do sym
            is_parameter(sys, sym) ||
                symbolic_type(sym) == ArraySymbolic() &&
                is_sized_array_symbolic(sym) &&
                all(Base.Fix1(is_parameter, sys), collect(sym)) ||
                iscall(sym) &&
                operation(sym) === getindex && is_parameter(sys, arguments(sym)[1])
        end
            # Everything in `varsbuf` is a parameter, so this is a cheap `is_parameter`
            # check.
            if !(eq.lhs in varsbuf)
                throw(ArgumentError("""
                LHS of parameter dependency equation must be a single parameter. Found \
                $(eq.lhs).
                """))
            end
            push!(pareq_idxs, i)
        end
    end

    pareqs = [get_parameter_dependencies(sys); eqs[pareq_idxs]]
    explicitpars = [eq.lhs for eq in pareqs]
    pareqs = topsort_equations(pareqs, explicitpars)

    eqs = eqs[setdiff(eachindex(eqs), pareq_idxs)]

    @set! sys.eqs = eqs
    @set! sys.parameter_dependencies = pareqs
    @set! sys.ps = setdiff(get_ps(sys), explicitpars)
    return sys
end

"""
    dump_parameters(sys::AbstractSystem)

Return an array of `NamedTuple`s containing the metadata associated with each parameter in
`sys`. Also includes the default value of the parameter, if provided.

```@example
using ModelingToolkit
using DynamicQuantities
using ModelingToolkit: t, D

@parameters p = 1.0, [description = "My parameter", tunable = false] q = 2.0, [description = "Other parameter"]
@variables x(t) = 3.0 [unit = u"m"]
@named sys = System(Equation[], t, [x], [p, q])
ModelingToolkit.dump_parameters(sys)
```

See also: [`ModelingToolkit.dump_variable_metadata`](@ref), [`ModelingToolkit.dump_unknowns`](@ref)
"""
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

"""
    dump_unknowns(sys::AbstractSystem)

Return an array of `NamedTuple`s containing the metadata associated with each unknown in
`sys`. Also includes the default value of the unknown, if provided.

```@example
using ModelingToolkit
using DynamicQuantities
using ModelingToolkit: t, D

@parameters p = 1.0, [description = "My parameter", tunable = false] q = 2.0, [description = "Other parameter"]
@variables x(t) = 3.0 [unit = u"m"]
@named sys = System(Equation[], t, [x], [p, q])
ModelingToolkit.dump_unknowns(sys)
```

See also: [`ModelingToolkit.dump_variable_metadata`](@ref), [`ModelingToolkit.dump_parameters`](@ref)
"""
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

"""
    $(TYPEDSIGNATURES)

Return the variable in `sys` referred to by its string representation `str`.
Roughly supports the following CFG:

```
varname                  = "D(" varname ")" | "Differential(" iv ")(" varname ")" | arrvar | maybe_dummy_var
arrvar                   = maybe_dummy_var "[idxs...]"
idxs                     = int | int "," idxs
maybe_dummy_var          = namespacedvar | namespacedvar "(" iv ")" |
                           namespacedvar "(" iv ")" "ˍ" ts | namespacedvar "ˍ" ts |
                           namespacedvar "ˍ" ts "(" iv ")"
ts                       = iv | iv ts
namespacedvar            = ident "₊" namespacedvar | ident "." namespacedvar | ident
```

Where `iv` is the independent variable, `int` is an integer and `ident` is an identifier.
"""
function parse_variable(sys::AbstractSystem, str::AbstractString)
    iv = has_iv(sys) ? string(getname(get_iv(sys))) : nothing

    # I'd write a regex to validate `str`, but https://xkcd.com/1171/
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

### Functions for accessing algebraic/differential equations in systems ###

"""
    is_diff_equation(eq)

Return `true` if the input is a differential equation, i.e. an equation that contains a
differential term.

Example:
```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
@parameters p d
@variables X(t)
eq1 = D(X) ~ p - d*X
eq2 = 0 ~ p - d*X

is_diff_equation(eq1) # true
is_diff_equation(eq2) # false
```
"""
function is_diff_equation(eq)
    (eq isa Equation) || (return false)
    isdefined(eq, :lhs) && recursive_hasoperator(Union{Differential, Shift}, eq.lhs) &&
        (return true)
    isdefined(eq, :rhs) && recursive_hasoperator(Union{Differential, Shift}, eq.rhs) &&
        (return true)
    return false
end

"""
    is_alg_equation(eq)

Return `true` if the input is an algebraic equation, i.e. an equation that does not contain
any differentials.

Example:
```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
@parameters p d
@variables X(t)
eq1 = D(X) ~ p - d*X
eq2 = 0 ~ p - d*X

is_alg_equation(eq1) # false
is_alg_equation(eq2) # true
```
"""
function is_alg_equation(eq)
    return (eq isa Equation) && !is_diff_equation(eq)
end

"""
    alg_equations(sys::AbstractSystem)

For a system, returns a vector of all its algebraic equations (i.e. that does not contain any
differentials).

Example:
```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
@parameters p d
@variables X(t)
eq1 = D(X) ~ p - d*X
eq2 = 0 ~ p - d*X
@named osys = System([eq1, eq2], t)

alg_equations(osys) # returns `[0 ~ p - d*X(t)]`.
"""
alg_equations(sys::AbstractSystem) = filter(is_alg_equation, equations(sys))

"""
    diff_equations(sys::AbstractSystem)

For a system, returns a vector of all its differential equations (i.e. that does contain a differential).

Example:
```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
@parameters p d
@variables X(t)
eq1 = D(X) ~ p - d*X
eq2 = 0 ~ p - d*X
@named osys = System([eq1, eq2], t)

diff_equations(osys) # returns `[Differential(t)(X(t)) ~ p - d*X(t)]`.
"""
diff_equations(sys::AbstractSystem) = filter(is_diff_equation, equations(sys))

"""
    has_alg_equations(sys::AbstractSystem)

For a system, returns true if it contain at least one algebraic equation (i.e. that does not contain any
differentials).

Example:
```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
@parameters p d
@variables X(t)
eq1 = D(X) ~ p - d*X
eq2 = 0 ~ p - d*X
@named osys1 = System([eq1], t)
@named osys2 = System([eq2], t)

has_alg_equations(osys1) # returns `false`.
has_alg_equations(osys2) # returns `true`.
```
"""
has_alg_equations(sys::AbstractSystem) = any(is_alg_equation, equations(sys))

"""
    has_diff_equations(sys::AbstractSystem)

For a system, returns true if it contain at least one differential equation (i.e. that contain a differential).

Example:
```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
@parameters p d
@variables X(t)
eq1 = D(X) ~ p - d*X
eq2 = 0 ~ p - d*X
@named osys1 = System([eq1], t)
@named osys2 = System([eq2], t)

has_diff_equations(osys1) # returns `true`.
has_diff_equations(osys2) # returns `false`.
```
"""
has_diff_equations(sys::AbstractSystem) = any(is_diff_equation, equations(sys))

"""
    get_alg_eqs(sys::AbstractSystem)

For a system, returns a vector of all algebraic equations (i.e. that does not contain any
differentials) in its *top-level system*.

Example:
```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
@parameters p d
@variables X(t)
eq1 = D(X) ~ p - d*X
eq2 = 0 ~ p - d*X
@named osys1 = ([eq1], t)
@named osys2 = ([eq2], t)
osys12 = compose(sys1, [osys2])
osys21 = compose(osys2, [osys1])

get_alg_eqs(osys12) # returns `Equation[]`.
get_alg_eqs(osys21) # returns `[0 ~ p - d*X(t)]`.
```
"""
get_alg_eqs(sys::AbstractSystem) = filter(is_alg_equation, get_eqs(sys))

"""
    get_diff_eqs(sys::AbstractSystem)

For a system, returns a vector of all differential equations (i.e. that does contain a differential)
in its *top-level system*.

Example:
```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
@parameters p d
@variables X(t)
eq1 = D(X) ~ p - d*X
eq2 = 0 ~ p - d*X
@named osys1 = tem([eq1], t)
@named osys2 = tem([eq2], t)
osys12 = compose(osys1, [osys2])
osys21 = compose(osys2, [osys1])

get_diff_eqs(osys12) # returns `[Differential(t)(X(t)) ~ p - d*X(t)]`.
get_diff_eqs(osys21) # returns `Equation[]``.
```
"""
get_diff_eqs(sys::AbstractSystem) = filter(is_diff_equation, get_eqs(sys))

"""
    has_alg_eqs(sys::AbstractSystem)

For a system, returns true if it contain at least one algebraic equation (i.e. that does not contain any
differentials) in its *top-level system*.

Example:
```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
@parameters p d
@variables X(t)
eq1 = D(X) ~ p - d*X
eq2 = 0 ~ p - d*X
@named osys1 = System([eq1], t)
@named osys2 = System([eq2], t)
osys12 = compose(osys1, [osys2])
osys21 = compose(osys2, [osys1])

has_alg_eqs(osys12) # returns `false`.
has_alg_eqs(osys21) # returns `true`.
```
"""
has_alg_eqs(sys::AbstractSystem) = any(is_alg_equation, get_eqs(sys))

"""
    has_diff_eqs(sys::AbstractSystem)

Return `true` if a system contains at least one differential equation (i.e. an equation with a
differential term). Note that this does not consider subsystems, and only takes into account
equations in the top-level system.

Example:
```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
@parameters p d
@variables X(t)
eq1 = D(X) ~ p - d*X
eq2 = 0 ~ p - d*X
@named osys1 = tem([eq1], t)
@named osys2 = tem([eq2], t)
osys12 = compose(osys1, [osys2])
osys21 = compose(osys2, [osys1])

has_diff_eqs(osys12) # returns `true`.
has_diff_eqs(osys21) # returns `false`.
```
"""
has_diff_eqs(sys::AbstractSystem) = any(is_diff_equation, get_eqs(sys))

"""
    $(TYPEDSIGNATURES)

Validate the rules for replacement of subcomponents as defined in `substitute_component`.
"""
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

"""
    $(TYPEDSIGNATURES)

Chain `getproperty` calls on `root` in the order given in `hierarchy`.

# Keyword Arguments

- `skip_namespace_first`: Whether to avoid namespacing in the first `getproperty` call.
"""
function recursive_getproperty(
        root::AbstractSystem, hierarchy::Vector{Symbol}; skip_namespace_first = true)
    cur = root
    for (i, name) in enumerate(hierarchy)
        cur = getproperty(cur, name; namespace = i > 1 || !skip_namespace_first)
    end
    return unwrap(cur)
end

"""
    $(TYPEDSIGNATURES)

Recursively descend through `sys`, finding all connection equations and re-creating them
using the names of the involved variables/systems and finding the required variables/
systems in the hierarchy.
"""
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

"""
    $(TYPEDSIGNATURES)

Given a hierarchical system `sys` and a rule `lhs => rhs`, replace the subsystem `lhs` in
`sys` by `rhs`. The `lhs` must be the namespaced version of a subsystem of `sys` (e.g.
obtained via `sys.inner.component`). The `rhs` must be valid as per the following
conditions:

1. `rhs` must not be namespaced.
2. The name of `rhs` must be the same as the unnamespaced name of `lhs`.
3. Neither one of `lhs` or `rhs` can be marked as complete.
4. Both `lhs` and `rhs` must share the same independent variable.
5. `rhs` must contain at least all of the unknowns and parameters present in
   `lhs`.
6. Corresponding unknowns in `rhs` must share the same connection and causality
   (input/output) metadata as their counterparts in `lhs`.
7. For each subsystem of `lhs`, there must be an identically named subsystem of `rhs`.
   These two corresponding subsystems must satisfy conditions 3, 4, 5, 6, 7. If the
   subsystem of `lhs` is a connector, the corresponding subsystem of `rhs` must also
   be a connector of the same type.

`sys` also cannot be marked as complete.
"""
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
