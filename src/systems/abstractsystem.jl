const SYSTEM_COUNT = Threads.Atomic{UInt}(0)

get_component_type(x::AbstractSystem) = get_gui_metadata(x).type
struct GUIMetadata
    type::GlobalRef
    layout::Any
end

GUIMetadata(type) = GUIMetadata(type, nothing)

"""
```julia
calculate_tgrad(sys::AbstractTimeDependentSystem)
```

Calculate the time gradient of a system.

Returns a vector of [`Num`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_tgrad end

"""
```julia
calculate_gradient(sys::AbstractSystem)
```

Calculate the gradient of a scalar system.

Returns a vector of [`Num`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_gradient end

"""
```julia
calculate_jacobian(sys::AbstractSystem)
```

Calculate the Jacobian matrix of a system.

Returns a matrix of [`Num`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_jacobian end

"""
```julia
calculate_control_jacobian(sys::AbstractSystem)
```

Calculate the Jacobian matrix of a system with respect to the system's controls.

Returns a matrix of [`Num`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_control_jacobian end

"""
```julia
calculate_factorized_W(sys::AbstractSystem)
```

Calculate the factorized W-matrix of a system.

Returns a matrix of [`Num`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_factorized_W end

"""
```julia
calculate_hessian(sys::AbstractSystem)
```

Calculate the hessian matrix of a scalar system.

Returns a matrix of [`Num`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_hessian end

"""
```julia
generate_tgrad(sys::AbstractTimeDependentSystem, dvs = unknowns(sys), ps = parameters(sys),
               expression = Val{true}; kwargs...)
```

Generates a function for the time gradient of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_tgrad end

"""
```julia
generate_gradient(sys::AbstractSystem, dvs = unknowns(sys), ps = parameters(sys),
                  expression = Val{true}; kwargs...)
```

Generates a function for the gradient of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_gradient end

"""
```julia
generate_jacobian(sys::AbstractSystem, dvs = unknowns(sys), ps = parameters(sys),
                  expression = Val{true}; sparse = false, kwargs...)
```

Generates a function for the Jacobian matrix of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_jacobian end

"""
```julia
generate_factorized_W(sys::AbstractSystem, dvs = unknowns(sys), ps = parameters(sys),
                      expression = Val{true}; sparse = false, kwargs...)
```

Generates a function for the factorized W matrix of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_factorized_W end

"""
```julia
generate_hessian(sys::AbstractSystem, dvs = unknowns(sys), ps = parameters(sys),
                 expression = Val{true}; sparse = false, kwargs...)
```

Generates a function for the hessian matrix of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_hessian end

"""
```julia
generate_function(sys::AbstractSystem, dvs = unknowns(sys), ps = parameters(sys),
                  expression = Val{true}; kwargs...)
```

Generate a function to evaluate the system's equations.
"""
function generate_function end

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
[`structural_simplify`](@ref) or [`@mtkbuild`](@ref), `p` is expected to be an `MTKParameters`
object.
"""
function generate_custom_function(sys::AbstractSystem, exprs, dvs = unknowns(sys),
        ps = parameters(sys);
        expression = Val{true}, eval_expression = false, eval_module = @__MODULE__,
        cachesyms::Tuple = (), kwargs...)
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system.")
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

function wrap_parameter_dependencies(sys::AbstractSystem, isscalar)
    wrap_assignments(isscalar, [eq.lhs ← eq.rhs for eq in parameter_dependencies(sys)])
end

"""
    $(TYPEDSIGNATURES)

Add the necessary assignment statements to allow use of unscalarized array variables
in the generated code. `expr` is the expression returned by the function. `dvs` and
`ps` are the unknowns and parameters of the system `sys` to use in the generated code.
`inputs` can be specified as an array of symbolics if the generated function has inputs.
If `history == true`, the generated function accepts a history function. `cachesyms` are
extra variables (arrays of variables) stored in the cache array(s) of the parameter
object. `extra_args` are extra arguments appended to the end of the argument list.

The function is assumed to have the signature `f(du, u, h, x, p, cache_syms..., t, extra_args...)`
Where:
- `du` is the optional buffer to write to for in-place functions.
- `u` is the list of unknowns. This argument is not present if `dvs === nothing`.
- `h` is the optional history function, present if `history == true`.
- `x` is the array of inputs, present only if `inputs !== nothing`. Values are assumed
  to be in the order of variables passed to `inputs`.
- `p` is the parameter object.
- `cache_syms` are the cache variables. These are part of the splatted parameter object.
- `t` is time, present only if the system is time dependent.
- `extra_args` are the extra arguments passed to the function, present only if
  `extra_args` is non-empty.
"""
function wrap_array_vars(
        sys::AbstractSystem, exprs; dvs = unknowns(sys), ps = parameters(sys),
        inputs = nothing, history = false, cachesyms::Tuple = (), extra_args::Tuple = ())
    isscalar = !(exprs isa AbstractArray)
    var_to_arridxs = Dict()

    if dvs === nothing
        uind = 0
    else
        uind = 1
        for (i, x) in enumerate(dvs)
            iscall(x) && operation(x) == getindex || continue
            arg = arguments(x)[1]
            inds = get!(() -> [], var_to_arridxs, arg)
            push!(inds, (uind, i))
        end
    end
    p_start = uind + 1 + history
    rps = (reorder_parameters(sys, ps)..., cachesyms...)
    if inputs !== nothing
        rps = (inputs, rps...)
    end
    if has_iv(sys)
        rps = (rps..., get_iv(sys))
    end
    rps = (rps..., extra_args...)
    for sym in reduce(vcat, rps; init = [])
        iscall(sym) && operation(sym) == getindex || continue
        arg = arguments(sym)[1]

        bufferidx = findfirst(buf -> any(isequal(sym), buf), rps)
        idxinbuffer = findfirst(isequal(sym), rps[bufferidx])
        inds = get!(() -> [], var_to_arridxs, arg)
        push!(inds, (p_start + bufferidx - 1, idxinbuffer))
    end

    viewsyms = Dict()
    splitsyms = Dict()
    for (arrsym, idxs) in var_to_arridxs
        length(idxs) == length(arrsym) || continue
        # allequal(first, idxs) is a 1.11 feature
        if allequal(Iterators.map(first, idxs))
            viewsyms[arrsym] = (first(first(idxs)), reshape(last.(idxs), size(arrsym)))
        else
            splitsyms[arrsym] = reshape(idxs, size(arrsym))
        end
    end
    if isscalar
        function (expr)
            Func(
                expr.args,
                [],
                Let(
                    vcat(
                        [sym ← :(view($(expr.args[i].name), $idxs))
                         for (sym, (i, idxs)) in viewsyms],
                        [sym ←
                         MakeArray([expr.args[bufi].elems[vali] for (bufi, vali) in idxs],
                             expr.args[idxs[1][1]]) for (sym, idxs) in splitsyms]
                    ),
                    expr.body,
                    false
                )
            )
        end
    else
        function (expr)
            Func(
                expr.args,
                [],
                Let(
                    vcat(
                        [sym ← :(view($(expr.args[i].name), $idxs))
                         for (sym, (i, idxs)) in viewsyms],
                        [sym ←
                         MakeArray([expr.args[bufi].elems[vali] for (bufi, vali) in idxs],
                             expr.args[idxs[1][1]]) for (sym, idxs) in splitsyms]
                    ),
                    expr.body,
                    false
                )
            )
        end,
        function (expr)
            Func(
                expr.args,
                [],
                Let(
                    vcat(
                        [sym ← :(view($(expr.args[i + 1].name), $idxs))
                         for (sym, (i, idxs)) in viewsyms],
                        [sym ← MakeArray(
                             [expr.args[bufi + 1].elems[vali] for (bufi, vali) in idxs],
                             expr.args[idxs[1][1] + 1]) for (sym, idxs) in splitsyms]
                    ),
                    expr.body,
                    false
                )
            )
        end
    end
end

const MTKPARAMETERS_ARG = Sym{Vector{Vector}}(:___mtkparameters___)

"""
    wrap_mtkparameters(sys::AbstractSystem, isscalar::Bool, p_start = 2, offset = Int(is_time_dependent(sys)))

Return function(s) to be passed to the `wrap_code` keyword of `build_function` which
allow the compiled function to be called as `f(u, p, t)` where `p isa MTKParameters`
instead of `f(u, p..., t)`. `isscalar` denotes whether the function expression being
wrapped is for a scalar value. `p_start` is the index of the argument containing
the first parameter vector in the out-of-place version of the function. For example,
if a history function (DDEs) was passed before `p`, then the function before wrapping
would have the signature `f(u, h, p..., t)` and hence `p_start` would need to be `3`.

`offset` is the number of arguments at the end of the argument list to ignore. Defaults
to 1 if the system is time-dependent (to ignore `t`) and 0 otherwise.

The returned function is `identity` if the system does not have an `IndexCache`.
"""
function wrap_mtkparameters(sys::AbstractSystem, isscalar::Bool, p_start = 2,
        offset = Int(is_time_dependent(sys)))
    if has_index_cache(sys) && get_index_cache(sys) !== nothing
        if isscalar
            function (expr)
                param_args = expr.args[p_start:(end - offset)]
                param_buffer_idxs = findall(x -> x isa DestructuredArgs, param_args)
                param_buffer_args = param_args[param_buffer_idxs]
                destructured_mtkparams = DestructuredArgs(
                    [x.name for x in param_buffer_args],
                    MTKPARAMETERS_ARG; inds = param_buffer_idxs)
                Func(
                    [
                        expr.args[begin:(p_start - 1)]...,
                        destructured_mtkparams,
                        expr.args[(end - offset + 1):end]...
                    ],
                    [],
                    Let(param_buffer_args, expr.body, false)
                )
            end
        else
            function (expr)
                param_args = expr.args[p_start:(end - offset)]
                param_buffer_idxs = findall(x -> x isa DestructuredArgs, param_args)
                param_buffer_args = param_args[param_buffer_idxs]
                destructured_mtkparams = DestructuredArgs(
                    [x.name for x in param_buffer_args],
                    MTKPARAMETERS_ARG; inds = param_buffer_idxs)
                Func(
                    [
                        expr.args[begin:(p_start - 1)]...,
                        destructured_mtkparams,
                        expr.args[(end - offset + 1):end]...
                    ],
                    [],
                    Let(param_buffer_args, expr.body, false)
                )
            end,
            function (expr)
                param_args = expr.args[(p_start + 1):(end - offset)]
                param_buffer_idxs = findall(x -> x isa DestructuredArgs, param_args)
                param_buffer_args = param_args[param_buffer_idxs]
                destructured_mtkparams = DestructuredArgs(
                    [x.name for x in param_buffer_args],
                    MTKPARAMETERS_ARG; inds = param_buffer_idxs)
                Func(
                    [
                        expr.args[begin:p_start]...,
                        destructured_mtkparams,
                        expr.args[(end - offset + 1):end]...
                    ],
                    [],
                    Let(param_buffer_args, expr.body, false)
                )
            end
        end
    else
        identity
    end
end

mutable struct Substitutions
    subs::Vector{Equation}
    deps::Vector{Vector{Int}}
    subed_eqs::Union{Nothing, Vector{Equation}}
end
Substitutions(subs, deps) = Substitutions(subs, deps, nothing)

Base.nameof(sys::AbstractSystem) = getfield(sys, :name)
description(sys::AbstractSystem) = has_description(sys) ? get_description(sys) : ""

#Deprecated
function independent_variable(sys::AbstractSystem)
    Base.depwarn(
        "`independent_variable` is deprecated. Use `get_iv` or `independent_variables` instead.",
        :independent_variable)
    isdefined(sys, :iv) ? getfield(sys, :iv) : nothing
end

function independent_variables(sys::AbstractTimeDependentSystem)
    return [getfield(sys, :iv)]
end

independent_variables(::AbstractTimeIndependentSystem) = []

function independent_variables(sys::AbstractMultivariateSystem)
    return getfield(sys, :ivs)
end

"""
$(TYPEDSIGNATURES)

Get the independent variable(s) of the system `sys`.

See also [`@independent_variables`](@ref) and [`ModelingToolkit.get_iv`](@ref).
"""
function independent_variables(sys::AbstractSystem)
    @warn "Please declare ($(typeof(sys))) as a subtype of `AbstractTimeDependentSystem`, `AbstractTimeIndependentSystem` or `AbstractMultivariateSystem`."
    if isdefined(sys, :iv)
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

SymbolicIndexingInterface.variable_symbols(sys::AbstractMultivariateSystem) = sys.dvs

function SymbolicIndexingInterface.variable_symbols(sys::AbstractSystem)
    return solved_unknowns(sys)
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
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        rawobs = build_explicit_observed_function(
            sys, sym; param_only = true, return_inplace = true)
        if rawobs isa Tuple
            if is_time_dependent(sys)
                obsfn = let oop = rawobs[1], iip = rawobs[2]
                    f1a(p, t) = oop(p, t)
                    f1a(out, p, t) = iip(out, p, t)
                end
            else
                obsfn = let oop = rawobs[1], iip = rawobs[2]
                    f1b(p) = oop(p)
                    f1b(out, p) = iip(out, p)
                end
            end
        else
            obsfn = rawobs
        end
    else
        obsfn = build_explicit_observed_function(sys, sym; param_only = true)
    end
    return obsfn
end

function has_observed_with_lhs(sys, sym)
    has_observed(sys) || return false
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        return haskey(ic.observed_syms_to_timeseries, sym)
    else
        return any(isequal(sym), [eq.lhs for eq in observed(sys)])
    end
end

function has_parameter_dependency_with_lhs(sys, sym)
    has_parameter_dependencies(sys) || return false
    if has_index_cache(sys) && (ic = get_index_cache(sys)) !== nothing
        return haskey(ic.dependent_pars_to_timeseries, unwrap(sym))
    else
        return any(isequal(sym), [eq.lhs for eq in parameter_dependencies(sys)])
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
                   is_variable(sys, operation(s)(get_iv(sys)))
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
           any(isequal(sym), [getname(eq.lhs) for eq in observed(sys)])
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
    return parameters(sys)
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
        sys::AbstractSystem, sym; eval_expression = false, eval_module = @__MODULE__, checkbounds = true)
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
    _fn = build_explicit_observed_function(
        sys, sym; eval_expression, eval_module, checkbounds)

    if is_time_dependent(sys)
        return _fn
    else
        return let _fn = _fn
            fn2(u, p) = _fn(u, p)
            fn2(::Nothing, p) = _fn([], p)
            fn2
        end
    end
end

function SymbolicIndexingInterface.default_values(sys::AbstractSystem)
    return merge(
        Dict(eq.lhs => eq.rhs for eq in observed(sys)),
        defaults(sys)
    )
end

SymbolicIndexingInterface.is_time_dependent(::AbstractTimeDependentSystem) = true
SymbolicIndexingInterface.is_time_dependent(::AbstractTimeIndependentSystem) = false

SymbolicIndexingInterface.is_markovian(sys::AbstractSystem) = !is_dde(sys)

SymbolicIndexingInterface.constant_structure(::AbstractSystem) = true

function SymbolicIndexingInterface.all_variable_symbols(sys::AbstractSystem)
    syms = variable_symbols(sys)
    obs = getproperty.(observed(sys), :lhs)
    return isempty(obs) ? syms : vcat(syms, obs)
end

function SymbolicIndexingInterface.all_symbols(sys::AbstractSystem)
    syms = all_variable_symbols(sys)
    for other in (full_parameters(sys), independent_variable_symbols(sys))
        isempty(other) || (syms = vcat(syms, other))
    end
    return syms
end

iscomplete(sys::AbstractSystem) = isdefined(sys, :complete) && getfield(sys, :complete)

"""
$(TYPEDSIGNATURES)

Mark a system as scheduled. It is only intended in compiler internals. A system
is scheduled after tearing based simplifications where equations are converted
into assignments.
"""
function schedule(sys::AbstractSystem)
    has_schedule(sys) ? sys : (@set! sys.isscheduled = true)
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
$(TYPEDSIGNATURES)

Mark a system as completed. A completed system is a system which is done being
defined/modified and is ready for structural analysis or other transformations.
This allows for analyses and optimizations to be performed which require knowing
the global structure of the system.

One property to note is that if a system is complete, the system will no longer
namespace its subsystems or variables, i.e. `isequal(complete(sys).v.i, v.i)`.
"""
function complete(sys::AbstractSystem; split = true, flatten = true)
    newunknowns = OrderedSet()
    newparams = OrderedSet()
    iv = has_iv(sys) ? get_iv(sys) : nothing
    collect_scoped_vars!(newunknowns, newparams, sys, iv; depth = -1)
    # don't update unknowns to not disturb `structural_simplify` order
    # `GlobalScope`d unknowns will be picked up and added there
    @set! sys.ps = unique!(vcat(get_ps(sys), collect(newparams)))

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
    end
    if split && has_index_cache(sys)
        @set! sys.index_cache = IndexCache(sys)
        all_ps = get_ps(sys)
        if !isempty(all_ps)
            # reorder parameters by portions
            ps_split = reorder_parameters(sys, all_ps)
            # if there are no tunables, vcat them
            if isempty(get_index_cache(sys).tunable_idx)
                ordered_ps = reduce(vcat, ps_split)
            else
                # if there are tunables, they will all be in `ps_split[1]`
                # and the arrays will have been scalarized
                ordered_ps = eltype(all_ps)[]
                i = 1
                # go through all the tunables
                while i <= length(ps_split[1])
                    sym = ps_split[1][i]
                    # if the sym is not a scalarized array symbolic OR it was already scalarized,
                    # just push it as-is
                    if !iscall(sym) || operation(sym) != getindex ||
                       any(isequal(sym), all_ps)
                        push!(ordered_ps, sym)
                        i += 1
                        continue
                    end
                    # the next `length(sym)` symbols should be scalarized versions of the same
                    # array symbolic
                    if !allequal(first(arguments(x))
                    for x in view(ps_split[1], i:(i + length(sym) - 1)))
                        error("This should not be possible. Please open an issue in ModelingToolkit.jl with an MWE and stacktrace.")
                    end
                    arrsym = first(arguments(sym))
                    push!(ordered_ps, arrsym)
                    i += length(arrsym)
                end
                ordered_ps = vcat(
                    ordered_ps, reduce(vcat, ps_split[2:end]; init = eltype(ordered_ps)[]))
            end
            @set! sys.ps = ordered_ps
        end
    elseif has_index_cache(sys)
        @set! sys.index_cache = nothing
    end
    if isdefined(sys, :initializesystem) && get_initializesystem(sys) !== nothing
        @set! sys.initializesystem = complete(get_initializesystem(sys); split)
    end
    isdefined(sys, :complete) ? (@set! sys.complete = true) : sys
end

for prop in [:eqs
             :tag
             :noiseeqs
             :iv
             :unknowns
             :ps
             :tspan
             :name
             :description
             :var_to_name
             :ctrls
             :defaults
             :guesses
             :observed
             :tgrad
             :jac
             :ctrl_jac
             :Wfact
             :Wfact_t
             :systems
             :structure
             :op
             :constraints
             :controls
             :loss
             :bcs
             :domain
             :ivs
             :dvs
             :connector_type
             :connections
             :preface
             :torn_matching
             :initializesystem
             :initialization_eqs
             :schedule
             :tearing_state
             :substitutions
             :metadata
             :gui_metadata
             :discrete_subsystems
             :parameter_dependencies
             :assertions
             :solved_unknowns
             :split_idxs
             :parent
             :is_dde
             :tstops
             :index_cache
             :is_scalar_noise
             :isscheduled]
    fname_get = Symbol(:get_, prop)
    fname_has = Symbol(:has_, prop)
    @eval begin
        """
        $(TYPEDSIGNATURES)

        Get the internal field `$($(QuoteNode(prop)))` of a system `sys`.
        It only includes `$($(QuoteNode(prop)))` local to `sys`; not those of its subsystems,
        like `unknowns(sys)`, `parameters(sys)` and `equations(sys)` does.
        This is equivalent to, but preferred over `sys.$($(QuoteNode(prop)))`.

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

const EMPTY_TGRAD = Vector{Num}(undef, 0)
const EMPTY_JAC = Matrix{Num}(undef, 0, 0)
function invalidate_cache!(sys::AbstractSystem)
    if has_tgrad(sys)
        get_tgrad(sys)[] = EMPTY_TGRAD
    end
    if has_jac(sys)
        get_jac(sys)[] = EMPTY_JAC
    end
    if has_ctrl_jac(sys)
        get_ctrl_jac(sys)[] = EMPTY_JAC
    end
    if has_Wfact(sys)
        get_Wfact(sys)[] = EMPTY_JAC
    end
    if has_Wfact_t(sys)
        get_Wfact_t(sys)[] = EMPTY_JAC
    end
    return sys
end

function Setfield.get(obj::AbstractSystem, ::Setfield.PropertyLens{field}) where {field}
    getfield(obj, field)
end
@generated function ConstructionBase.setproperties(obj::AbstractSystem, patch::NamedTuple)
    if issubset(fieldnames(patch), fieldnames(obj))
        args = map(fieldnames(obj)) do fn
            if fn in fieldnames(patch)
                :(patch.$fn)
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

rename(x, name) = @set x.name = name

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
            push!(names, getname(s))
        end
        has_observed(sys) && for s in get_observed(sys)
            push!(names, getname(s.lhs))
        end
        has_iv(sys) && push!(names, getname(get_iv(sys)))
        return names
    end
end

function Base.getproperty(sys::AbstractSystem, name::Symbol; namespace = !iscomplete(sys))
    if has_parent(sys) && (parent = get_parent(sys); parent !== nothing)
        return getproperty(parent, name; namespace)
    end
    wrap(getvar(sys, name; namespace = namespace))
end
function getvar(sys::AbstractSystem, name::Symbol; namespace = !iscomplete(sys))
    systems = get_systems(sys)
    if isdefined(sys, name)
        Base.depwarn(
            "`sys.name` like `sys.$name` is deprecated. Use getters like `get_$name` instead.",
            "sys.$name")
        return getfield(sys, name)
    elseif !isempty(systems)
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
            if eq.lhs isa AnalysisPoint && nameof(eq.rhs) == name
                return namespace ? renamespace(sys, eq.rhs) : eq.rhs
            end
        end
    end

    throw(ArgumentError("System $(nameof(sys)): variable $name does not exist"))
end

function Base.setproperty!(sys::AbstractSystem, prop::Symbol, val)
    # We use this weird syntax because `parameters` and `unknowns` calls are
    # potentially expensive.
    if (params = parameters(sys);
    idx = findfirst(s -> getname(s) == prop, params);
    idx !== nothing)
        get_defaults(sys)[params[idx]] = value(val)
    elseif (sts = unknowns(sys);
    idx = findfirst(s -> getname(s) == prop, sts);
    idx !== nothing)
        get_defaults(sys)[sts[idx]] = value(val)
    else
        setfield!(sys, prop, val)
    end
end

apply_to_variables(f::F, ex) where {F} = _apply_to_variables(f, ex)
apply_to_variables(f::F, ex::Num) where {F} = wrap(_apply_to_variables(f, unwrap(ex)))
function _apply_to_variables(f::F, ex) where {F}
    if isvariable(ex)
        return f(ex)
    end
    iscall(ex) || return ex
    maketerm(typeof(ex), _apply_to_variables(f, operation(ex)),
        map(Base.Fix1(_apply_to_variables, f), arguments(ex)),
        metadata(ex))
end

abstract type SymScope end

struct LocalScope <: SymScope end
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

struct ParentScope <: SymScope
    parent::SymScope
end
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

struct DelayParentScope <: SymScope
    parent::SymScope
    N::Int
end
function DelayParentScope(sym::Union{Num, Symbolic, Symbolics.Arr{Num}}, N)
    apply_to_variables(sym) do sym
        if iscall(sym) && operation(sym) == getindex
            args = arguments(sym)
            a1 = setmetadata(args[1], SymScope,
                DelayParentScope(getmetadata(value(args[1]), SymScope, LocalScope()), N))
            maketerm(typeof(sym), operation(sym), [a1, args[2:end]...],
                metadata(sym))
        else
            setmetadata(sym, SymScope,
                DelayParentScope(getmetadata(value(sym), SymScope, LocalScope()), N))
        end
    end
end
DelayParentScope(sym::Union{Num, Symbolic, Symbolics.Arr{Num}}) = DelayParentScope(sym, 1)

struct GlobalScope <: SymScope end
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
            elseif scope isa DelayParentScope
                if scope.N > 0
                    x = setmetadata(x, SymScope,
                        DelayParentScope(scope.parent, scope.N - 1))
                    rename(x, renamespace(getname(sys), getname(x)))::T
                else
                    #rename(x, renamespace(getname(sys), getname(x)))
                    setmetadata(x, SymScope, scope.parent)::T
                end
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
namespace_controls(sys::AbstractSystem) = controls(sys, controls(sys))

function namespace_defaults(sys)
    defs = defaults(sys)
    Dict((isparameter(k) ? parameters(sys, k) : unknowns(sys, k)) => namespace_expr(v, sys)
    for (k, v) in pairs(defs))
end

function namespace_guesses(sys)
    guess = guesses(sys)
    Dict(unknowns(sys, k) => namespace_expr(v, sys) for (k, v) in guess)
end

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

function namespace_equation(eq::Equation,
        sys,
        n = nameof(sys);
        ivs = independent_variables(sys))
    _lhs = namespace_expr(eq.lhs, sys, n; ivs)
    _rhs = namespace_expr(eq.rhs, sys, n; ivs)
    (_lhs ~ _rhs)::Equation
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
        O, sys, n = nameof(sys); ivs = independent_variables(sys))
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
$(TYPEDSIGNATURES)

Get the parameters of the system `sys` and its subsystems.

See also [`@parameters`](@ref) and [`ModelingToolkit.get_ps`](@ref).
"""
function parameters(sys::AbstractSystem)
    ps = get_ps(sys)
    if ps == SciMLBase.NullParameters()
        return []
    end
    if eltype(ps) <: Pair
        ps = first.(ps)
    end
    systems = get_systems(sys)
    unique(isempty(systems) ? ps : [ps; reduce(vcat, namespace_parameters.(systems))])
end

function dependent_parameters(sys::AbstractSystem)
    return map(eq -> eq.lhs, parameter_dependencies(sys))
end

"""
$(TYPEDSIGNATURES)
Get the parameter dependencies of the system `sys` and its subsystems.

See also [`defaults`](@ref) and [`ModelingToolkit.get_parameter_dependencies`](@ref).
"""
function parameter_dependencies(sys::AbstractSystem)
    if !has_parameter_dependencies(sys)
        return Equation[]
    end
    pdeps = get_parameter_dependencies(sys)
    systems = get_systems(sys)
    # put pdeps after those of subsystems to maintain topological sorted order
    namespaced_deps = mapreduce(
        s -> map(eq -> namespace_equation(eq, s), parameter_dependencies(s)), vcat,
        systems; init = Equation[])

    return vcat(namespaced_deps, pdeps)
end

function full_parameters(sys::AbstractSystem)
    vcat(parameters(sys), dependent_parameters(sys))
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

function controls(sys::AbstractSystem)
    ctrls = get_ctrls(sys)
    systems = get_systems(sys)
    isempty(systems) ? ctrls : [ctrls; reduce(vcat, namespace_controls.(systems))]
end

function observed(sys::AbstractSystem)
    obs = get_observed(sys)
    systems = get_systems(sys)
    [obs;
     reduce(vcat,
         (map(o -> namespace_equation(o, s), observed(s)) for s in systems),
         init = Equation[])]
end

Base.@deprecate default_u0(x) defaults(x) false
Base.@deprecate default_p(x) defaults(x) false

"""
$(TYPEDSIGNATURES)

Get the default values of the system sys and its subsystems.
If they are not explicitly provided, variables and parameters are initialized to these values.

See also [`initialization_equations`](@ref), [`parameter_dependencies`](@ref) and [`ModelingToolkit.get_defaults`](@ref).
"""
function defaults(sys::AbstractSystem)
    systems = get_systems(sys)
    defs = get_defaults(sys)
    # `mapfoldr` is really important!!! We should prefer the base model for
    # defaults, because people write:
    #
    # `compose(ODESystem(...; defaults=defs), ...)`
    #
    # Thus, right associativity is required and crucial for correctness.
    isempty(systems) ? defs : mapfoldr(namespace_defaults, merge, systems; init = defs)
end

function defaults_and_guesses(sys::AbstractSystem)
    merge(guesses(sys), defaults(sys))
end

unknowns(sys::Union{AbstractSystem, Nothing}, v) = renamespace(sys, v)
for vType in [Symbolics.Arr, Symbolics.Symbolic{<:AbstractArray}]
    @eval unknowns(sys::AbstractSystem, v::$vType) = renamespace(sys, v)
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
$(TYPEDSIGNATURES)

Get the initialization equations of the system `sys` and its subsystems.

See also [`guesses`](@ref), [`defaults`](@ref), [`parameter_dependencies`](@ref) and [`ModelingToolkit.get_initialization_eqs`](@ref).
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

function islinear(sys::AbstractSystem)
    rhs = [eq.rhs for eq in equations(sys)]

    all(islinear(r, unknowns(sys)) for r in rhs)
end

function isaffine(sys::AbstractSystem)
    rhs = [eq.rhs for eq in equations(sys)]

    all(isaffine(r, unknowns(sys)) for r in rhs)
end

function time_varying_as_func(x, sys::AbstractTimeDependentSystem)
    # if something is not x(t) (the current unknown)
    # but is `x(t-1)` or something like that, pass in `x` as a callable function rather
    # than pass in a value in place of x(t).
    #
    # This is done by just making `x` the argument of the function.
    if iscall(x) &&
       issym(operation(x)) &&
       !(length(arguments(x)) == 1 && isequal(arguments(x)[1], get_iv(sys)))
        return operation(x)
    end
    return x
end

"""
$(SIGNATURES)

Return a list of actual unknowns needed to be solved by solvers.
"""
function solved_unknowns(sys::AbstractSystem)
    sts = unknowns(sys)
    if has_solved_unknowns(sys)
        sts = something(get_solved_unknowns(sys), sts)
    end
    return sts
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
end

function ObservedFunctionCache(
        sys; steady_state = false, eval_expression = false,
        eval_module = @__MODULE__, checkbounds = true)
    return ObservedFunctionCache(
        sys, Dict(), steady_state, eval_expression, eval_module, checkbounds)
end

# This is hit because ensemble problems do a deepcopy
function Base.deepcopy_internal(ofc::ObservedFunctionCache, stackdict::IdDict)
    sys = deepcopy(ofc.sys)
    dict = deepcopy(ofc.dict)
    steady_state = ofc.steady_state
    eval_expression = ofc.eval_expression
    eval_module = ofc.eval_module
    checkbounds = ofc.checkbounds
    newofc = ObservedFunctionCache(
        sys, dict, steady_state, eval_expression, eval_module, checkbounds)
    stackdict[ofc] = newofc
    return newofc
end

function (ofc::ObservedFunctionCache)(obsvar, args...)
    obs = get!(ofc.dict, value(obsvar)) do
        SymbolicIndexingInterface.observed(
            ofc.sys, obsvar; eval_expression = ofc.eval_expression,
            eval_module = ofc.eval_module, checkbounds = ofc.checkbounds)
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

function push_defaults!(stmt, defs, var2name)
    defs_name = gensym(:defs)
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
    defs_name = push_defaults!(stmt, defaults(sys), var2name)
    obs_name = push_eqs!(stmt, obs, var2name)

    if sys isa ODESystem
        iv = get_iv(sys)
        ivname = gensym(:iv)
        push!(stmt, :($ivname = (@variables $(getname(iv)))[1]))
        push!(stmt,
            :($ODESystem($eqs_name, $ivname, $stsname, $psname; defaults = $defs_name,
                observed = $obs_name,
                name = $name, checks = false)))
    elseif sys isa NonlinearSystem
        push!(stmt,
            :($NonlinearSystem($eqs_name, $stsname, $psname; defaults = $defs_name,
                observed = $obs_name,
                name = $name, checks = false)))
    end

    expr = :(let
        $expr
    end)
    Base.remove_linenums!(expr) # keeping the line numbers is never helpful
end

Base.write(io::IO, sys::AbstractSystem) = write(io, readable_code(toexpr(sys)))

function get_or_construct_tearing_state(sys)
    if has_tearing_state(sys)
        state = get_tearing_state(sys)
        if state === nothing
            state = TearingState(sys)
        end
    else
        state = nothing
    end
    state
end

"""
    n_expanded_connection_equations(sys::AbstractSystem)

Returns the number of equations that the connections in `sys` expands to.
Equivalent to `length(equations(expand_connections(sys))) - length(filter(eq -> !(eq.lhs isa Connection), equations(sys)))`.
"""
function n_expanded_connection_equations(sys::AbstractSystem)
    # TODO: what about inputs?
    isconnector(sys) && return length(get_unknowns(sys))
    sys = remove_analysis_points(sys)
    n_variable_connect_eqs = 0
    for eq in equations(sys)
        is_causal_variable_connection(eq.rhs) || continue
        n_variable_connect_eqs += length(get_systems(eq.rhs)) - 1
    end

    sys, (csets, _) = generate_connection_set(sys)
    ceqs, instream_csets = generate_connection_equations_and_stream_connections(csets)
    n_outer_stream_variables = 0
    for cset in instream_csets
        n_outer_stream_variables += count(x -> x.isouter, cset.set)
    end

    #n_toplevel_unused_flows = 0
    #toplevel_flows = Set()
    #for cset in csets
    #    e1 = first(cset.set)
    #    e1.sys.namespace === nothing || continue
    #    for e in cset.set
    #        get_connection_type(e.v) === Flow || continue
    #        push!(toplevel_flows, e.v)
    #    end
    #end
    #for m in get_systems(sys)
    #    isconnector(m) || continue
    #    n_toplevel_unused_flows += count(x->get_connection_type(x) === Flow && !(x in toplevel_flows), get_unknowns(m))
    #end

    nextras = n_outer_stream_variables + length(ceqs) + n_variable_connect_eqs
end

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

    # Print observed
    nobs = has_observed(sys) ? length(observed(sys)) : 0
    nobs > 0 && printstyled(io, "\nObserved ($nobs):"; bold)
    nobs > 0 && hint && print(io, " see observed($name)")

    return nothing
end

function Graphs.incidence_matrix(sys)
    if has_torn_matching(sys) && has_tearing_state(sys)
        state = get_tearing_state(sys)
        incidence_matrix(state.structure.graph, Num(Sym{Real}(:×)))
    else
        return nothing
    end
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

function single_named_expr(expr)
    name, call = split_assign(expr)
    if Meta.isexpr(name, :ref)
        name, idxs = name.args
        check_name(name)
        var = gensym(name)
        ex = quote
            $var = $(_named(name, call))
            $name = map(i -> $rename($var, Symbol($(Meta.quot(name)), :_, i)), $idxs)
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
        if !hasmetadata(uv, SymScope)
            ParentScope(sym)
        else
            sym
        end
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
                $Setfield.@set!(res.gui_metadata=$GUIMetadata($GlobalRef(@__MODULE__, name)))
            else
                res
            end
        end
    end
end

macro component(expr)
    esc(component_post_processing(expr, false))
end

macro mtkbuild(exprs...)
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
    call_expr = Expr(:call, structural_simplify, kwargs..., name)
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
        error("debug_system(sys) only works on systems with no sub-systems! Consider flattening it with flatten(sys) or structural_simplify(sys) first.")
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

function eliminate_constants(sys::AbstractSystem)
    if has_eqs(sys)
        eqs = get_eqs(sys)
        eq_cs = collect_constants(eqs)
        if !isempty(eq_cs)
            new_eqs = eliminate_constants(eqs, eq_cs)
            @set! sys.eqs = new_eqs
        end
    end
    return sys
end

function io_preprocessing(sys::AbstractSystem, inputs,
        outputs; simplify = false, kwargs...)
    sys, input_idxs = structural_simplify(sys, (inputs, outputs); simplify, kwargs...)

    eqs = equations(sys)
    alg_start_idx = findfirst(!isdiffeq, eqs)
    if alg_start_idx === nothing
        alg_start_idx = length(eqs) + 1
    end
    diff_idxs = 1:(alg_start_idx - 1)
    alge_idxs = alg_start_idx:length(eqs)

    sys, diff_idxs, alge_idxs, input_idxs
end

"""
    lin_fun, simplified_sys = linearization_function(sys::AbstractSystem, inputs, outputs; simplify = false, initialize = true, initialization_solver_alg = TrustRegion(), kwargs...)

Return a function that linearizes the system `sys`. The function [`linearize`](@ref) provides a higher-level and easier to use interface.

`lin_fun` is a function `(variables, p, t) -> (; f_x, f_z, g_x, g_z, f_u, g_u, h_x, h_z, h_u)`, i.e., it returns a NamedTuple with the Jacobians of `f,g,h` for the nonlinear `sys` (technically for `simplified_sys`) on the form

```math
\\begin{aligned}
ẋ &= f(x, z, u) \\\\
0 &= g(x, z, u) \\\\
y &= h(x, z, u)
\\end{aligned}
```

where `x` are differential unknown variables, `z` algebraic variables, `u` inputs and `y` outputs. To obtain a linear statespace representation, see [`linearize`](@ref). The input argument `variables` is a vector defining the operating point, corresponding to `unknowns(simplified_sys)` and `p` is a vector corresponding to the parameters of `simplified_sys`. Note: all variables in `inputs` have been converted to parameters in `simplified_sys`.

The `simplified_sys` has undergone [`structural_simplify`](@ref) and had any occurring input or output variables replaced with the variables provided in arguments `inputs` and `outputs`. The unknowns of this system also indicate the order of the unknowns that holds for the linearized matrices.

# Arguments:

  - `sys`: An [`ODESystem`](@ref). This function will automatically apply simplification passes on `sys` and return the resulting `simplified_sys`.
  - `inputs`: A vector of variables that indicate the inputs of the linearized input-output model.
  - `outputs`: A vector of variables that indicate the outputs of the linearized input-output model.
  - `simplify`: Apply simplification in tearing.
  - `initialize`: If true, a check is performed to ensure that the operating point is consistent (satisfies algebraic equations). If the op is not consistent, initialization is performed.
  - `initialization_solver_alg`: A NonlinearSolve algorithm to use for solving for a feasible set of state and algebraic variables that satisfies the specified operating point.
  - `kwargs`: Are passed on to `find_solvables!`

See also [`linearize`](@ref) which provides a higher-level interface.
"""
function linearization_function(sys::AbstractSystem, inputs,
        outputs; simplify = false,
        initialize = true,
        op = Dict(),
        p = DiffEqBase.NullParameters(),
        zero_dummy_der = false,
        initialization_solver_alg = TrustRegion(),
        eval_expression = false, eval_module = @__MODULE__,
        warn_initialize_determined = true,
        kwargs...)
    op = Dict(op)
    inputs isa AbstractVector || (inputs = [inputs])
    outputs isa AbstractVector || (outputs = [outputs])
    inputs = mapreduce(vcat, inputs; init = []) do var
        symbolic_type(var) == ArraySymbolic() ? collect(var) : [var]
    end
    outputs = mapreduce(vcat, outputs; init = []) do var
        symbolic_type(var) == ArraySymbolic() ? collect(var) : [var]
    end
    ssys, diff_idxs, alge_idxs, input_idxs = io_preprocessing(sys, inputs, outputs;
        simplify,
        kwargs...)
    if zero_dummy_der
        dummyder = setdiff(unknowns(ssys), unknowns(sys))
        defs = Dict(x => 0.0 for x in dummyder)
        @set! ssys.defaults = merge(defs, defaults(ssys))
        op = merge(defs, op)
    end
    sys = ssys
    u0map = Dict(k => v for (k, v) in op if is_variable(ssys, k))
    initsys = structural_simplify(
        generate_initializesystem(
            sys, u0map = u0map, guesses = guesses(sys), algebraic_only = true),
        fully_determined = false)

    # HACK: some unknowns may not be involved in any initialization equations, and are
    # thus removed from the system during `structural_simplify`.
    # This causes `getu(initsys, unknowns(sys))` to fail, so we add them back as parameters
    # for now.
    missing_unknowns = setdiff(unknowns(sys), all_symbols(initsys))
    if !isempty(missing_unknowns)
        if warn_initialize_determined
            @warn "Initialization system is underdetermined. No equations for $(missing_unknowns). Initialization will default to using least squares. To suppress this warning pass warn_initialize_determined = false."
        end
        new_parameters = [parameters(initsys); missing_unknowns]
        @set! initsys.ps = new_parameters
        initsys = complete(initsys)
    end

    if p isa SciMLBase.NullParameters
        p = Dict()
    else
        p = todict(p)
    end
    x0 = merge(defaults_and_guesses(sys), op)
    if has_index_cache(sys) && get_index_cache(sys) !== nothing
        sys_ps = MTKParameters(sys, p, x0)
    else
        sys_ps = varmap_to_vars(p, parameters(sys); defaults = x0)
    end
    p[get_iv(sys)] = NaN
    if has_index_cache(initsys) && get_index_cache(initsys) !== nothing
        oldps = MTKParameters(initsys, p, merge(guesses(sys), defaults(sys), op))
        initsys_ps = parameters(initsys)
        p_getter = build_explicit_observed_function(
            sys, initsys_ps; eval_expression, eval_module)

        u_getter = isempty(unknowns(initsys)) ? (_...) -> nothing :
                   build_explicit_observed_function(
            sys, unknowns(initsys); eval_expression, eval_module)
        get_initprob_u_p = let p_getter = p_getter,
            p_setter! = setp(initsys, initsys_ps),
            u_getter = u_getter

            function (u, p, t)
                p_setter!(oldps, p_getter(u, p, t))
                newu = u_getter(u, p, t)
                return newu, oldps
            end
        end
    else
        get_initprob_u_p = let p_getter = getu(sys, parameters(initsys)),
            u_getter = build_explicit_observed_function(
                sys, unknowns(initsys); eval_expression, eval_module)

            function (u, p, t)
                state = ProblemState(; u, p, t)
                return u_getter(
                    state_values(state), parameter_values(state), current_time(state)),
                p_getter(state)
            end
        end
    end
    initfn = NonlinearFunction(initsys; eval_expression, eval_module)
    initprobmap = build_explicit_observed_function(
        initsys, unknowns(sys); eval_expression, eval_module)
    ps = parameters(sys)
    h = build_explicit_observed_function(sys, outputs; eval_expression, eval_module)
    lin_fun = let diff_idxs = diff_idxs,
        alge_idxs = alge_idxs,
        input_idxs = input_idxs,
        sts = unknowns(sys),
        get_initprob_u_p = get_initprob_u_p,
        fun = ODEFunction{true, SciMLBase.FullSpecialize}(
            sys, unknowns(sys), ps; eval_expression, eval_module),
        initfn = initfn,
        initprobmap = initprobmap,
        h = h,
        chunk = ForwardDiff.Chunk(input_idxs),
        sys_ps = sys_ps,
        initialize = initialize,
        initialization_solver_alg = initialization_solver_alg,
        sys = sys

        function (u, p, t)
            if !isa(p, MTKParameters)
                p = todict(p)
                newps = deepcopy(sys_ps)
                for (k, v) in p
                    if is_parameter(sys, k)
                        v = fixpoint_sub(v, p)
                        setp(sys, k)(newps, v)
                    end
                end
                p = newps
            end

            if u !== nothing # Handle systems without unknowns
                length(sts) == length(u) ||
                    error("Number of unknown variables ($(length(sts))) does not match the number of input unknowns ($(length(u)))")
                if initialize && !isempty(alge_idxs) # This is expensive and can be omitted if the user knows that the system is already initialized
                    residual = fun(u, p, t)
                    if norm(residual[alge_idxs]) > √(eps(eltype(residual)))
                        initu0, initp = get_initprob_u_p(u, p, t)
                        initprob = NonlinearLeastSquaresProblem(initfn, initu0, initp)
                        nlsol = solve(initprob, initialization_solver_alg)
                        u = initprobmap(state_values(nlsol), parameter_values(nlsol))
                    end
                end
                uf = SciMLBase.UJacobianWrapper(fun, t, p)
                fg_xz = ForwardDiff.jacobian(uf, u)
                h_xz = ForwardDiff.jacobian(
                    let p = p, t = t
                        xz -> h(xz, p, t)
                    end, u)
                pf = SciMLBase.ParamJacobianWrapper(fun, t, u)
                fg_u = jacobian_wrt_vars(pf, p, input_idxs, chunk)
            else
                length(sts) == 0 ||
                    error("Number of unknown variables (0) does not match the number of input unknowns ($(length(u)))")
                fg_xz = zeros(0, 0)
                h_xz = fg_u = zeros(0, length(inputs))
            end
            hp = let u = u, t = t
                _hp(p) = h(u, p, t)
                _hp
            end
            h_u = jacobian_wrt_vars(hp, p, input_idxs, chunk)
            (f_x = fg_xz[diff_idxs, diff_idxs],
                f_z = fg_xz[diff_idxs, alge_idxs],
                g_x = fg_xz[alge_idxs, diff_idxs],
                g_z = fg_xz[alge_idxs, alge_idxs],
                f_u = fg_u[diff_idxs, :],
                g_u = fg_u[alge_idxs, :],
                h_x = h_xz[:, diff_idxs],
                h_z = h_xz[:, alge_idxs],
                h_u = h_u)
        end
    end
    return lin_fun, sys
end

"""
    (; A, B, C, D), simplified_sys = linearize_symbolic(sys::AbstractSystem, inputs, outputs; simplify = false, allow_input_derivatives = false, kwargs...)

Similar to [`linearize`](@ref), but returns symbolic matrices `A,B,C,D` rather than numeric. While `linearize` uses ForwardDiff to perform the linearization, this function uses `Symbolics.jacobian`.

See [`linearize`](@ref) for a description of the arguments.

# Extended help
The named tuple returned as the first argument additionally contains the jacobians `f_x, f_z, g_x, g_z, f_u, g_u, h_x, h_z, h_u` of
```math
\\begin{aligned}
ẋ &= f(x, z, u) \\\\
0 &= g(x, z, u) \\\\
y &= h(x, z, u)
\\end{aligned}
```
where `x` are differential unknown variables, `z` algebraic variables, `u` inputs and `y` outputs.
"""
function linearize_symbolic(sys::AbstractSystem, inputs,
        outputs; simplify = false, allow_input_derivatives = false,
        eval_expression = false, eval_module = @__MODULE__,
        kwargs...)
    sys, diff_idxs, alge_idxs, input_idxs = io_preprocessing(
        sys, inputs, outputs; simplify,
        kwargs...)
    sts = unknowns(sys)
    t = get_iv(sys)
    ps = parameters(sys)
    p = reorder_parameters(sys, ps)

    fun_expr = generate_function(sys, sts, ps; expression = Val{true})[1]
    fun = eval_or_rgf(fun_expr; eval_expression, eval_module)
    dx = fun(sts, p, t)

    h = build_explicit_observed_function(sys, outputs; eval_expression, eval_module)
    y = h(sts, p, t)

    fg_xz = Symbolics.jacobian(dx, sts)
    fg_u = Symbolics.jacobian(dx, inputs)
    h_xz = Symbolics.jacobian(y, sts)
    h_u = Symbolics.jacobian(y, inputs)
    f_x = fg_xz[diff_idxs, diff_idxs]
    f_z = fg_xz[diff_idxs, alge_idxs]
    g_x = fg_xz[alge_idxs, diff_idxs]
    g_z = fg_xz[alge_idxs, alge_idxs]
    f_u = fg_u[diff_idxs, :]
    g_u = fg_u[alge_idxs, :]
    h_x = h_xz[:, diff_idxs]
    h_z = h_xz[:, alge_idxs]

    nx, nu = size(f_u)
    nz = size(f_z, 2)
    ny = size(h_x, 1)

    D = h_u

    if isempty(g_z) # ODE
        A = f_x
        B = f_u
        C = h_x
    else
        gz = lu(g_z; check = false)
        issuccess(gz) ||
            error("g_z not invertible, this indicates that the DAE is of index > 1.")
        gzgx = -(gz \ g_x)
        A = [f_x f_z
             gzgx*f_x gzgx*f_z]
        B = [f_u
             gzgx * f_u] # The cited paper has zeros in the bottom block, see derivation in https://github.com/SciML/ModelingToolkit.jl/pull/1691 for the correct formula

        C = [h_x h_z]
        Bs = -(gz \ g_u) # This equation differ from the cited paper, the paper is likely wrong since their equaiton leads to a dimension mismatch.
        if !iszero(Bs)
            if !allow_input_derivatives
                der_inds = findall(vec(any(!iszero, Bs, dims = 1)))
                @show typeof(der_inds)
                error("Input derivatives appeared in expressions (-g_z\\g_u != 0), the following inputs appeared differentiated: $(ModelingToolkit.inputs(sys)[der_inds]). Call `linearize_symbolic` with keyword argument `allow_input_derivatives = true` to allow this and have the returned `B` matrix be of double width ($(2nu)), where the last $nu inputs are the derivatives of the first $nu inputs.")
            end
            B = [B [zeros(nx, nu); Bs]]
            D = [D zeros(ny, nu)]
        end
    end

    (; A, B, C, D, f_x, f_z, g_x, g_z, f_u, g_u, h_x, h_z, h_u), sys
end

function markio!(state, orig_inputs, inputs, outputs; check = true)
    fullvars = get_fullvars(state)
    inputset = Dict{Any, Bool}(i => false for i in inputs)
    outputset = Dict{Any, Bool}(o => false for o in outputs)
    for (i, v) in enumerate(fullvars)
        if v in keys(inputset)
            if v in keys(outputset)
                v = setio(v, true, true)
                outputset[v] = true
            else
                v = setio(v, true, false)
            end
            inputset[v] = true
            fullvars[i] = v
        elseif v in keys(outputset)
            v = setio(v, false, true)
            outputset[v] = true
            fullvars[i] = v
        else
            if isinput(v)
                push!(orig_inputs, v)
            end
            v = setio(v, false, false)
            fullvars[i] = v
        end
    end
    if check
        ikeys = keys(filter(!last, inputset))
        if !isempty(ikeys)
            error(
                "Some specified inputs were not found in system. The following variables were not found ",
                ikeys)
        end
    end
    check && (all(values(outputset)) ||
     error(
        "Some specified outputs were not found in system. The following Dict indicates the found variables ",
        outputset))
    state, orig_inputs
end

"""
    (; A, B, C, D), simplified_sys = linearize(sys, inputs, outputs;    t=0.0, op = Dict(), allow_input_derivatives = false, zero_dummy_der=false, kwargs...)
    (; A, B, C, D)                 = linearize(simplified_sys, lin_fun; t=0.0, op = Dict(), allow_input_derivatives = false, zero_dummy_der=false)

Linearize `sys` between `inputs` and `outputs`, both vectors of variables. Return a NamedTuple with the matrices of a linear statespace representation
on the form

```math
\\begin{aligned}
ẋ &= Ax + Bu\\\\
y &= Cx + Du
\\end{aligned}
```

The first signature automatically calls [`linearization_function`](@ref) internally,
while the second signature expects the outputs of [`linearization_function`](@ref) as input.

`op` denotes the operating point around which to linearize. If none is provided,
the default values of `sys` are used.

If `allow_input_derivatives = false`, an error will be thrown if input derivatives (``u̇``) appear as inputs in the linearized equations. If input derivatives are allowed, the returned `B` matrix will be of double width, corresponding to the input `[u; u̇]`.

`zero_dummy_der` can be set to automatically set the operating point to zero for all dummy derivatives.

See also [`linearization_function`](@ref) which provides a lower-level interface, [`linearize_symbolic`](@ref) and [`ModelingToolkit.reorder_unknowns`](@ref).

See extended help for an example.

The implementation and notation follows that of
["Linear Analysis Approach for Modelica Models", Allain et al. 2009](https://ep.liu.se/ecp/043/075/ecp09430097.pdf)

# Extended help

This example builds the following feedback interconnection and linearizes it from the input of `F` to the output of `P`.

```

  r ┌─────┐       ┌─────┐     ┌─────┐
───►│     ├──────►│     │  u  │     │
    │  F  │       │  C  ├────►│  P  │ y
    └─────┘     ┌►│     │     │     ├─┬─►
                │ └─────┘     └─────┘ │
                │                     │
                └─────────────────────┘
```

```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
function plant(; name)
    @variables x(t) = 1
    @variables u(t)=0 y(t)=0
    eqs = [D(x) ~ -x + u
           y ~ x]
    ODESystem(eqs, t; name = name)
end

function ref_filt(; name)
    @variables x(t)=0 y(t)=0
    @variables u(t)=0 [input = true]
    eqs = [D(x) ~ -2 * x + u
           y ~ x]
    ODESystem(eqs, t, name = name)
end

function controller(kp; name)
    @variables y(t)=0 r(t)=0 u(t)=0
    @parameters kp = kp
    eqs = [
        u ~ kp * (r - y),
    ]
    ODESystem(eqs, t; name = name)
end

@named f = ref_filt()
@named c = controller(1)
@named p = plant()

connections = [f.y ~ c.r # filtered reference to controller reference
               c.u ~ p.u # controller output to plant input
               p.y ~ c.y]

@named cl = ODESystem(connections, t, systems = [f, c, p])

lsys0, ssys = linearize(cl, [f.u], [p.x])
desired_order = [f.x, p.x]
lsys = ModelingToolkit.reorder_unknowns(lsys0, unknowns(ssys), desired_order)

@assert lsys.A == [-2 0; 1 -2]
@assert lsys.B == [1; 0;;]
@assert lsys.C == [0 1]
@assert lsys.D[] == 0

## Symbolic linearization
lsys_sym, _ = ModelingToolkit.linearize_symbolic(cl, [f.u], [p.x])

@assert substitute(lsys_sym.A, ModelingToolkit.defaults(cl)) == lsys.A
```
"""
function linearize(sys, lin_fun; t = 0.0, op = Dict(), allow_input_derivatives = false,
        p = DiffEqBase.NullParameters())
    x0 = merge(defaults(sys), Dict(missing_variable_defaults(sys)), op)
    u0, defs = get_u0(sys, x0, p)
    if has_index_cache(sys) && get_index_cache(sys) !== nothing
        if p isa SciMLBase.NullParameters
            p = op
        elseif p isa Dict
            p = merge(p, op)
        elseif p isa Vector && eltype(p) <: Pair
            p = merge(Dict(p), op)
        elseif p isa Vector
            p = merge(Dict(parameters(sys) .=> p), op)
        end
    end
    linres = lin_fun(u0, p, t)
    f_x, f_z, g_x, g_z, f_u, g_u, h_x, h_z, h_u = linres

    nx, nu = size(f_u)
    nz = size(f_z, 2)
    ny = size(h_x, 1)

    D = h_u

    if isempty(g_z)
        A = f_x
        B = f_u
        C = h_x
        @assert iszero(g_x)
        @assert iszero(g_z)
        @assert iszero(g_u)
    else
        gz = lu(g_z; check = false)
        issuccess(gz) ||
            error("g_z not invertible, this indicates that the DAE is of index > 1.")
        gzgx = -(gz \ g_x)
        A = [f_x f_z
             gzgx*f_x gzgx*f_z]
        B = [f_u
             gzgx * f_u] # The cited paper has zeros in the bottom block, see derivation in https://github.com/SciML/ModelingToolkit.jl/pull/1691 for the correct formula

        C = [h_x h_z]
        Bs = -(gz \ g_u) # This equation differ from the cited paper, the paper is likely wrong since their equaiton leads to a dimension mismatch.
        if !iszero(Bs)
            if !allow_input_derivatives
                der_inds = findall(vec(any(!=(0), Bs, dims = 1)))
                error("Input derivatives appeared in expressions (-g_z\\g_u != 0), the following inputs appeared differentiated: $(inputs(sys)[der_inds]). Call `linearize` with keyword argument `allow_input_derivatives = true` to allow this and have the returned `B` matrix be of double width ($(2nu)), where the last $nu inputs are the derivatives of the first $nu inputs.")
            end
            B = [B [zeros(nx, nu); Bs]]
            D = [D zeros(ny, nu)]
        end
    end

    (; A, B, C, D)
end

function linearize(sys, inputs, outputs; op = Dict(), t = 0.0,
        allow_input_derivatives = false,
        zero_dummy_der = false,
        kwargs...)
    lin_fun, ssys = linearization_function(sys,
        inputs,
        outputs;
        zero_dummy_der,
        op,
        kwargs...)
    linearize(ssys, lin_fun; op, t, allow_input_derivatives), ssys
end

"""
    (; Ã, B̃, C̃, D̃) = similarity_transform(sys, T; unitary=false)

Perform a similarity transform `T : Tx̃ = x` on linear system represented by matrices in NamedTuple `sys` such that

```
Ã = T⁻¹AT
B̃ = T⁻¹ B
C̃ = CT
D̃ = D
```

If `unitary=true`, `T` is assumed unitary and the matrix adjoint is used instead of the inverse.
"""
function similarity_transform(sys::NamedTuple, T; unitary = false)
    if unitary
        A = T'sys.A * T
        B = T'sys.B
    else
        Tf = lu(T)
        A = Tf \ sys.A * T
        B = Tf \ sys.B
    end
    C = sys.C * T
    D = sys.D
    (; A, B, C, D)
end

"""
    reorder_unknowns(sys::NamedTuple, old, new)

Permute the state representation of `sys` obtained from [`linearize`](@ref) so that the state unknown is changed from `old` to `new`
Example:

```
lsys, ssys = linearize(pid, [reference.u, measurement.u], [ctr_output.u])
desired_order = [int.x, der.x] # Unknowns that are present in unknowns(ssys)
lsys = ModelingToolkit.reorder_unknowns(lsys, unknowns(ssys), desired_order)
```

See also [`ModelingToolkit.similarity_transform`](@ref)
"""
function reorder_unknowns(sys::NamedTuple, old, new)
    nx = length(old)
    length(new) == nx || error("old and new must have the same length")
    perm = [findfirst(isequal(n), old) for n in new]
    issorted(perm) && return sys # shortcut return, no reordering
    P = zeros(Int, nx, nx)
    for i in 1:nx # Build permutation matrix
        P[i, perm[i]] = 1
    end
    similarity_transform(sys, P; unitary = true)
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
        throw(ArgumentError("The system has array equations. Call `structural_simplify` to handle such equations or scalarize them manually."))
    end
    if any(x -> Symbolics.isarraysymbolic(x), dvs)
        throw(ArgumentError("The system has array unknowns. Call `structural_simplify` to handle this or scalarize them manually."))
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
function Base.hash(sys::AbstractSystem, s::UInt)
    s = hash(nameof(sys), s)
    s = foldr(hash, get_systems(sys), init = s)
    s = foldr(hash, get_unknowns(sys), init = s)
    s = foldr(hash, get_ps(sys), init = s)
    if sys isa OptimizationSystem
        s = hash(get_op(sys), s)
    else
        s = foldr(hash, get_eqs(sys), init = s)
    end
    s = foldr(hash, get_observed(sys), init = s)
    s = foldr(hash, get_continuous_events(sys), init = s)
    s = foldr(hash, get_discrete_events(sys), init = s)
    s = hash(independent_variables(sys), s)
    return s
end

"""
$(TYPEDSIGNATURES)

Extend `basesys` with `sys`.
By default, the resulting system inherits `sys`'s name and description.

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
    dep_ps = union(parameter_dependencies(basesys), parameter_dependencies(sys))
    obs = union(get_observed(basesys), get_observed(sys))
    cevs = union(get_continuous_events(basesys), get_continuous_events(sys))
    devs = union(get_discrete_events(basesys), get_discrete_events(sys))
    defs = merge(get_defaults(basesys), get_defaults(sys)) # prefer `sys`
    meta = union_nothing(get_metadata(basesys), get_metadata(sys))
    syss = union(get_systems(basesys), get_systems(sys))
    args = length(ivs) == 0 ? (eqs, sts, ps) : (eqs, ivs[1], sts, ps)
    kwargs = (parameter_dependencies = dep_ps, observed = obs, continuous_events = cevs,
        discrete_events = devs, defaults = defs, systems = syss, metadata = meta,
        name = name, description = description, gui_metadata = gui_metadata)

    # collect fields specific to some system types
    if basesys isa ODESystem
        ieqs = union(get_initialization_eqs(basesys), get_initialization_eqs(sys))
        guesses = merge(get_guesses(basesys), get_guesses(sys)) # prefer `sys`
        kwargs = merge(kwargs, (initialization_eqs = ieqs, guesses = guesses))
    end

    if has_assertions(basesys)
        kwargs = merge(
            kwargs, (; assertions = merge(get_assertions(basesys), get_assertions(sys))))
    end

    return T(args...; kwargs...)
end

function extend(sys, basesys::Vector{T}) where {T <: AbstractSystem}
    foldl(extend, basesys, init = sys)
end

function Base.:(&)(sys::AbstractSystem, basesys::AbstractSystem; kwargs...)
    extend(sys, basesys; kwargs...)
end

function Base.:(&)(
        sys::AbstractSystem, basesys::Vector{T}; kwargs...) where {T <: AbstractSystem}
    extend(sys, basesys; kwargs...)
end

"""
$(SIGNATURES)

Compose multiple systems together. The resulting system would inherit the first
system's name.

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
    return sys
end
function compose(syss...; name = nameof(first(syss)))
    compose(first(syss), collect(syss[2:end]); name = name)
end
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
    elseif sys isa ODESystem
        rules = todict(map(r -> Symbolics.unwrap(r[1]) => Symbolics.unwrap(r[2]),
            collect(rules)))
        eqs = fast_substitute(get_eqs(sys), rules)
        pdeps = fast_substitute(get_parameter_dependencies(sys), rules)
        defs = Dict(fast_substitute(k, rules) => fast_substitute(v, rules)
        for (k, v) in get_defaults(sys))
        guess = Dict(fast_substitute(k, rules) => fast_substitute(v, rules)
        for (k, v) in get_guesses(sys))
        subsys = map(s -> substitute(s, rules), get_systems(sys))
        ODESystem(eqs, get_iv(sys); name = nameof(sys), defaults = defs,
            guesses = guess, parameter_dependencies = pdeps, systems = subsys)
    else
        error("substituting symbols is not supported for $(typeof(sys))")
    end
end

struct InvalidParameterDependenciesType
    got::Any
end

function Base.showerror(io::IO, err::InvalidParameterDependenciesType)
    print(
        io, "Parameter dependencies must be a `Dict`, or an array of `Pair` or `Equation`.")
    if err.got !== nothing
        print(io, " Got ", err.got)
    end
end

function process_parameter_dependencies(pdeps, ps)
    if pdeps === nothing || isempty(pdeps)
        return Equation[], ps
    end
    if pdeps isa Dict
        pdeps = [k ~ v for (k, v) in pdeps]
    else
        pdeps isa AbstractArray || throw(InvalidParameterDependenciesType(pdeps))
        pdeps = [if p isa Pair
                     p[1] ~ p[2]
                 elseif p isa Equation
                     p
                 else
                     error("Parameter dependencies must be a `Dict`, `Vector{Pair}` or `Vector{Equation}`")
                 end
                 for p in pdeps]
    end
    lhss = []
    for p in pdeps
        if !isparameter(p.lhs)
            error("LHS of parameter dependency must be a single parameter. Found $(p.lhs).")
        end
        syms = vars(p.rhs)
        if !all(isparameter, syms)
            error("RHS of parameter dependency must only include parameters. Found $(p.rhs)")
        end
        push!(lhss, p.lhs)
    end
    lhss = map(identity, lhss)
    pdeps = topsort_equations(pdeps, union(ps, lhss))
    ps = filter(ps) do p
        !any(isequal(p), lhss)
    end
    return pdeps, ps
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
@named sys = ODESystem(Equation[], t, [x], [p, q])
ModelingToolkit.dump_parameters(sys)
```

See also: [`ModelingToolkit.dump_variable_metadata`](@ref), [`ModelingToolkit.dump_unknowns`](@ref)
"""
function dump_parameters(sys::AbstractSystem)
    defs = defaults(sys)
    pdeps = parameter_dependencies(sys)
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
@named sys = ODESystem(Equation[], t, [x], [p, q])
ModelingToolkit.dump_unknowns(sys)
```

See also: [`ModelingToolkit.dump_variable_metadata`](@ref), [`ModelingToolkit.dump_parameters`](@ref)
"""
function dump_unknowns(sys::AbstractSystem)
    defs = varmap_with_toterm(defaults(sys))
    gs = varmap_with_toterm(guesses(sys))
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
    isdefined(eq, :lhs) && hasnode(is_derivative, wrap(eq.lhs)) &&
        (return true)
    isdefined(eq, :rhs) && hasnode(is_derivative, wrap(eq.rhs)) &&
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
@named osys = ODESystem([eq1, eq2], t)

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
@named osys = ODESystem([eq1, eq2], t)

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
@named osys1 = ODESystem([eq1], t)
@named osys2 = ODESystem([eq2], t)

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
@named osys1 = ODESystem([eq1], t)
@named osys2 = ODESystem([eq2], t)

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
@named osys1 = ODESystem([eq1], t)
@named osys2 = ODESystem([eq2], t)
osys12 = compose(osys1, [osys2])
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
@named osys1 = ODESystem([eq1], t)
@named osys2 = ODESystem([eq2], t)
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
@named osys1 = ODESystem([eq1], t)
@named osys2 = ODESystem([eq2], t)
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
@named osys1 = ODESystem([eq1], t)
@named osys2 = ODESystem([eq2], t)
osys12 = compose(osys1, [osys2])
osys21 = compose(osys2, [osys1])

has_diff_eqs(osys12) # returns `true`.
has_diff_eqs(osys21) # returns `false`.
```
"""
has_diff_eqs(sys::AbstractSystem) = any(is_diff_equation, get_eqs(sys))
