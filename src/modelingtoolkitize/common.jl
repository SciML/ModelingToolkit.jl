
"""
    $(TYPEDSIGNATURES)

Check if the length of variables `vars` matches the number of names for those variables,
given by `names`. `is_unknowns` denotes whether the variable are unknowns or parameters.
"""
function varnames_length_check(vars, names; is_unknowns = false)
    length(names) == length(vars) && return
    throw(ArgumentError("""
        Number of $(is_unknowns ? "unknowns" : "parameters") ($(length(vars))) \
        does not match number of names ($(length(names))).
    """))
end

"""
    $(TYPEDSIGNATURES)

Define a subscripted time-dependent variable with name `x` and subscript `i`. Equivalent
to `@variables \$name(..)`. `T` is the desired symtype of the variable when called with
the independent variable.
"""
_defvaridx(x, i; T = Real) = variable(x, i, T = SymbolicUtils.FnType{Tuple, T})
"""
    $(TYPEDSIGNATURES)

Define a time-dependent variable with name `x`. Equivalent to `@variables \$x(..)`.
`T` is the desired symtype of the variable when called with the independent variable.
"""
_defvar(x; T = Real) = variable(x, T = SymbolicUtils.FnType{Tuple, T})

"""
    $(TYPEDSIGNATURES)

Define an array of symbolic unknowns of the appropriate type and size for `u` with
independent variable `t`.
"""
function define_vars(u, t)
    [_defvaridx(:x, i)(t) for i in eachindex(u)]
end

function define_vars(u, ::Nothing)
    [variable(:x, i) for i in eachindex(u)]
end

"""
    $(TYPEDSIGNATURES)

Return a symbolic state for the given problem `prob.`. `t` is the independent variable.
`u_names` optionally contains the names to use for the created symbolic variables.
"""
function construct_vars(prob, t, u_names = nothing)
    if prob.u0 === nothing
        return []
    end
    # construct `_vars`, AbstractSciMLFunction, AbstractSciMLFunction, a list of MTK variables for `prob.u0`.
    if u_names !== nothing
        # explicitly provided names
        varnames_length_check(state_values(prob), u_names; is_unknowns = true)
        if t === nothing
            _vars = [variable(name) for name in u_names]
        else
            _vars = [_defvar(name)(t) for name in u_names]
        end
    elseif SciMLBase.has_sys(prob.f)
        # get names from the system
        varnames = getname.(variable_symbols(prob.f.sys))
        varidxs = variable_index.((prob.f.sys,), varnames)
        invpermute!(varnames, varidxs)
        if t === nothing
            _vars = [variable(name) for name in varnames]
        else
            _vars = [_defvar(name)(t) for name in varnames]
        end
        if prob.f.sys isa System
            for (i, sym) in enumerate(variable_symbols(prob.f.sys))
                if hasbounds(sym)
                    _vars[i] = Symbolics.setmetadata(
                        _vars[i], VariableBounds, getbounds(sym))
                end
            end
        end
    else
        # auto-generate names
        _vars = define_vars(state_values(prob), t)
    end

    # Handle different types of arrays
    return prob.u0 isa Number ? _vars : ArrayInterface.restructure(prob.u0, _vars)
end

"""
    $(METHODLIST)

Define symbolic names for each value in parameter object `p`. `t` is the independent
variable of the system. `names` is a collection mapping indexes of `p` to their
names, or `nothing` to automatically generate names.

The returned value has the same structure as `p`, but symbolic variables instead of
values.
"""
function define_params(p, t, _ = nothing)
    throw(ModelingtoolkitizeParametersNotSupportedError(typeof(p)))
end

function define_params(p::AbstractArray, t, names = nothing)
    if names === nothing
        [toparam(variable(:α, i)) for i in eachindex(p)]
    else
        varnames_length_check(p, names)
        [toparam(variable(names[i])) for i in eachindex(p)]
    end
end

function define_params(p::Number, t, names = nothing)
    if names === nothing
        [toparam(variable(:α))]
    elseif names isa Union{AbstractArray, AbstractDict}
        varnames_length_check(p, names)
        [toparam(variable(names[i])) for i in eachindex(p)]
    else
        [toparam(variable(names))]
    end
end

function define_params(p::AbstractDict, t, names = nothing)
    if names === nothing
        OrderedDict(k => toparam(variable(:α, i)) for (i, k) in zip(1:length(p), keys(p)))
    else
        varnames_length_check(p, names)
        OrderedDict(k => toparam(variable(names[k])) for k in keys(p))
    end
end

function define_params(p::Tuple, t, names = nothing)
    if names === nothing
        tuple((toparam(variable(:α, i)) for i in eachindex(p))...)
    else
        varnames_length_check(p, names)
        tuple((toparam(variable(names[i])) for i in eachindex(p))...)
    end
end

function define_params(p::NamedTuple, t, names = nothing)
    if names === nothing
        NamedTuple(x => toparam(variable(x)) for x in keys(p))
    else
        varnames_length_check(p, names)
        NamedTuple(x => toparam(variable(names[x])) for x in keys(p))
    end
end

function define_params(p::MTKParameters, t, names = nothing)
    if names === nothing
        ps = []
        i = 1
        # tunables are all treated as scalar reals
        for x in p.tunable
            push!(ps, toparam(variable(:α, i)))
            i += 1
        end
        # ignore initials
        # discretes should be time-dependent
        for buf in p.discrete
            T = eltype(buf)
            for val in buf
                # respect array sizes
                shape = val isa AbstractArray ? axes(val) : nothing
                push!(ps, declare_timevarying_parameter(:α, i, t; T, shape))
                i += 1
            end
        end
        # handle constants
        for buf in p.constant
            T = eltype(buf)
            for val in buf
                # respect array sizes
                shape = val isa AbstractArray ? axes(val) : nothing
                push!(ps, declare_parameter(:α, i; T, shape))
                i += 1
            end
        end
        # handle nonnumerics
        for buf in p.nonnumeric
            T = eltype(buf)
            for val in buf
                # respect array sizes
                shape = val isa AbstractArray ? axes(val) : nothing
                push!(ps, declare_parameter(:α, i; T, shape))
                i += 1
            end
        end
        return identity.(ps)
    else
        new_p = as_any_buffer(p)
        @set! new_p.initials = []
        for (k, v) in names
            val = p[k]
            shape = val isa AbstractArray ? axes(val) : nothing
            T = typeof(val)
            if k.portion == SciMLStructures.Initials()
                continue
            end
            if k.portion == SciMLStructures.Tunable()
                T = Real
            end
            if k.portion == SciMLStructures.Discrete()
                var = declare_timevarying_parameter(getname(v), nothing, t; T, shape)
            else
                var = declare_parameter(getname(v), nothing; T, shape)
            end
            new_p[k] = var
        end
        return new_p
    end
end

"""
    $(TYPEDSIGNATURES)

Given a parameter object `p` containing symbolic variables instead of values, return
a vector of the symbolic variables.
"""
function to_paramvec(p)
    vec(collect(values(p)))
end

function to_paramvec(p::MTKParameters)
    reduce(vcat, collect(p); init = [])
end

"""
    $(TYPEDSIGNATURES)

Create a time-varying parameter with name `x`, subscript `i`, independent variable `t`
which stores values of type `T`. `shape` denotes the shape of array values, or `nothing`
for scalars.

To ignore the subscript, pass `nothing` for `i`.
"""
function declare_timevarying_parameter(x::Symbol, i, t; T, shape = nothing)
    # turn specific floating point numbers to `Real`
    if T <: Union{AbstractFloat, ForwardDiff.Dual}
        T = Real
    end
    if T <: Array{<:Union{AbstractFloat, ForwardDiff.Dual}, N} where {N}
        T = Array{Real, ndims(T)}
    end

    if i === nothing
        var = _defvar(x; T)
    else
        var = _defvaridx(x, i; T)
    end
    var = toparam(unwrap(var(t)))
    if shape !== nothing
        var = setmetadata(var, Symbolics.ArrayShapeCtx, shape)
    end
    return var
end

"""
    $(TYPEDSIGNATURES)

Create a time-varying parameter with name `x` and subscript `i`, which stores values of
type `T`. `shape` denotes the shape of array values, or `nothing` for scalars.

To ignore the subscript, pass `nothing` for `i`.
"""
function declare_parameter(x::Symbol, i; T, shape = nothing)
    # turn specific floating point numbers to `Real`
    if T <: Union{AbstractFloat, ForwardDiff.Dual}
        T = Real
    end
    if T <: Array{<:Union{AbstractFloat, ForwardDiff.Dual}, N} where {N}
        T = Array{Real, ndims(T)}
    end

    i = i === nothing ? () : (i,)
    var = toparam(unwrap(variable(x, i...; T)))
    if shape !== nothing
        var = setmetadata(var, Symbolics.ArrayShapeCtx, shape)
    end
    return var
end

"""
    $(TYPEDSIGNATURES)

Return a symbolic parameter object for the given problem `prob.`. `t` is the independent
variable. `p_names` optionally contains the names to use for the created symbolic
variables.
"""
function construct_params(prob, t, p_names = nothing)
    p = parameter_values(prob)
    has_p = !(p isa Union{DiffEqBase.NullParameters, Nothing})

    # Get names of parameters
    if has_p
        if p_names === nothing && SciMLBase.has_sys(prob.f)
            # get names from the system
            p_names = Dict(parameter_index(prob.f.sys, sym) => sym
            for sym in parameter_symbols(prob.f.sys))
        end
        params = define_params(p, t, p_names)
        if p isa Number
            params = params[1]
        elseif p isa AbstractArray
            params = ArrayInterface.restructure(p, params)
        end
    else
        params = []
    end

    return params
end

"""
    $(TYPEDSIGNATURES)

Given the differential operator `D`, mass matrix `mm` and ordered list of unknowns `vars`,
return the list of
"""
function lhs_from_mass_matrix(D, mm, vars)
    var_set = Set(vars)
    # calculate equation LHS from mass matrix
    if mm === I
        lhs = map(v -> D(v), vars)
    else
        lhs = map(mm * vars) do v
            if iszero(v)
                0
            elseif v in var_set
                D(v)
            else
                error("Non-permutation mass matrix is not supported.")
            end
        end
    end
    return lhs
end

"""
    $(TYPEDSIGNATURES)

Given a problem `prob`, the symbolic unknowns and params and the independent variable,
trace through `prob.f` and return the resultant expression.
"""
function trace_rhs(prob, vars, params, t; prototype = nothing)
    args = (vars, params)
    if t !== nothing
        args = (args..., t)
    end
    # trace prob.f to get equation RHS
    if SciMLBase.isinplace(prob.f)
        if prototype === nothing
            rhs = ArrayInterface.restructure(prob.u0, similar(vars, Num))
        else
            rhs = similar(prototype, Num)
        end
        fill!(rhs, 0)
        if prob.f isa SciMLBase.AbstractSciMLFunction &&
           prob.f.f isa FunctionWrappersWrappers.FunctionWrappersWrapper
            prob.f.f.fw[1].obj[](rhs, args...)
        else
            prob.f(rhs, args...)
        end
    else
        rhs = prob.f(args...)
    end
    return rhs
end

"""
    $(TYPEDSIGNATURES)

Obtain default values for unknowns `vars` and parameters `paramvec`
given the problem `prob` and symbolic parameter object `paramobj`.
"""
function defaults_from_u0_p(prob, vars, paramobj, paramvec)
    u0 = state_values(prob)
    p = parameter_values(prob)
    defaults = Dict{Any, Any}(vec(vars) .=> vec(collect(u0)))
    if !(p isa Union{SciMLBase.NullParameters, Nothing})
        if p isa Union{NamedTuple, AbstractDict}
            merge!(defaults, Dict(v => p[k] for (k, v) in pairs(paramobj)))
        elseif p isa MTKParameters
            pvals = [p.tunable; reduce(vcat, p.discrete; init = []);
                     reduce(vcat, p.constant; init = []);
                     reduce(vcat, p.nonnumeric; init = [])]
            merge!(defaults, Dict(paramvec .=> pvals))
        else
            merge!(defaults, Dict(paramvec .=> vec(collect(p))))
        end
    end
    return defaults
end
