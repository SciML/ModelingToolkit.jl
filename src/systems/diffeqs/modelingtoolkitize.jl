"""
$(TYPEDSIGNATURES)

Generate `ODESystem`, dependent variables, and parameters from an `ODEProblem`.
"""
function modelingtoolkitize(prob::DiffEqBase.ODEProblem; kwargs...)
    prob.f isa DiffEqBase.AbstractParameterizedFunction &&
        return prob.f.sys
    @parameters t

    p = prob.p
    has_p = !(p isa Union{DiffEqBase.NullParameters, Nothing})

    _vars = define_vars(prob.u0, t)

    vars = prob.u0 isa Number ? _vars : ArrayInterface.restructure(prob.u0, _vars)
    params = if has_p
        _params = define_params(p)
        p isa Number ? _params[1] :
        (p isa Tuple || p isa NamedTuple || p isa AbstractDict ? _params :
         ArrayInterface.restructure(p, _params))
    else
        []
    end

    var_set = Set(vars)

    D = Differential(t)
    mm = prob.f.mass_matrix

    if mm === I
        lhs = map(v -> D(v), vars)
    else
        lhs = map(mm * vars) do v
            if iszero(v)
                0
            elseif v in var_set
                D(v)
            else
                error("Non-permuation mass matrix is not supported.")
            end
        end
    end

    if DiffEqBase.isinplace(prob)
        rhs = ArrayInterface.restructure(prob.u0, similar(vars, Num))
        fill!(rhs, 0)
        if prob.f isa ODEFunction &&
           prob.f.f isa FunctionWrappersWrappers.FunctionWrappersWrapper
            prob.f.f.fw[1].obj[](rhs, vars, params, t)
        else
            prob.f(rhs, vars, params, t)
        end
    else
        rhs = prob.f(vars, params, t)
    end

    eqs = vcat([lhs[i] ~ rhs[i] for i in eachindex(prob.u0)]...)

    sts = vec(collect(vars))

    _params = params
    params = values(params)
    params = if params isa Number || (params isa Array && ndims(params) == 0)
        [params[1]]
    else
        vec(collect(params))
    end
    default_u0 = Dict(sts .=> vec(collect(prob.u0)))
    default_p = if has_p
        if prob.p isa AbstractDict
            Dict(v => prob.p[k] for (k, v) in pairs(_params))
        else
            Dict(params .=> vec(collect(prob.p)))
        end
    else
        Dict()
    end

    de = ODESystem(eqs, t, sts, params,
                   defaults = merge(default_u0, default_p);
                   name = gensym(:MTKizedODE),
                   tspan = prob.tspan,
                   kwargs...)

    de
end

_defvaridx(x, i) = variable(x, i, T = SymbolicUtils.FnType{Tuple, Real})
_defvar(x) = variable(x, T = SymbolicUtils.FnType{Tuple, Real})

function define_vars(u, t)
    [_defvaridx(:x, i)(t) for i in eachindex(u)]
end

function define_vars(u::Union{SLArray, LArray}, t)
    [_defvar(x)(t) for x in LabelledArrays.symnames(typeof(u))]
end

function define_vars(u::NTuple{<:Number}, t)
    tuple((_defvaridx(:x, i)(ModelingToolkit.value(t)) for i in eachindex(u))...)
end

function define_vars(u::NamedTuple, t)
    NamedTuple(x => _defvar(x)(ModelingToolkit.value(t)) for x in keys(u))
end

const PARAMETERS_NOT_SUPPORTED_MESSAGE = """
                                         The chosen parameter type is currently not supported by `modelingtoolkitize`. The
                                         current supported types are:

                                         - AbstractArrays
                                         - AbstractDicts
                                         - LabelledArrays (SLArray, LArray)
                                         - Flat tuples (tuples of numbers)
                                         - Flat named tuples (namedtuples of numbers)
                                         """

struct ModelingtoolkitizeParametersNotSupportedError <: Exception
    type::Any
end

function Base.showerror(io::IO, e::ModelingtoolkitizeParametersNotSupportedError)
    println(io, PARAMETERS_NOT_SUPPORTED_MESSAGE)
    print(io, "Parameter type: ")
    println(io, e.type)
end

function define_params(p)
    throw(ModelingtoolkitizeParametersNotSupportedError(typeof(p)))
end

function define_params(p::AbstractArray)
    [toparam(variable(:α, i)) for i in eachindex(p)]
end

function define_params(p::Number)
    [toparam(variable(:α))]
end

function define_params(p::AbstractDict)
    OrderedDict(k => toparam(variable(:α, i)) for (i, k) in zip(1:length(p), keys(p)))
end

function define_params(p::Union{SLArray, LArray})
    [toparam(variable(x)) for x in LabelledArrays.symnames(typeof(p))]
end

function define_params(p::Tuple)
    tuple((toparam(variable(:α, i)) for i in eachindex(p))...)
end

function define_params(p::NamedTuple)
    NamedTuple(x => toparam(variable(x)) for x in keys(p))
end

"""
$(TYPEDSIGNATURES)

Generate `SDESystem`, dependent variables, and parameters from an `SDEProblem`.
"""
function modelingtoolkitize(prob::DiffEqBase.SDEProblem; kwargs...)
    prob.f isa DiffEqBase.AbstractParameterizedFunction &&
        return (prob.f.sys, prob.f.sys.states, prob.f.sys.ps)
    @parameters t
    p = prob.p
    has_p = !(p isa Union{DiffEqBase.NullParameters, Nothing})

    _vars = define_vars(prob.u0, t)

    vars = prob.u0 isa Number ? _vars : ArrayInterface.restructure(prob.u0, _vars)
    params = if has_p
        _params = define_params(p)
        p isa Number ? _params[1] :
        (p isa Tuple || p isa NamedTuple ? _params :
         ArrayInterface.restructure(p, _params))
    else
        []
    end

    D = Differential(t)

    rhs = [D(var) for var in vars]

    if DiffEqBase.isinplace(prob)
        lhs = similar(vars, Any)

        prob.f(lhs, vars, params, t)

        if DiffEqBase.is_diagonal_noise(prob)
            neqs = similar(vars, Any)
            prob.g(neqs, vars, params, t)
        else
            neqs = similar(vars, Any, size(prob.noise_rate_prototype))
            prob.g(neqs, vars, params, t)
        end
    else
        lhs = prob.f(vars, params, t)
        if DiffEqBase.is_diagonal_noise(prob)
            neqs = prob.g(vars, params, t)
        else
            neqs = prob.g(vars, params, t)
        end
    end
    deqs = vcat([rhs[i] ~ lhs[i] for i in eachindex(prob.u0)]...)

    params = if ndims(params) == 0
        [params[1]]
    else
        Vector(vec(params))
    end

    de = SDESystem(deqs, neqs, t, Vector(vec(vars)), params;
                   name = gensym(:MTKizedSDE),
                   tspan = prob.tspan,
                   kwargs...)

    de
end
