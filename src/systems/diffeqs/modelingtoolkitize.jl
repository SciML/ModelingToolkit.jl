"""
$(TYPEDSIGNATURES)

Generate `ODESystem`, dependent variables, and parameters from an `ODEProblem`.
"""
function modelingtoolkitize(prob::DiffEqBase.ODEProblem; kwargs...)
    prob.f isa DiffEqBase.AbstractParameterizedFunction &&
                            return prob.f.sys
    @parameters t

    p = prob.p
    has_p = !(p isa Union{DiffEqBase.NullParameters,Nothing})

    _vars = define_vars(prob.u0,t)

    vars = prob.u0 isa Number ? _vars : ArrayInterface.restructure(prob.u0,_vars)
    params = if has_p
        _params = define_params(p)
        p isa Number ? _params[1] : (p isa Tuple || p isa NamedTuple ? _params : ArrayInterface.restructure(p,_params))
    else
        []
    end

    var_set = Set(vars)

    D = Differential(t)
    mm = prob.f.mass_matrix

    if mm === I
        lhs = map(v->D(v), vars)
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
        rhs = ArrayInterface.restructure(prob.u0,similar(vars, Num))
        prob.f(rhs, vars, params, t)
    else
        rhs = prob.f(vars, params, t)
    end

    eqs = vcat([lhs[i] ~ rhs[i] for i in eachindex(prob.u0)]...)

    sts = vec(collect(vars))

    params = if params isa Array && ndims(params) == 0
        [params[1]]
    else
        vec(collect(params))
    end
    default_u0 = Dict(sts .=> vec(collect(prob.u0)))
    default_p = has_p ? Dict(params .=> vec(collect(prob.p))) : Dict()

    de = ODESystem(
        eqs, t, sts, params,
        defaults=merge(default_u0, default_p);
        kwargs...
    )

    de
end

_defvaridx(x, i, t) = variable(x, i, T=SymbolicUtils.FnType{Tuple,Real})
_defvar(x, t) = variable(x, T=SymbolicUtils.FnType{Tuple,Real})

function define_vars(u,t)
    _vars = [_defvaridx(:x, i, t)(t) for i in eachindex(u)]
end

function define_vars(u::Union{SLArray,LArray},t)
    _vars = [_defvar(x, t)(t) for x in LabelledArrays.symnames(typeof(u))]
end

function define_vars(u::Tuple,t)
    _vars = tuple((_defvaridx(:x, i, t)(ModelingToolkit.value(t)) for i in eachindex(u))...)
end

function define_vars(u::NamedTuple,t)
    _vars = NamedTuple(x=>_defvar(x, t)(ModelingToolkit.value(t)) for x in keys(u))
end

function define_params(p)
    [toparam(variable(:α, i)) for i in eachindex(p)]
end

function define_params(p::Union{SLArray,LArray})
    [toparam(variable(x)) for x in LabelledArrays.symnames(typeof(p))]
end

function define_params(p::Tuple)
    tuple((toparam(variable(:α, i)) for i in eachindex(p))...)
end

function define_params(p::NamedTuple)
    NamedTuple(x=>toparam(variable(x)) for x in keys(p))
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
    has_p = !(p isa Union{DiffEqBase.NullParameters,Nothing})

    _vars = define_vars(prob.u0,t)

    vars = prob.u0 isa Number ? _vars : ArrayInterface.restructure(prob.u0,_vars)
    params = if has_p
        _params = define_params(p)
        p isa Number ? _params[1] : (p isa Tuple || p isa NamedTuple ? _params : ArrayInterface.restructure(p,_params))
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

    de = SDESystem(deqs,neqs,t,Vector(vec(vars)),params; kwargs...)

    de
end


"""
$(TYPEDSIGNATURES)

Generate `OptimizationSystem`, dependent variables, and parameters from an `OptimizationProblem`.
"""
function modelingtoolkitize(prob::DiffEqBase.OptimizationProblem; kwargs...)

    if prob.p isa Tuple || prob.p isa NamedTuple
        p = [x for x in prob.p]
    else
        p = prob.p
    end

    vars = reshape([variable(:x, i) for i in eachindex(prob.u0)],size(prob.u0))
    params = p isa DiffEqBase.NullParameters ? [] :
        reshape([variable(:α, i) for i in eachindex(p)],size(Array(p)))

    eqs = prob.f(vars, params)
    de = OptimizationSystem(eqs,vec(vars),vec(params); kwargs...)
    de
end
