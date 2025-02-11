"""
$(TYPEDSIGNATURES)

Generate `NonlinearSystem`, dependent variables, and parameters from an `NonlinearProblem`.
"""
function modelingtoolkitize(
        prob::Union{NonlinearProblem, NonlinearLeastSquaresProblem};
        u_names = nothing, p_names = nothing, kwargs...)
    p = prob.p
    has_p = !(p isa Union{DiffEqBase.NullParameters, Nothing})

    if u_names !== nothing
        varnames_length_check(prob.u0, u_names; is_unknowns = true)
        _vars = [variable(name) for name in u_names]
    elseif SciMLBase.has_sys(prob.f)
        varnames = getname.(variable_symbols(prob.f.sys))
        varidxs = variable_index.((prob.f.sys,), varnames)
        invpermute!(varnames, varidxs)
        _vars = [variable(name) for name in varnames]
    else
        _vars = [variable(:x, i) for i in eachindex(prob.u0)]
    end
    _vars = reshape(_vars, size(prob.u0))

    vars = prob.u0 isa Number ? _vars : ArrayInterface.restructure(prob.u0, _vars)
    params = if has_p
        if p_names === nothing && SciMLBase.has_sys(prob.f)
            p_names = Dict(parameter_index(prob.f.sys, sym) => sym
            for sym in parameter_symbols(prob.f.sys))
        end
        _params = define_params(p, p_names)
        p isa Number ? _params[1] :
        (p isa Tuple || p isa NamedTuple || p isa AbstractDict || p isa MTKParameters ?
         _params :
         ArrayInterface.restructure(p, _params))
    else
        []
    end

    if DiffEqBase.isinplace(prob)
        if prob isa NonlinearLeastSquaresProblem
            rhs = ArrayInterface.restructure(
                prob.f.resid_prototype, similar(prob.f.resid_prototype, Num))
            prob.f(rhs, vars, params)
            eqs = vcat([0.0 ~ rhs[i] for i in 1:length(prob.f.resid_prototype)]...)
        else
            rhs = ArrayInterface.restructure(prob.u0, similar(vars, Num))
            prob.f(rhs, vars, p isa MTKParameters ? (params,) : params)
            eqs = vcat([0.0 ~ rhs[i] for i in 1:length(rhs)]...)
        end

    else
        rhs = prob.f(vars, params)
        out_def = prob.f(prob.u0, prob.p)
        eqs = vcat([0.0 ~ rhs[i] for i in 1:length(out_def)]...)
    end

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
        elseif prob.p isa MTKParameters
            Dict(params .=> reduce(vcat, prob.p))
        else
            Dict(params .=> vec(collect(prob.p)))
        end
    else
        Dict()
    end

    de = NonlinearSystem(eqs, sts, params,
        defaults = merge(default_u0, default_p);
        name = gensym(:MTKizedNonlinProb),
        kwargs...)

    de
end
