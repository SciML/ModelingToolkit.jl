"""
$(TYPEDSIGNATURES)

Generate `NonlinearSystem`, dependent variables, and parameters from an `NonlinearProblem`.
"""
function modelingtoolkitize(prob::NonlinearProblem; kwargs...)
    p = prob.p
    has_p = !(p isa Union{DiffEqBase.NullParameters, Nothing})

    _vars = reshape([variable(:x, i) for i in eachindex(prob.u0)], size(prob.u0))

    vars = prob.u0 isa Number ? _vars : ArrayInterface.restructure(prob.u0, _vars)
    params = if has_p
        _params = define_params(p)
        p isa Number ? _params[1] :
        (p isa Tuple || p isa NamedTuple ? _params :
         ArrayInterface.restructure(p, _params))
    else
        []
    end

    if DiffEqBase.isinplace(prob)
        rhs = ArrayInterface.restructure(prob.u0, similar(vars, Num))
        prob.f(rhs, vars, params)
    else
        rhs = prob.f(vars, params)
    end
    out_def = prob.f(prob.u0, prob.p)
    eqs = vcat([0.0 ~ rhs[i] for i in 1:length(out_def)]...)

    sts = vec(collect(vars))

    params = if params isa Number || (params isa Array && ndims(params) == 0)
        [params[1]]
    else
        vec(collect(params))
    end
    default_u0 = Dict(sts .=> vec(collect(prob.u0)))
    default_p = has_p ? Dict(params .=> vec(collect(prob.p))) : Dict()

    de = NonlinearSystem(eqs, sts, params,
                         defaults = merge(default_u0, default_p);
                         name = gensym(:MTKizedNonlinProb),
                         kwargs...)

    de
end
