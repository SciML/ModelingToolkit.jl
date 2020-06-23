"""
$(TYPEDSIGNATURES)

Generate `ODESystem`, dependent variables, and parameters from an `ODEProblem`.
"""
function modelingtoolkitize(prob::DiffEqBase.ODEProblem)
    prob.f isa DiffEqBase.AbstractParameterizedFunction &&
                            return (prob.f.sys, prob.f.sys.states, prob.f.sys.ps)
    @parameters t
    vars = reshape([Variable(:x, i)(t) for i in eachindex(prob.u0)],size(prob.u0))
    params = prob.p isa DiffEqBase.NullParameters ? [] :
             reshape([Variable(:α,i)() for i in eachindex(prob.p)],size(prob.p))
    @derivatives D'~t

    rhs = [D(var) for var in vars]

    if DiffEqBase.isinplace(prob)
        lhs = similar(vars, Any)
        prob.f(lhs, vars, params, t)
    else
        lhs = prob.f(vars, params, t)
    end

    eqs = vcat([rhs[i] ~ lhs[i] for i in eachindex(prob.u0)]...)
    de = ODESystem(eqs,t,vec(vars),vec(params))

    de
end
