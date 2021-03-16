"""
$(TYPEDSIGNATURES)

Generate `ODESystem`, dependent variables, and parameters from an `ODEProblem`.
"""
function modelingtoolkitize(prob::DiffEqBase.ODEProblem)
    prob.f isa DiffEqBase.AbstractParameterizedFunction &&
                            return prob.f.sys
    @parameters t

    if prob.p isa Tuple || prob.p isa NamedTuple
        p = [x for x in prob.p]
    else
        p = prob.p
    end

    var(x, i) = Num(Sym{FnType{Tuple{symtype(t)}, Real}}(nameof(Variable(x, i))))
    vars = ArrayInterface.restructure(prob.u0,[var(:x, i)(ModelingToolkit.value(t)) for i in eachindex(prob.u0)])
    params = p isa DiffEqBase.NullParameters ? [] :
             reshape([Num(Sym{Real}(nameof(Variable(:α, i)))) for i in eachindex(p)],size(p))
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
        rhs = similar(vars, Num)
        prob.f(rhs, vars, params, t)
    else
        rhs = prob.f(vars, params, t)
    end

    eqs = vcat([lhs[i] ~ rhs[i] for i in eachindex(prob.u0)]...)

    sts = vec(collect(vars))
    params = if ndims(params) == 0
        [params[1]]
    else
        vec(collect(params))
    end

    de = ODESystem(
        eqs, t, sts, params,
        defaults=merge(Dict(sts .=> vec(collect(prob.u0))), Dict(params .=> vec(collect(prob.p)))),
    )

    de
end



"""
$(TYPEDSIGNATURES)

Generate `SDESystem`, dependent variables, and parameters from an `SDEProblem`.
"""
function modelingtoolkitize(prob::DiffEqBase.SDEProblem)
    prob.f isa DiffEqBase.AbstractParameterizedFunction &&
                            return (prob.f.sys, prob.f.sys.states, prob.f.sys.ps)
    @parameters t
    if prob.p isa Tuple || prob.p isa NamedTuple
        p = [x for x in prob.p]
    else
        p = prob.p
    end
    var(x, i) = Num(Sym{FnType{Tuple{symtype(t)}, Real}}(nameof(Variable(x, i))))
    vars = ArrayInterface.restructure(prob.u0,[var(:x, i)(ModelingToolkit.value(t)) for i in eachindex(prob.u0)])
    params = p isa DiffEqBase.NullParameters ? [] :
             reshape([Num(Sym{Real}(nameof(Variable(:α, i)))) for i in eachindex(p)],size(p))

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

    de = SDESystem(deqs,neqs,t,Vector(vec(vars)),params)

    de
end


"""
$(TYPEDSIGNATURES)

Generate `OptimizationSystem`, dependent variables, and parameters from an `OptimizationProblem`.
"""
function modelingtoolkitize(prob::DiffEqBase.OptimizationProblem)

    if prob.p isa Tuple || prob.p isa NamedTuple
        p = [x for x in prob.p]
    else
        p = prob.p
    end

    vars = reshape([Num(Sym{Real}(nameof(Variable(:x, i)))) for i in eachindex(prob.u0)],size(prob.u0))
    params = p isa DiffEqBase.NullParameters ? [] :
             reshape([Num(Sym{Real}(nameof(Variable(:α, i)))) for i in eachindex(p)],size(Array(p)))

    eqs = prob.f(vars, params)
    de = OptimizationSystem(eqs,vec(vars),vec(params))
    de
end
