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



"""
$(TYPEDSIGNATURES)

Generate `SDESystem`, dependent variables, and parameters from an `SDEProblem`.
Choose correction_factor=-1//2 (1//2) to converte Ito -> Stratonovich (Stratonovich->Ito).
The default correction_factor is `nothing`.
"""
function modelingtoolkitize(prob::DiffEqBase.SDEProblem; correction_factor=nothing)
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

    if correction_factor!=nothing
        # use the general interface
        if DiffEqBase.is_diagonal_noise(prob)
            eqs = vcat([rhs[i] ~ neqs[i] for i in eachindex(prob.u0)]...)
            de = ODESystem(eqs,t,vec(vars),vec(params))

            jac = calculate_jacobian(de, sparse=false, simplify=true)
            ∇σσ′ = simplify.(jac*neqs)

            deqs = vcat([rhs[i] ~ lhs[i] + correction_factor*∇σσ′[i] for i in eachindex(prob.u0)]...)
        else
            dimstate, m = size(prob.noise_rate_prototype)
            eqs = vcat([rhs[i] ~ neqs[i] for i in eachindex(prob.u0)]...)
            de = ODESystem(eqs,t,vec(vars),vec(params))

            jac = calculate_jacobian(de, sparse=false, simplify=true)
            ∇σσ′ = simplify.(jac*neqs[:,1])
            for k = 2:m
                eqs = vcat([rhs[i] ~ neqs[Int(i+(k-1)*dimstate)] for i in eachindex(prob.u0)]...)
                de = ODESystem(eqs,t,vec(vars),vec(params))

                jac = calculate_jacobian(de, sparse=false, simplify=true)
                ∇σσ′ = ∇σσ′ + simplify.(jac*neqs[:,k])
            end

            deqs = vcat([rhs[i] ~ lhs[i] + correction_factor*∇σσ′[i] for i in eachindex(prob.u0)]...)
        end
    else
        deqs = vcat([rhs[i] ~ lhs[i] for i in eachindex(prob.u0)]...)
    end

    de = SDESystem(deqs,neqs,t,vec(vars),vec(params))

    de
end
