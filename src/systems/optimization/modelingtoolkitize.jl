"""
$(TYPEDSIGNATURES)

Generate `OptimizationSystem`, dependent variables, and parameters from an `OptimizationProblem`.
"""
function modelingtoolkitize(prob::DiffEqBase.OptimizationProblem; kwargs...)
    num_cons = isnothing(prob.lcons) ? 0 : length(prob.lcons)
    if prob.p isa Tuple || prob.p isa NamedTuple
        p = [x for x in prob.p]
    else
        p = prob.p
    end

    vars = ArrayInterface.restructure(prob.u0,
        [variable(:x, i) for i in eachindex(prob.u0)])
    params = p isa DiffEqBase.NullParameters ? [] :
             ArrayInterface.restructure(p, [variable(:α, i) for i in eachindex(p)])

    eqs = prob.f(vars, params)

    if DiffEqBase.isinplace(prob) && !isnothing(prob.f.cons)
        lhs = Array{Num}(undef, num_cons)
        prob.f.cons(lhs, vars, params)
        cons = Union{Equation, Inequality}[]

        if !isnothing(prob.lcons)
            for i in 1:num_cons
                if !isinf(prob.lcons[i])
                    if prob.lcons[i] != prob.ucons[i]
                        push!(cons, prob.lcons[i] ≲ lhs[i])
                    else
                        push!(cons, lhs[i] ~ prob.ucons[i])
                    end
                end
            end
        end

        if !isnothing(prob.ucons)
            for i in 1:num_cons
                if !isinf(prob.ucons[i]) && prob.lcons[i] != prob.ucons[i]
                    push!(cons, lhs[i] ≲ prob.ucons[i])
                end
            end
        end
        if (isnothing(prob.lcons) || all(isinf, prob.lcons)) &&
           (isnothing(prob.ucons) || all(isinf, prob.ucons))
            throw(ArgumentError("Constraints passed have no proper bounds defined.
            Ensure you pass equal bounds (the scalar that the constraint should evaluate to) for equality constraints
            or pass the lower and upper bounds for inequality constraints."))
        end
    elseif !isnothing(prob.f.cons)
        cons = prob.f.cons(vars, params)
    else
        cons = []
    end

    de = OptimizationSystem(eqs, vec(vars), vec(toparam.(params));
        name = gensym(:MTKizedOpt),
        constraints = cons,
        kwargs...)
    de
end
