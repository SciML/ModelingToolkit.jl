"""
    $(TYPEDSIGNATURES)

Convert an `OptimizationProblem` to a `ModelingToolkit.System`.

# Keyword arguments

- `u_names`: An array of names of the same size as `prob.u0` to use as the names of the
  unknowns of the system. The names should be given as `Symbol`s.
- `p_names`: A collection of names to use for parameters of the system. The collection
  should have keys corresponding to indexes of `prob.p`. For example, if `prob.p` is an
  associative container like `NamedTuple`, then `p_names` should map keys of `prob.p` to
  the name that the corresponding parameter should have in the returned system. The names
  should be given as `Symbol`s.

All other keyword arguments are forwarded to the created `System`.
"""
function modelingtoolkitize(
        prob::OptimizationProblem; u_names = nothing, p_names = nothing,
        kwargs...)
    num_cons = isnothing(prob.lcons) ? 0 : length(prob.lcons)
    p = prob.p
    has_p = !(p isa Union{DiffEqBase.NullParameters, Nothing})

    vars = construct_vars(prob, nothing, u_names)
    if prob.ub !== nothing # lb is also !== nothing
        vars = map(vars, prob.lb, prob.ub) do sym, lb, ub
            if iszero(lb) && iszero(ub) || isinf(lb) && lb < 0 && isinf(ub) && ub > 0
                sym
            else
                Symbolics.setmetadata(sym, VariableBounds, (lb, ub))
            end
        end
    end
    params = construct_params(prob, nothing, p_names)

    objective = prob.f(vars, params)

    if prob.f.cons === nothing
        cons = []
    else
        if DiffEqBase.isinplace(prob.f)
            lhs = Array{Num}(undef, num_cons)
            prob.f.cons(lhs, vars, params)
        else
            lhs = prob.f.cons(vars, params)
        end
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
    end

    # turn `params` into a list of symbolic variables as opposed to
    # a parameter object containing symbolic variables.
    _params = params
    params = to_paramvec(params)

    defaults = defaults_from_u0_p(prob, vars, _params, params)
    # In case initials crept in, specifically from when we constructed parameters
    # using prob.f.sys
    filter!(x -> !iscall(x) || !(operation(x) isa Initial), params)
    filter!(x -> !iscall(x[1]) || !(operation(x[1]) isa Initial), defaults)

    sts = vec(collect(vars))
    sys = OptimizationSystem(objective, sts, params;
        defaults,
        constraints = cons,
        name = gensym(:MTKizedOpt),
        kwargs...)
end
