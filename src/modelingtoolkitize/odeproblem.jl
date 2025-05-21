"""
    $(TYPEDSIGNATURES)

Convert an `ODEProblem` to a `ModelingToolkit.System`.

# Keyword arguments

- `u_names`: An array of names of the same size as `prob.u0` to use as the names of the
  unknowns of the system. The names should be given as `Symbol`s.
- `p_names`: A collection of names to use for parameters of the system. The collection
  should have keys corresponding to indexes of `prob.p`. For example, if `prob.p` is an
  associative container like `NamedTuple`, then `p_names` should map keys of `prob.p` to
  the name that the corresponding parameter should have in the returned system. The names
  should be given as `Symbol`s.
- INTERNAL `return_symbolic_u0_p`: Also return the symbolic state and parameter objects.

All other keyword arguments are forwarded to the created `System`.
"""
function modelingtoolkitize(prob::ODEProblem; u_names = nothing, p_names = nothing,
        return_symbolic_u0_p = false, kwargs...)
    if prob.f isa DiffEqBase.AbstractParameterizedFunction
        return prob.f.sys
    end

    t = t_nounits
    p = prob.p
    has_p = !(p isa Union{DiffEqBase.NullParameters, Nothing})

    vars = construct_vars(prob, t, u_names)
    params = construct_params(prob, t, p_names)

    lhs = lhs_from_mass_matrix(D_nounits, prob.f.mass_matrix, vars)
    rhs = trace_rhs(prob, vars, params, t)
    eqs = vcat([lhs[i] ~ rhs[i] for i in eachindex(prob.u0)]...)

    sts = vec(collect(vars))

    # turn `params` into a list of symbolic variables as opposed to
    # a parameter object containing symbolic variables.
    _params = params
    params = to_paramvec(params)

    defaults = defaults_from_u0_p(prob, vars, _params, params)
    # In case initials crept in, specifically from when we constructed parameters
    # using prob.f.sys
    filter!(x -> !iscall(x) || !(operation(x) isa Initial), params)
    filter!(x -> !iscall(x[1]) || !(operation(x[1]) isa Initial), defaults)

    sys = System(eqs, t, sts, params;
        defaults,
        name = gensym(:MTKizedODE),
        kwargs...)

    if return_symbolic_u0_p
        return sys, vars, _params
    else
        return sys
    end
end
