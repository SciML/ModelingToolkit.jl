"""
    $(TYPEDSIGNATURES)

Convert an `SDEProblem` to a `ModelingToolkit.System`.

# Keyword arguments

- `u_names`: an array of names of the same size as `prob.u0` to use as the names of the
  unknowns of the system. The names should be given as `Symbol`s.
- `p_names`: a collection of names to use for parameters of the system. The collection
  should have keys corresponding to indexes of `prob.p`. For example, if `prob.p` is an
  associative container like `NamedTuple`, then `p_names` should map keys of `prob.p` to
  the name that the corresponding parameter should have in the returned system. The names
  should be given as `Symbol`s.

All other keyword arguments are forwarded to the created `System`.
"""
function modelingtoolkitize(
        prob::SDEProblem; u_names = nothing, p_names = nothing, kwargs...)
    if prob.f isa DiffEqBase.AbstractParameterizedFunction
        return prob.f.sys
    end

    # just create the equivalent ODEProblem, `modelingtoolkitize` that
    # and add on the noise
    odefn = ODEFunction{SciMLBase.isinplace(prob)}(
        prob.f.f; mass_matrix = prob.f.mass_matrix, sys = prob.f.sys)
    odeprob = ODEProblem(odefn, prob.u0, prob.tspan, prob.p)
    sys, vars,
    params = modelingtoolkitize(
        odeprob; u_names, p_names, return_symbolic_u0_p = true,
        name = gensym(:MTKizedSDE), kwargs...)
    t = get_iv(sys)

    if SciMLBase.isinplace(prob)
        if SciMLBase.is_diagonal_noise(prob)
            neqs = similar(vars, Any)
            prob.g(neqs, vars, params, t)
        else
            neqs = similar(prob.noise_rate_prototype, Any)
            prob.g(neqs, vars, params, t)
        end
    else
        if SciMLBase.is_diagonal_noise(prob)
            neqs = prob.g(vars, params, t)
        else
            neqs = prob.g(vars, params, t)
        end
    end

    @set! sys.noise_eqs = neqs

    return sys
end
