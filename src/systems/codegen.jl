"""
    $(TYPEDSIGNATURES)

Generate the RHS function for the `equations` of a `System`.

# Arguments

# Keyword Arguments

"""
function generate_rhs(sys::System, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true); implicit_dae = false,
        scalar = false, kwargs...)
    eqs = equations(sys)
    obs = observed(sys)
    u = dvs
    p = reorder_parameters(sys, ps)
    t = get_iv(sys)
    ddvs = nothing
    extra_assignments = Assignment[]

    # used for DAEProblem and ImplicitDiscreteProblem
    if implicit_dae
        if is_discrete_system(sys)
            # ImplicitDiscrete case
            D = Shift(t, 1)
            rhss = map(eqs) do eq
                # Algebraic equations get shifted forward 1, to match with differential
                # equations
                _iszero(eq.lhs) ? distribute_shift(D(eq.rhs)) : (eq.rhs - eq.lhs)
            end
            # Handle observables in algebraic equations, since they are shifted
            shifted_obs = Equation[distribute_shift(D(eq)) for eq in obs]
            obsidxs = observed_equations_used_by(sys, rhss; obs = shifted_obs)
            extra_assignments = [Assignment(shifted_obs[i].lhs, shifted_obs[i].rhs)
                                 for i in obsidxs]
        else
            D = Differential(t)
            rhss = [_iszero(eq.lhs) ? eq.rhs : eq.rhs - eq.lhs for eq in eqs]
        end
        ddvs = map(D, dvs)
    else
        check_operator_variables(eqs, Differential)
        check_lhs(eqs, Differential, Set(dvs))
        rhss = [eq.rhs for eq in eqs]
    end

    if !isempty(assertions(sys))
        rhss[end] += unwrap(get_assertions_expr(sys))
    end

    # TODO: add an optional check on the ordering of observed equations
    if scalar
        rhss = only(rhss)
        u = only(u)
    end

    args = (u, p...)
    p_start = 2
    if t !== nothing
        args = (args..., t)
    end
    if implicit_dae
        args = (ddvs, args...)
        p_start += 1
    end

    build_function_wrapper(sys, rhss, args...; p_start, extra_assignments, kwargs...)
end

function calculate_jacobian(sys::System;
        sparse = false, simplify = false, dvs = unknowns(sys))
    obs = Dict(eq.lhs => eq.rhs for eq in observed(sys))
    rhs = map(eq -> fixpoint_sub(eq.rhs, obs), equations(sys))

    if sparse
        jac = sparsejacobian(rhs, dvs, simplify)
        W_s = W_sparsity(sys)
        (Is, Js, Vs) = findnz(W_s)
        # Add nonzeros of W as non-structural zeros of the Jacobian (to ensure equal
        # results for oop and iip Jacobian)
        for (i, j) in zip(Is, Js)
            iszero(jac[i, j]) && begin
                jac[i, j] = 1
                jac[i, j] = 0
            end
        end
    else
        jac = jacobian(rhs, dvs; simplify)
    end

    return jac
end

