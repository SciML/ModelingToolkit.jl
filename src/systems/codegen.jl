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

function calculate_tgrad(sys::System; simplify = false)
    # We need to remove explicit time dependence on the unknown because when we
    # have `u(t) * t` we want to have the tgrad to be `u(t)` instead of `u'(t) *
    # t + u(t)`.
    rhs = [detime_dvs(eq.rhs) for eq in full_equations(sys)]
    iv = get_iv(sys)
    xs = unknowns(sys)
    rule = Dict(map((x, xt) -> xt => x, detime_dvs.(xs), xs))
    rhs = substitute.(rhs, Ref(rule))
    tgrad = [expand_derivatives(Differential(iv)(r), simplify) for r in rhs]
    reverse_rule = Dict(map((x, xt) -> x => xt, detime_dvs.(xs), xs))
    tgrad = Num.(substitute.(tgrad, Ref(reverse_rule)))
    return tgrad
end

function calculate_jacobian(sys::System;
        sparse = false, simplify = false, dvs = unknowns(sys))
    obs = Dict(eq.lhs => eq.rhs for eq in observed(sys))
    rhs = map(eq -> fixpoint_sub(eq.rhs, obs), equations(sys))

    if sparse
        jac = sparsejacobian(rhs, dvs; simplify)
        if get_iv(sys) !== nothing
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
        end
    else
        jac = jacobian(rhs, dvs; simplify)
    end

    return jac
end

function generate_jacobian(sys::System, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true);
        simplify = false, sparse = false, kwargs...)
    jac = calculate_jacobian(sys; simplify, sparse, dvs)
    p = reorder_parameters(sys, ps)
    t = get_iv(sys)
    if t !== nothing
        wrap_code = sparse ? assert_jac_length_header(sys) : (identity, identity)
    end
    return build_function_wrapper(sys, jac, dvs, p..., t; wrap_code, kwargs...)
end

function assert_jac_length_header(sys)
    W = W_sparsity(sys)
    identity,
    function add_header(expr)
        Func(expr.args, [], expr.body,
            [:(@assert $(SymbolicUtils.Code.toexpr(term(findnz, expr.args[1])))[1:2] ==
                       $(findnz(W)[1:2]))])
    end
end

function generate_tgrad(
        sys::System, dvs = unknowns(sys), ps = parameters(
            sys; initial_parameters = true);
        simplify = false, kwargs...)
    tgrad = calculate_tgrad(sys, simplify = simplify)
    p = reorder_parameters(sys, ps)
    return build_function_wrapper(sys, tgrad,
        dvs,
        p...,
        get_iv(sys);
        kwargs...)
end

const W_GAMMA = only(@variables ˍ₋gamma)

function generate_W(sys::System, γ = 1.0, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true);
        simplify = false, sparse = false, kwargs...)
    M = calculate_massmatrix(sys; simplify)
    if sparse
        M = SparseArrays.sparse(M)
    end
    J = calculate_jacobian(sys; simplify, sparse, dvs)
    W = W_GAMMA * M + J
    t = get_iv(sys)
    if t !== nothing
        wrap_code = sparse ? assert_jac_length_header(sys) : (identity, identity)
    end

    p = reorder_parameters(sys, ps)
    return build_function_wrapper(sys, W, dvs, p..., W_GAMMA, t; wrap_code,
        p_end = 1 + length(p), kwargs...)
end

function generate_dae_jacobian(sys::System, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true); simplify = false, sparse = false,
        kwargs...)
    jac_u = calculate_jacobian(sys; simplify = simplify, sparse = sparse)
    t = get_iv(sys)
    derivatives = Differential(t).(unknowns(sys))
    jac_du = calculate_jacobian(sys; simplify = simplify, sparse = sparse,
        dvs = derivatives)
    dvs = unknowns(sys)
    jac = W_GAMMA * jac_du + jac_u
    p = reorder_parameters(sys, ps)
    return build_function_wrapper(sys, jac, derivatives, dvs, p..., W_GAMMA, t;
        p_start = 3, p_end = 2 + length(p), kwargs...)
end

function calculate_massmatrix(sys::System; simplify = false)
    eqs = [eq for eq in equations(sys)]
    M = zeros(length(eqs), length(eqs))
    for (i, eq) in enumerate(eqs)
        if iscall(eq.lhs) && operation(eq.lhs) isa Differential
            st = var_from_nested_derivative(eq.lhs)[1]
            j = variable_index(sys, st)
            M[i, j] = 1
        else
            _iszero(eq.lhs) ||
                error("Only semi-explicit constant mass matrices are currently supported. Faulty equation: $eq.")
        end
    end
    M = simplify ? simplify.(M) : M
    # M should only contain concrete numbers
    M == I ? I : M
end

function jacobian_sparsity(sys::System)
    sparsity = torn_system_jacobian_sparsity(sys)
    sparsity === nothing || return sparsity

    Symbolics.jacobian_sparsity([eq.rhs for eq in full_equations(sys)],
        [dv for dv in unknowns(sys)])
end

function jacobian_dae_sparsity(sys::System)
    J1 = jacobian_sparsity([eq.rhs for eq in full_equations(sys)],
        [dv for dv in unknowns(sys)])
    derivatives = Differential(get_iv(sys)).(unknowns(sys))
    J2 = jacobian_sparsity([eq.rhs for eq in full_equations(sys)],
        [dv for dv in derivatives])
    J1 + J2
end

function W_sparsity(sys::System)
    jac_sparsity = jacobian_sparsity(sys)
    (n, n) = size(jac_sparsity)
    M = calculate_massmatrix(sys)
    M_sparsity = M isa UniformScaling ? sparse(I(n)) :
                 SparseMatrixCSC{Bool, Int64}((!iszero).(M))
    jac_sparsity .| M_sparsity
end

function isautonomous(sys::System)
    tgrad = calculate_tgrad(sys; simplify = true)
    all(iszero, tgrad)
end
