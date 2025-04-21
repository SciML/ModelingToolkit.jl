"""
    $(TYPEDSIGNATURES)

Generate the RHS function for the `equations` of a `System`.

# Arguments

# Keyword Arguments

"""
function generate_rhs(sys::System, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true); implicit_dae = false,
        scalar = false, expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, override_discrete = false,
        kwargs...)
    eqs = equations(sys)
    obs = observed(sys)
    u = dvs
    p = reorder_parameters(sys, ps)
    t = get_iv(sys)
    ddvs = nothing
    extra_assignments = Assignment[]

    # used for DAEProblem and ImplicitDiscreteProblem
    if implicit_dae
        if override_discrete || is_discrete_system(sys)
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
        if !override_discrete && !is_discrete_system(sys)
            check_operator_variables(eqs, Differential)
            check_lhs(eqs, Differential, Set(dvs))
        end
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

    res = build_function_wrapper(sys, rhss, args...; p_start, extra_assignments,
        expression = Val{true}, expression_module = eval_module, kwargs...)
    nargs = length(args) - length(p) + 1
    if is_dde(sys)
        p_start += 1
        nargs += 1
    end
    return maybe_compile_function(
        expression, wrap_gfw, (p_start, nargs, is_split(sys)),
        res; eval_expression, eval_module)
end

function generate_diffusion_function(sys::System, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true); expression = Val{true},
        wrap_gfw = Val{false}, eval_expression = false,
        eval_module = @__MODULE__, kwargs...)
    eqs = get_noise_eqs(sys)
    if ndims(eqs) == 2 && size(eqs, 2) == 1
        # scalar noise
        eqs = vec(eqs)
    end
    p = reorder_parameters(sys, ps)
    res = build_function_wrapper(sys, eqs, dvs, p..., get_iv(sys); kwargs...)
    if expression == Val{true}
        return res
    end
    f_oop, f_iip = eval_or_rgf.(res; eval_expression, eval_module)
    p_start = 2
    nargs = 3
    if is_dde(sys)
        p_start += 1
        nargs += 1
    end
    return maybe_compile_function(
        expression, wrap_gfw, (p_start, nargs, is_split(sys)), res; eval_expression, eval_module)
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
    rhs = map(eq -> fixpoint_sub(eq.rhs - eq.lhs, obs), equations(sys))

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

function generate_jacobian(sys::System;
        simplify = false, sparse = false, eval_expression = false,
        eval_module = @__MODULE__, expression = Val{true}, wrap_gfw = Val{false},
        kwargs...)
    dvs = unknowns(sys)
    jac = calculate_jacobian(sys; simplify, sparse, dvs)
    p = reorder_parameters(sys)
    t = get_iv(sys)
    if t === nothing
        wrap_code = (identity, identity)
    else
        wrap_code = sparse ? assert_jac_length_header(sys) : (identity, identity)
    end
    args = (dvs, p...)
    nargs = 2
    if is_time_dependent(sys)
        args = (args..., t)
        nargs = 3
    end
    res = build_function_wrapper(sys, jac, args...; wrap_code, expression = Val{true},
        expression_module = eval_module, kwargs...)
    return maybe_compile_function(
        expression, wrap_gfw, (2, nargs, is_split(sys)), res; eval_expression, eval_module)
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
        simplify = false, eval_expression = false, eval_module = @__MODULE__,
        expression = Val{true}, wrap_gfw = Val{false}, kwargs...)
    tgrad = calculate_tgrad(sys, simplify = simplify)
    p = reorder_parameters(sys, ps)
    res = build_function_wrapper(sys, tgrad,
        dvs,
        p...,
        get_iv(sys);
        expression = Val{true},
        expression_module = eval_module,
        kwargs...)

    return maybe_compile_function(
        expression, wrap_gfw, (2, 3, is_split(sys)), res; eval_expression, eval_module)
end

const W_GAMMA = only(@variables ˍ₋gamma)

function generate_W(sys::System, γ = 1.0, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true);
        simplify = false, sparse = false, expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, kwargs...)
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
    res = build_function_wrapper(sys, W, dvs, p..., W_GAMMA, t; wrap_code,
        p_end = 1 + length(p), kwargs...)
    return maybe_compile_function(
        expression, wrap_gfw, (2, 4, is_split(sys)), res; eval_expression, eval_module)
end

function generate_dae_jacobian(sys::System, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true); simplify = false, sparse = false,
        expression = Val{true}, wrap_gfw = Val{false}, eval_expression = false,
        eval_module = @__MODULE__, kwargs...)
    jac_u = calculate_jacobian(sys; simplify = simplify, sparse = sparse)
    t = get_iv(sys)
    derivatives = Differential(t).(unknowns(sys))
    jac_du = calculate_jacobian(sys; simplify = simplify, sparse = sparse,
        dvs = derivatives)
    dvs = unknowns(sys)
    jac = W_GAMMA * jac_du + jac_u
    p = reorder_parameters(sys, ps)
    res = build_function_wrapper(sys, jac, derivatives, dvs, p..., W_GAMMA, t;
        p_start = 3, p_end = 2 + length(p), kwargs...)
    return maybe_compile_function(
        expression, wrap_gfw, (3, 5, is_split(sys)), res; eval_expression, eval_module)
end

function generate_history(sys::System, u0; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, kwargs...)
    p = reorder_parameters(sys)
    res = build_function_wrapper(sys, u0, p..., get_iv(sys); expression = Val{true},
        expression_module = eval_module, p_start = 1, p_end = length(p),
        similarto = typeof(u0), wrap_delays = false, kwargs...)
    return maybe_compile_function(
        expression, wrap_gfw, (1, 2, is_split(sys)), res; eval_expression, eval_module)
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
    if isdiag(M)
        M = Diagonal(M)
    end
    # M should only contain concrete numbers
    M == I ? I : M
end

function concrete_massmatrix(M; sparse = false, u0 = nothing)
    if sparse && !(u0 === nothing || M === I)
        SparseArrays.sparse(M)
    elseif u0 === nothing || M === I
        M
    elseif M isa Diagonal
        Diagonal(ArrayInterface.restructure(u0, diag(M)))
    else
        ArrayInterface.restructure(u0 .* u0', M)
    end
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

function calculate_W_prototype(W_sparsity; u0 = nothing, sparse = false)
    sparse || return nothing
    uElType = u0 === nothing ? Float64 : eltype(u0)
    return similar(W_sparsity, uElType)
end

function isautonomous(sys::System)
    tgrad = calculate_tgrad(sys; simplify = true)
    all(iszero, tgrad)
end

function get_bv_solution_symbol(ns)
    only(@variables BV_SOLUTION(..)[1:ns])
end

function get_constraint_unknown_subs!(subs::Dict, cons::Vector, stidxmap::Dict, iv, sol)
    vs = vars(cons)
    for v in vs
        iscall(v) || continue
        op = operation(v)
        args = arguments(v)
        issym(op) && length(args) == 1 || continue
        newv = op(iv)
        haskey(stidxmap, newv) || continue
        subs[v] = sol(args[1])[stidxmap[newv]]
    end
end

function generate_boundary_conditions(sys::System, u0, u0_idxs, t0; expression = Val{true},
        wrap_gfw = Val{false}, eval_expression = false, eval_module = @__MODULE__,
        kwargs...)
    iv = get_iv(sys)
    sts = unknowns(sys)
    ps = parameters(sys)
    np = length(ps)
    ns = length(sts)
    stidxmap = Dict([v => i for (i, v) in enumerate(sts)])
    pidxmap = Dict([v => i for (i, v) in enumerate(ps)])

    # sol = get_bv_solution_symbol(ns)

    cons = [con.lhs - con.rhs for con in constraints(sys)]
    # conssubs = Dict()
    # get_constraint_unknown_subs!(conssubs, cons, stidxmap, iv, sol)
    # cons = map(x -> fast_substitute(x, conssubs), cons)

    init_conds = Any[]
    for i in u0_idxs
        expr = BVP_SOLUTION(t0)[i] - u0[i]
        push!(init_conds, expr)
    end

    exprs = vcat(init_conds, cons)
    _p = reorder_parameters(sys, ps)

    res = build_function_wrapper(sys, exprs, _p..., iv; output_type = Array,
        p_start = 1, histfn = (p, t) -> BVP_SOLUTION(t),
        histfn_symbolic = BVP_SOLUTION, wrap_delays = true, kwargs...)
    return maybe_compile_function(
        expression, wrap_gfw, (2, 3, is_split(sys)), res; eval_expression, eval_module)
end

function generate_cost(sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, kwargs...)
    obj = cost(sys)
    dvs = unknowns(sys)
    ps = reorder_parameters(sys)
    res = build_function_wrapper(sys, obj, dvs, ps...; expression = Val{true}, kwargs...)
    if expression == Val{true}
        return res
    end
    f_oop = eval_or_rgf(res; eval_expression, eval_module)
    return maybe_compile_function(
        expression, wrap_gfw, (2, 2, is_split(sys)), res; eval_expression, eval_module)
end

function calculate_cost_gradient(sys::System; simplify = false)
    obj = cost(sys)
    dvs = unknowns(sys)
    return Symbolics.gradient(obj, dvs; simplify)
end

function generate_cost_gradient(
        sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, simplify = false, kwargs...)
    obj = cost(sys)
    dvs = unknowns(sys)
    ps = reorder_parameters(sys)
    exprs = calculate_cost_gradient(sys; simplify)
    res = build_function_wrapper(sys, exprs, dvs, ps...; expression = Val{true}, kwargs...)
    return maybe_compile_function(
        expression, wrap_gfw, (2, 2, is_split(sys)), res; eval_expression, eval_module)
end

function calculate_cost_hessian(sys::System; sparse = false, simplify = false)
    obj = cost(sys)
    dvs = unknowns(sys)
    if sparse
        exprs = Symbolics.sparsehessian(obj, dvs; simplify)::AbstractSparseArray
        sparsity = similar(exprs, Float64)
    else
        exprs = Symbolics.hessian(obj, dvs; simplify)
    end
end

function cost_hessian_sparsity(sys::System)
    return similar(calculate_cost_hessian(sys; sparse = true), Float64)
end

function generate_cost_hessian(
        sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, simplify = false,
        sparse = false, return_sparsity = false, kwargs...)
    obj = cost(sys)
    dvs = unknowns(sys)
    ps = reorder_parameters(sys)
    sparsity = nothing
    exprs = calculate_cost_hessian(sys; sparse, simplify)
    if sparse
        sparsity = similar(exprs, Float64)
    end
    res = build_function_wrapper(sys, exprs, dvs, ps...; expression = Val{true}, kwargs...)
    fn = maybe_compile_function(
        expression, wrap_gfw, (2, 2, is_split(sys)), res; eval_expression, eval_module)

    return return_sparsity ? (fn, sparsity) : fn
end

function canonical_constraints(sys::System)
    return map(constraints(sys)) do cstr
        Symbolics.canonical_form(cstr).lhs
    end
end

function generate_cons(sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, kwargs...)
    cons = canonical_constraints(sys)
    dvs = unknowns(sys)
    ps = reorder_parameters(sys)
    res = build_function_wrapper(sys, cons, dvs, ps...; expression = Val{true}, kwargs...)
    return maybe_compile_function(
        expression, wrap_gfw, (2, 2, is_split(sys)), res; eval_expression, eval_module)
end

function generate_constraint_jacobian(
        sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, return_sparsity = false,
        simplify = false, sparse = false, kwargs...)
    cons = canonical_constraints(sys)
    dvs = unknowns(sys)
    ps = reorder_parameters(sys)
    sparsity = nothing
    if sparse
        jac = Symbolics.sparsejacobian(cons, dvs; simplify)::AbstractSparseArray
        sparsity = similar(jac, Float64)
    else
        jac = Symbolics.jacobian(cons, dvs; simplify)
    end
    res = build_function_wrapper(sys, jac, dvs, ps...; expression = Val{true}, kwargs...)
    fn = maybe_compile_function(
        expression, wrap_gfw, (2, 2, is_split(sys)), res; eval_expression, eval_module)
    return return_sparsity ? (fn, sparsity) : fn
end

function generate_constraint_hessian(
        sys::System; expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, return_sparsity = false,
        simplify = false, sparse = false, kwargs...)
    cons = canonical_constraints(sys)
    dvs = unknowns(sys)
    ps = reorder_parameters(sys)
    sparsity = nothing
    if sparse
        hess = map(cons) do cstr
            Symbolics.sparsehessian(cstr, dvs; simplify)::AbstractSparseArray
        end
        sparsity = similar.(hess, Float64)
    else
        hess = [Symbolics.hessian(cstr, dvs; simplify) for cstr in cons]
    end
    res = build_function_wrapper(sys, hess, dvs, ps...; expression = Val{true}, kwargs...)
    fn = maybe_compile_function(
        expression, wrap_gfw, (2, 2, is_split(sys)), res; eval_expression, eval_module)
    return return_sparsity ? (fn, sparsity) : fn
end

function calculate_control_jacobian(sys::AbstractSystem;
        sparse = false, simplify = false)
    rhs = [eq.rhs for eq in full_equations(sys)]
    ctrls = unbound_inputs(sys)

    if sparse
        jac = sparsejacobian(rhs, ctrls, simplify = simplify)
    else
        jac = jacobian(rhs, ctrls, simplify = simplify)
    end

    return jac
end

function generate_control_jacobian(sys::AbstractSystem, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true);
        expression = Val{true}, wrap_gfw = Val{false}, eval_expression = false,
        eval_module = @__MODULE__, simplify = false, sparse = false, kwargs...)
    jac = calculate_control_jacobian(sys; simplify = simplify, sparse = sparse)
    p = reorder_parameters(sys, ps)
    res = build_function_wrapper(sys, jac, dvs, p..., get_iv(sys); kwargs...)
    return maybe_compile_function(
        expression, wrap_gfw, (2, 3, is_split(sys)), res; eval_expression, eval_module)
end

function generate_rate_function(js::System, rate)
    p = reorder_parameters(js)
    build_function_wrapper(js, rate, unknowns(js), p...,
        get_iv(js),
        expression = Val{true})
end

function generate_affect_function(js::System, affect; kwargs...)
    compile_equational_affect(affect, js; checkvars = false, kwargs...)
end

function assemble_vrj(
        js, vrj, unknowntoid; eval_expression = false, eval_module = @__MODULE__)
    rate = eval_or_rgf(generate_rate_function(js, vrj.rate); eval_expression, eval_module)
    rate = GeneratedFunctionWrapper{(2, 3, is_split(js))}(rate, nothing)
    outputvars = (value(affect.lhs) for affect in vrj.affect!)
    outputidxs = [unknowntoid[var] for var in outputvars]
    affect = generate_affect_function(js, vrj.affect!; eval_expression, eval_module)
    VariableRateJump(rate, affect; save_positions = vrj.save_positions)
end

function assemble_crj(
        js, crj, unknowntoid; eval_expression = false, eval_module = @__MODULE__)
    rate = eval_or_rgf(generate_rate_function(js, crj.rate); eval_expression, eval_module)
    rate = GeneratedFunctionWrapper{(2, 3, is_split(js))}(rate, nothing)
    outputvars = (value(affect.lhs) for affect in crj.affect!)
    outputidxs = [unknowntoid[var] for var in outputvars]
    affect = generate_affect_function(js, crj.affect!; eval_expression, eval_module)
    ConstantRateJump(rate, affect)
end

# assemble a numeric MassActionJump from a MT symbolics MassActionJumps
function assemble_maj(majv::Vector{U}, unknowntoid, pmapper) where {U <: MassActionJump}
    rs = [numericrstoich(maj.reactant_stoch, unknowntoid) for maj in majv]
    ns = [numericnstoich(maj.net_stoch, unknowntoid) for maj in majv]
    MassActionJump(rs, ns; param_mapper = pmapper, nocopy = true)
end

function numericrstoich(mtrs::Vector{Pair{V, W}}, unknowntoid) where {V, W}
    rs = Vector{Pair{Int, W}}()
    for (wspec, stoich) in mtrs
        spec = value(wspec)
        if !iscall(spec) && _iszero(spec)
            push!(rs, 0 => stoich)
        else
            push!(rs, unknowntoid[spec] => stoich)
        end
    end
    sort!(rs)
    rs
end

function numericnstoich(mtrs::Vector{Pair{V, W}}, unknowntoid) where {V, W}
    ns = Vector{Pair{Int, W}}()
    for (wspec, stoich) in mtrs
        spec = value(wspec)
        !iscall(spec) && _iszero(spec) &&
            error("Net stoichiometry can not have a species labelled 0.")
        push!(ns, unknowntoid[spec] => stoich)
    end
    sort!(ns)
end

"""
    build_explicit_observed_function(sys, ts; kwargs...) -> Function(s)

Generates a function that computes the observed value(s) `ts` in the system `sys`, while making the assumption that there are no cycles in the equations.

## Arguments 
- `sys`: The system for which to generate the function
- `ts`: The symbolic observed values whose value should be computed

## Keywords
- `return_inplace = false`: If true and the observed value is a vector, then return both the in place and out of place methods.
- `expression = false`: Generates a Julia `Expr`` computing the observed value if `expression` is true
- `eval_expression = false`: If true and `expression = false`, evaluates the returned function in the module `eval_module`
- `output_type = Array` the type of the array generated by a out-of-place vector-valued function
- `param_only = false` if true, only allow the generated function to access system parameters
- `inputs = nothing` additinoal symbolic variables that should be provided to the generated function
- `checkbounds = true` checks bounds if true when destructuring parameters
- `op = Operator` sets the recursion terminator for the walk done by `vars` to identify the variables that appear in `ts`. See the documentation for `vars` for more detail.
- `throw = true` if true, throw an error when generating a function for `ts` that reference variables that do not exist.
- `mkarray`: only used if the output is an array (that is, `!isscalar(ts)`  and `ts` is not a tuple, in which case the result will always be a tuple). Called as `mkarray(ts, output_type)` where `ts` are the expressions to put in the array and `output_type` is the argument of the same name passed to build_explicit_observed_function.
- `cse = true`: Whether to use Common Subexpression Elimination (CSE) to generate a more efficient function.

## Returns

The return value will be either:
* a single function `f_oop` if the input is a scalar or if the input is a Vector but `return_inplace` is false
* the out of place and in-place functions `(f_ip, f_oop)` if `return_inplace` is true and the input is a `Vector`

The function(s) `f_oop` (and potentially `f_ip`) will be:
* `RuntimeGeneratedFunction`s by default,
* A Julia `Expr` if `expression` is true,
* A directly evaluated Julia function in the module `eval_module` if `eval_expression` is true and `expression` is false.

The signatures will be of the form `g(...)` with arguments:

- `output` for in-place functions
- `unknowns` if `param_only` is `false`
- `inputs` if `inputs` is an array of symbolic inputs that should be available in `ts` 
- `p...` unconditionally; note that in the case of `MTKParameters` more than one parameters argument may be present, so it must be splatted
- `t` if the system is time-dependent; for example systems of nonlinear equations will not have `t`

For example, a function `g(op, unknowns, p..., inputs, t)` will be the in-place function generated if `return_inplace` is true, `ts` is a vector, 
an array of inputs `inputs` is given, and `param_only` is false for a time-dependent system.
"""
function build_explicit_observed_function(sys, ts;
        inputs = nothing,
        disturbance_inputs = nothing,
        disturbance_argument = false,
        expression = false,
        eval_expression = false,
        eval_module = @__MODULE__,
        output_type = Array,
        checkbounds = true,
        ps = parameters(sys; initial_parameters = true),
        return_inplace = false,
        param_only = false,
        op = Operator,
        throw = true,
        cse = true,
        mkarray = nothing)
    # TODO: cleanup
    is_tuple = ts isa Tuple
    if is_tuple
        ts = collect(ts)
        output_type = Tuple
    end

    allsyms = all_symbols(sys)
    if symbolic_type(ts) == NotSymbolic() && ts isa AbstractArray
        ts = map(x -> symbol_to_symbolic(sys, x; allsyms), ts)
    else
        ts = symbol_to_symbolic(sys, ts; allsyms)
    end

    vs = ModelingToolkit.vars(ts; op)
    namespace_subs = Dict()
    ns_map = Dict{Any, Any}(renamespace(sys, eq.lhs) => eq.lhs for eq in observed(sys))
    for sym in unknowns(sys)
        ns_map[renamespace(sys, sym)] = sym
        if iscall(sym) && operation(sym) === getindex
            ns_map[renamespace(sys, arguments(sym)[1])] = arguments(sym)[1]
        end
    end
    for sym in full_parameters(sys)
        ns_map[renamespace(sys, sym)] = sym
        if iscall(sym) && operation(sym) === getindex
            ns_map[renamespace(sys, arguments(sym)[1])] = arguments(sym)[1]
        end
    end
    allsyms = Set(all_symbols(sys))
    iv = has_iv(sys) ? get_iv(sys) : nothing
    for var in vs
        var = unwrap(var)
        newvar = get(ns_map, var, nothing)
        if newvar !== nothing
            namespace_subs[var] = newvar
            var = newvar
        end
        if throw && !var_in_varlist(var, allsyms, iv)
            Base.throw(ArgumentError("Symbol $var is not present in the system."))
        end
    end
    ts = fast_substitute(ts, namespace_subs)

    obsfilter = if param_only
        if is_split(sys)
            let ic = get_index_cache(sys)
                eq -> !(ContinuousTimeseries() in ic.observed_syms_to_timeseries[eq.lhs])
            end
        else
            Returns(false)
        end
    else
        Returns(true)
    end
    dvs = if param_only
        ()
    else
        (unknowns(sys),)
    end
    if inputs === nothing
        inputs = ()
    else
        ps = setdiff(ps, inputs) # Inputs have been converted to parameters by io_preprocessing, remove those from the parameter list
        inputs = (inputs,)
    end
    if disturbance_inputs !== nothing
        # Disturbance inputs may or may not be included as inputs, depending on disturbance_argument
        ps = setdiff(ps, disturbance_inputs)
    end
    if disturbance_argument
        disturbance_inputs = (disturbance_inputs,)
    else
        disturbance_inputs = ()
    end
    ps = reorder_parameters(sys, ps)
    iv = if is_time_dependent(sys)
        (get_iv(sys),)
    else
        ()
    end
    args = (dvs..., inputs..., ps..., iv..., disturbance_inputs...)
    p_start = length(dvs) + length(inputs) + 1
    p_end = length(dvs) + length(inputs) + length(ps)
    fns = build_function_wrapper(
        sys, ts, args...; p_start, p_end, filter_observed = obsfilter,
        output_type, mkarray, try_namespaced = true, expression = Val{true}, cse)
    if fns isa Tuple
        if expression
            return return_inplace ? fns : fns[1]
        end
        oop, iip = eval_or_rgf.(fns; eval_expression, eval_module)
        f = GeneratedFunctionWrapper{(
            p_start + is_dde(sys), length(args) - length(ps) + 1 + is_dde(sys), is_split(sys))}(
            oop, iip)
        return return_inplace ? (f, f) : f
    else
        if expression
            return fns
        end
        f = eval_or_rgf(fns; eval_expression, eval_module)
        f = GeneratedFunctionWrapper{(
            p_start + is_dde(sys), length(args) - length(ps) + 1 + is_dde(sys), is_split(sys))}(
            f, nothing)
        return f
    end
end
