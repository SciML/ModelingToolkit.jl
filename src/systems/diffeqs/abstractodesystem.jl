struct Schedule
    var_eq_matching::Any
    dummy_sub::Any
end

"""
    is_dde(sys::AbstractSystem)

Return a boolean indicating whether a system represents a set of delay
differential equations.
"""
is_dde(sys::AbstractSystem) = has_is_dde(sys) && get_is_dde(sys)

function _check_if_dde(eqs, iv, subsystems)
    is_dde = any(ModelingToolkit.is_dde, subsystems)
    if !is_dde
        vs = Set()
        for eq in eqs
            vars!(vs, eq)
            is_dde = any(vs) do sym
                isdelay(unwrap(sym), iv)
            end
            is_dde && break
        end
    end
    return is_dde
end

function filter_kwargs(kwargs)
    kwargs = Dict(kwargs)
    for key in keys(kwargs)
        key in DiffEqBase.allowedkeywords || delete!(kwargs, key)
    end
    pairs(NamedTuple(kwargs))
end
function gen_quoted_kwargs(kwargs)
    kwargparam = Expr(:parameters)
    for kw in kwargs
        push!(kwargparam.args, Expr(:kw, kw[1], kw[2]))
    end
    kwargparam
end

function calculate_tgrad(sys::AbstractODESystem;
        simplify = false)
    isempty(get_tgrad(sys)[]) || return get_tgrad(sys)[]  # use cached tgrad, if possible

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
    get_tgrad(sys)[] = tgrad
    return tgrad
end

function calculate_jacobian(sys::AbstractODESystem;
        sparse = false, simplify = false, dvs = unknowns(sys))
    if isequal(dvs, unknowns(sys))
        cache = get_jac(sys)[]
        if cache isa Tuple && cache[2] == (sparse, simplify)
            return cache[1]
        end
    end

    rhs = [eq.rhs - eq.lhs for eq in full_equations(sys)] #need du terms on rhs for differentiating wrt du

    if sparse
        jac = sparsejacobian(rhs, dvs, simplify = simplify)
        W_s = W_sparsity(sys)
        (Is, Js, Vs) = findnz(W_s)
        # Add nonzeros of W as non-structural zeros of the Jacobian (to ensure equal results for oop and iip Jacobian.)
        for (i, j) in zip(Is, Js)
            iszero(jac[i, j]) && begin
                jac[i, j] = 1
                jac[i, j] = 0
            end
        end
    else
        jac = jacobian(rhs, dvs, simplify = simplify)
    end

    if isequal(dvs, unknowns(sys))
        get_jac(sys)[] = jac, (sparse, simplify) # cache Jacobian
    end

    return jac
end

function calculate_control_jacobian(sys::AbstractODESystem;
        sparse = false, simplify = false)
    cache = get_ctrl_jac(sys)[]
    if cache isa Tuple && cache[2] == (sparse, simplify)
        return cache[1]
    end

    rhs = [eq.rhs for eq in full_equations(sys)]
    ctrls = controls(sys)

    if sparse
        jac = sparsejacobian(rhs, ctrls, simplify = simplify)
    else
        jac = jacobian(rhs, ctrls, simplify = simplify)
    end

    get_ctrl_jac(sys)[] = jac, (sparse, simplify) # cache Jacobian
    return jac
end

function generate_tgrad(
        sys::AbstractODESystem, dvs = unknowns(sys), ps = parameters(
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

function generate_jacobian(sys::AbstractODESystem, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true);
        simplify = false, sparse = false, kwargs...)
    jac = calculate_jacobian(sys; simplify = simplify, sparse = sparse)
    p = reorder_parameters(sys, ps)
    return build_function_wrapper(sys, jac,
        dvs,
        p...,
        get_iv(sys);
        wrap_code = sparse ? assert_jac_length_header(sys) : (identity, identity),
        kwargs...)
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

function generate_W(sys::AbstractODESystem, γ = 1., dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true); 
        simplify = false, sparse = false, kwargs...)
    @variables ˍ₋gamma
    M = calculate_massmatrix(sys; simplify)
    sparse && (M = SparseArrays.sparse(M))
    J = calculate_jacobian(sys; simplify, sparse, dvs)
    W = ˍ₋gamma*M + J

    p = reorder_parameters(sys, ps)
    return build_function_wrapper(sys, W, 
        dvs,
        p...,
        ˍ₋gamma,
        get_iv(sys);
        wrap_code = sparse ? assert_jac_length_header(sys) : (identity, identity),
        p_end = 1 + length(p),
        kwargs...)
end

function generate_control_jacobian(sys::AbstractODESystem, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true);
        simplify = false, sparse = false, kwargs...)
    jac = calculate_control_jacobian(sys; simplify = simplify, sparse = sparse)
    p = reorder_parameters(sys, ps)
    return build_function_wrapper(sys, jac, dvs, p..., get_iv(sys); kwargs...)
end

function generate_dae_jacobian(sys::AbstractODESystem, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true); simplify = false, sparse = false,
        kwargs...)
    jac_u = calculate_jacobian(sys; simplify = simplify, sparse = sparse)
    derivatives = Differential(get_iv(sys)).(unknowns(sys))
    jac_du = calculate_jacobian(sys; simplify = simplify, sparse = sparse,
        dvs = derivatives)
    dvs = unknowns(sys)
    @variables ˍ₋gamma
    jac = ˍ₋gamma * jac_du + jac_u
    pre = get_preprocess_constants(jac)
    p = reorder_parameters(sys, ps)
    return build_function_wrapper(sys, jac, derivatives, dvs, p..., ˍ₋gamma, get_iv(sys);
        p_start = 3, p_end = 2 + length(p), kwargs...)
end

function generate_function(sys::AbstractODESystem, dvs = unknowns(sys),
        ps = parameters(sys; initial_parameters = true);
        implicit_dae = false,
        ddvs = implicit_dae ? map(Differential(get_iv(sys)), dvs) :
               nothing,
        isdde = false,
        kwargs...)
    eqs = [eq for eq in equations(sys)]
    if !implicit_dae
        check_operator_variables(eqs, Differential)
        check_lhs(eqs, Differential, Set(dvs))
    end

    rhss = implicit_dae ? [_iszero(eq.lhs) ? eq.rhs : eq.rhs - eq.lhs for eq in eqs] :
           [eq.rhs for eq in eqs]

    if !isempty(assertions(sys))
        rhss[end] += unwrap(get_assertions_expr(sys))
    end

    # TODO: add an optional check on the ordering of observed equations
    u = dvs
    p = reorder_parameters(sys, ps)
    t = get_iv(sys)

    if implicit_dae
        build_function_wrapper(sys, rhss, ddvs, u, p..., t; p_start = 3, kwargs...)
    else
        build_function_wrapper(sys, rhss, u, p..., t; kwargs...)
    end
end

function isdelay(var, iv)
    iv === nothing && return false
    isvariable(var) || return false
    isparameter(var) && return false
    if iscall(var) && !ModelingToolkit.isoperator(var, Symbolics.Operator)
        args = arguments(var)
        length(args) == 1 || return false
        isequal(args[1], iv) || return true
    end
    return false
end
const DDE_HISTORY_FUN = Sym{Symbolics.FnType{Tuple{Any, <:Real}, Vector{Real}}}(:___history___)
const DEFAULT_PARAMS_ARG = Sym{Any}(:ˍ₋arg3)
function delay_to_function(
        sys::AbstractODESystem, eqs = full_equations(sys); history_arg = DEFAULT_PARAMS_ARG)
    delay_to_function(eqs,
        get_iv(sys),
        Dict{Any, Int}(operation(s) => i for (i, s) in enumerate(unknowns(sys))),
        parameters(sys),
        DDE_HISTORY_FUN; history_arg)
end
function delay_to_function(eqs::Vector, iv, sts, ps, h; history_arg = DEFAULT_PARAMS_ARG)
    delay_to_function.(eqs, (iv,), (sts,), (ps,), (h,); history_arg)
end
function delay_to_function(eq::Equation, iv, sts, ps, h; history_arg = DEFAULT_PARAMS_ARG)
    delay_to_function(eq.lhs, iv, sts, ps, h; history_arg) ~ delay_to_function(
        eq.rhs, iv, sts, ps, h; history_arg)
end
function delay_to_function(expr, iv, sts, ps, h; history_arg = DEFAULT_PARAMS_ARG)
    if isdelay(expr, iv)
        v = operation(expr)
        time = arguments(expr)[1]
        idx = sts[v]
        return term(getindex, h(history_arg, time), idx, type = Real) # BIG BIG HACK
    elseif iscall(expr)
        return maketerm(typeof(expr),
            operation(expr),
            map(x -> delay_to_function(x, iv, sts, ps, h; history_arg), arguments(expr)),
            metadata(expr))
    else
        return expr
    end
end

function calculate_massmatrix(sys::AbstractODESystem; simplify = false)
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
    M = simplify ? ModelingToolkit.simplify.(M) : M
    # M should only contain concrete numbers
    M == I ? I : M
end

function jacobian_sparsity(sys::AbstractODESystem)
    sparsity = torn_system_jacobian_sparsity(sys)
    sparsity === nothing || return sparsity

    jacobian_sparsity([eq.rhs for eq in full_equations(sys)],
        [dv for dv in unknowns(sys)])
end

function jacobian_dae_sparsity(sys::AbstractODESystem)
    J1 = jacobian_sparsity([eq.rhs for eq in full_equations(sys)],
        [dv for dv in unknowns(sys)])
    derivatives = Differential(get_iv(sys)).(unknowns(sys))
    J2 = jacobian_sparsity([eq.rhs for eq in full_equations(sys)],
        [dv for dv in derivatives])
    J1 + J2
end

function W_sparsity(sys::AbstractODESystem) 
    jac_sparsity = jacobian_sparsity(sys)
    (n, n) = size(jac_sparsity)
    M = calculate_massmatrix(sys)
    M_sparsity = M isa UniformScaling ? sparse(I(n)) : SparseMatrixCSC{Bool, Int64}((!iszero).(M))
    jac_sparsity .| M_sparsity
end

function isautonomous(sys::AbstractODESystem)
    tgrad = calculate_tgrad(sys; simplify = true)
    all(iszero, tgrad)
end

"""
```julia
DiffEqBase.ODEFunction{iip}(sys::AbstractODESystem, dvs = unknowns(sys),
                            ps = parameters(sys);
                            version = nothing, tgrad = false,
                            jac = false,
                            sparse = false,
                            kwargs...) where {iip}
```

Create an `ODEFunction` from the [`ODESystem`](@ref). The arguments `dvs` and `ps`
are used to set the order of the dependent variable and parameter vectors,
respectively.
"""
function DiffEqBase.ODEFunction(sys::AbstractODESystem, args...; kwargs...)
    ODEFunction{true}(sys, args...; kwargs...)
end

function DiffEqBase.ODEFunction{true}(sys::AbstractODESystem, args...;
        kwargs...)
    ODEFunction{true, SciMLBase.AutoSpecialize}(sys, args...; kwargs...)
end

function DiffEqBase.ODEFunction{false}(sys::AbstractODESystem, args...;
        kwargs...)
    ODEFunction{false, SciMLBase.FullSpecialize}(sys, args...; kwargs...)
end

function DiffEqBase.ODEFunction{iip, specialize}(sys::AbstractODESystem,
        dvs = unknowns(sys),
        ps = parameters(sys), u0 = nothing;
        version = nothing, tgrad = false,
        jac = false, p = nothing,
        t = nothing,
        eval_expression = false,
        sparse = false, simplify = false,
        eval_module = @__MODULE__,
        steady_state = false,
        checkbounds = false,
        sparsity = false,
        analytic = nothing,
        split_idxs = nothing,
        initialization_data = nothing,
        cse = true,
        kwargs...) where {iip, specialize}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system before creating an `ODEFunction`")
    end
    f_gen = generate_function(sys, dvs, ps; expression = Val{true},
        expression_module = eval_module, checkbounds = checkbounds, cse,
        kwargs...)
    f_oop, f_iip = eval_or_rgf.(f_gen; eval_expression, eval_module)
    f = GeneratedFunctionWrapper{(2, 3, is_split(sys))}(f_oop, f_iip)

    if specialize === SciMLBase.FunctionWrapperSpecialize && iip
        if u0 === nothing || p === nothing || t === nothing
            error("u0, p, and t must be specified for FunctionWrapperSpecialize on ODEFunction.")
        end
        f = SciMLBase.wrapfun_iip(f, (u0, u0, p, t))
    end

    if tgrad
        tgrad_gen = generate_tgrad(sys, dvs, ps;
            simplify = simplify,
            expression = Val{true},
            expression_module = eval_module, cse,
            checkbounds = checkbounds, kwargs...)
        tgrad_oop, tgrad_iip = eval_or_rgf.(tgrad_gen; eval_expression, eval_module)
        _tgrad = GeneratedFunctionWrapper{(2, 3, is_split(sys))}(tgrad_oop, tgrad_iip)
    else
        _tgrad = nothing
    end

    if jac
        jac_gen = generate_jacobian(sys, dvs, ps;
            simplify = simplify, sparse = sparse,
            expression = Val{true},
            expression_module = eval_module, cse,
            checkbounds = checkbounds, kwargs...)
        jac_oop, jac_iip = eval_or_rgf.(jac_gen; eval_expression, eval_module)

        _jac = GeneratedFunctionWrapper{(2, 3, is_split(sys))}(jac_oop, jac_iip)
    else
        _jac = nothing
    end

    M = calculate_massmatrix(sys)

    _M = if sparse && !(u0 === nothing || M === I)
        SparseArrays.sparse(M)
    elseif u0 === nothing || M === I
        M
    else
        ArrayInterface.restructure(u0 .* u0', M)
    end

    observedfun = ObservedFunctionCache(
        sys; steady_state, eval_expression, eval_module, checkbounds, cse)

    if sparse
        uElType = u0 === nothing ? Float64 : eltype(u0)
        W_prototype = similar(W_sparsity(sys), uElType)
    else
        W_prototype = nothing
    end

    @set! sys.split_idxs = split_idxs

    ODEFunction{iip, specialize}(f;
        sys = sys,
        jac = _jac === nothing ? nothing : _jac,
        tgrad = _tgrad === nothing ? nothing : _tgrad,
        mass_matrix = _M,
        jac_prototype = W_prototype,
        observed = observedfun,
        sparsity = sparsity ? W_sparsity(sys) : nothing,
        analytic = analytic,
        initialization_data)
end

"""
```julia
DiffEqBase.DAEFunction{iip}(sys::AbstractODESystem, dvs = unknowns(sys),
                            ps = parameters(sys);
                            version = nothing, tgrad = false,
                            jac = false,
                            sparse = false,
                            kwargs...) where {iip}
```

Create an `DAEFunction` from the [`ODESystem`](@ref). The arguments `dvs` and
`ps` are used to set the order of the dependent variable and parameter vectors,
respectively.
"""
function DiffEqBase.DAEFunction(sys::AbstractODESystem, args...; kwargs...)
    DAEFunction{true}(sys, args...; kwargs...)
end

function DiffEqBase.DAEFunction{iip}(sys::AbstractODESystem, dvs = unknowns(sys),
        ps = parameters(sys), u0 = nothing;
        ddvs = map(Base.Fix2(diff2term, get_iv(sys)) ∘ Differential(get_iv(sys)), dvs),
        version = nothing, p = nothing,
        jac = false,
        eval_expression = false,
        sparse = false, simplify = false,
        eval_module = @__MODULE__,
        checkbounds = false,
        initialization_data = nothing,
        cse = true,
        kwargs...) where {iip}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system before creating a `DAEFunction`")
    end
    f_gen = generate_function(sys, dvs, ps; implicit_dae = true,
        expression = Val{true}, cse,
        expression_module = eval_module, checkbounds = checkbounds,
        kwargs...)
    f_oop, f_iip = eval_or_rgf.(f_gen; eval_expression, eval_module)
    f = GeneratedFunctionWrapper{(3, 4, is_split(sys))}(f_oop, f_iip)

    if jac
        jac_gen = generate_dae_jacobian(sys, dvs, ps;
            simplify = simplify, sparse = sparse,
            expression = Val{true},
            expression_module = eval_module, cse,
            checkbounds = checkbounds, kwargs...)
        jac_oop, jac_iip = eval_or_rgf.(jac_gen; eval_expression, eval_module)

        _jac = GeneratedFunctionWrapper{(3, 5, is_split(sys))}(jac_oop, jac_iip)
    else
        _jac = nothing
    end

    observedfun = ObservedFunctionCache(
        sys; eval_expression, eval_module, checkbounds = get(kwargs, :checkbounds, false), cse)

    jac_prototype = if sparse
        uElType = u0 === nothing ? Float64 : eltype(u0)
        if jac
            J1 = calculate_jacobian(sys, sparse = sparse)
            derivatives = Differential(get_iv(sys)).(unknowns(sys))
            J2 = calculate_jacobian(sys; sparse = sparse, dvs = derivatives)
            similar(J1 + J2, uElType)
        else
            similar(jacobian_dae_sparsity(sys), uElType)
        end
    else
        nothing
    end

    DAEFunction{iip}(f;
        sys = sys,
        jac = _jac === nothing ? nothing : _jac,
        jac_prototype = jac_prototype,
        observed = observedfun,
        initialization_data)
end

function DiffEqBase.DDEFunction(sys::AbstractODESystem, args...; kwargs...)
    DDEFunction{true}(sys, args...; kwargs...)
end

function DiffEqBase.DDEFunction{iip}(sys::AbstractODESystem, dvs = unknowns(sys),
        ps = parameters(sys), u0 = nothing;
        eval_expression = false,
        eval_module = @__MODULE__,
        checkbounds = false,
        initialization_data = nothing,
        cse = true,
        kwargs...) where {iip}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system before creating an `DDEFunction`")
    end
    f_gen = generate_function(sys, dvs, ps; isdde = true,
        expression = Val{true},
        expression_module = eval_module, checkbounds = checkbounds,
        cse, kwargs...)
    f_oop, f_iip = eval_or_rgf.(f_gen; eval_expression, eval_module)
    f = GeneratedFunctionWrapper{(3, 4, is_split(sys))}(f_oop, f_iip)

    DDEFunction{iip}(f; sys = sys, initialization_data)
end

function DiffEqBase.SDDEFunction(sys::AbstractODESystem, args...; kwargs...)
    SDDEFunction{true}(sys, args...; kwargs...)
end

function DiffEqBase.SDDEFunction{iip}(sys::AbstractODESystem, dvs = unknowns(sys),
        ps = parameters(sys), u0 = nothing;
        eval_expression = false,
        eval_module = @__MODULE__,
        checkbounds = false,
        initialization_data = nothing,
        cse = true,
        kwargs...) where {iip}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system before creating an `SDDEFunction`")
    end
    f_gen = generate_function(sys, dvs, ps; isdde = true,
        expression = Val{true},
        expression_module = eval_module, checkbounds = checkbounds,
        cse, kwargs...)
    f_oop, f_iip = eval_or_rgf.(f_gen; eval_expression, eval_module)
    f = GeneratedFunctionWrapper{(3, 4, is_split(sys))}(f_oop, f_iip)

    g_gen = generate_diffusion_function(sys, dvs, ps; expression = Val{true},
        isdde = true, cse, kwargs...)
    g_oop, g_iip = eval_or_rgf.(g_gen; eval_expression, eval_module)
    g = GeneratedFunctionWrapper{(3, 4, is_split(sys))}(g_oop, g_iip)

    SDDEFunction{iip}(f, g; sys = sys, initialization_data)
end

"""
```julia
ODEFunctionExpr{iip}(sys::AbstractODESystem, dvs = unknowns(sys),
                     ps = parameters(sys);
                     version = nothing, tgrad = false,
                     jac = false,
                     sparse = false,
                     kwargs...) where {iip}
```

Create a Julia expression for an `ODEFunction` from the [`ODESystem`](@ref).
The arguments `dvs` and `ps` are used to set the order of the dependent
variable and parameter vectors, respectively.
"""
struct ODEFunctionExpr{iip, specialize} end

function ODEFunctionExpr{iip, specialize}(sys::AbstractODESystem, dvs = unknowns(sys),
        ps = parameters(sys), u0 = nothing;
        version = nothing, tgrad = false,
        jac = false, p = nothing,
        linenumbers = false,
        sparse = false, simplify = false,
        steady_state = false,
        sparsity = false,
        observedfun_exp = nothing,
        kwargs...) where {iip, specialize}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system before creating an `ODEFunctionExpr`")
    end
    f_oop, f_iip = generate_function(sys, dvs, ps; expression = Val{true}, kwargs...)

    fsym = gensym(:f)
    _f = :($fsym = $(GeneratedFunctionWrapper{(2, 3, is_split(sys))})($f_oop, $f_iip))
    tgradsym = gensym(:tgrad)
    if tgrad
        tgrad_oop, tgrad_iip = generate_tgrad(sys, dvs, ps;
            simplify = simplify,
            expression = Val{true}, kwargs...)
        _tgrad = :($tgradsym = $(GeneratedFunctionWrapper{(2, 3, is_split(sys))})(
            $tgrad_oop, $tgrad_iip))
    else
        _tgrad = :($tgradsym = nothing)
    end

    jacsym = gensym(:jac)
    if jac
        jac_oop, jac_iip = generate_jacobian(sys, dvs, ps;
            sparse = sparse, simplify = simplify,
            expression = Val{true}, kwargs...)
        _jac = :($jacsym = $(GeneratedFunctionWrapper{(2, 3, is_split(sys))})(
            $jac_oop, $jac_iip))
    else
        _jac = :($jacsym = nothing)
    end

    Msym = gensym(:M)
    M = calculate_massmatrix(sys)
    if sparse && !(u0 === nothing || M === I)
        _M = :($Msym = $(SparseArrays.sparse(M)))
    elseif u0 === nothing || M === I
        _M = :($Msym = $M)
    else
        _M = :($Msym = $(ArrayInterface.restructure(u0 .* u0', M)))
    end

    jp_expr = sparse ? :($similar($(get_jac(sys)[]), Float64)) : :nothing
    ex = quote
        let $_f, $_tgrad, $_jac, $_M
            ODEFunction{$iip, $specialize}($fsym,
                sys = $sys,
                jac = $jacsym,
                tgrad = $tgradsym,
                mass_matrix = $Msym,
                jac_prototype = $jp_expr,
                sparsity = $(sparsity ? jacobian_sparsity(sys) : nothing),
                observed = $observedfun_exp)
        end
    end
    !linenumbers ? Base.remove_linenums!(ex) : ex
end

function ODEFunctionExpr(sys::AbstractODESystem, args...; kwargs...)
    ODEFunctionExpr{true}(sys, args...; kwargs...)
end

function ODEFunctionExpr{true}(sys::AbstractODESystem, args...; kwargs...)
    return ODEFunctionExpr{true, SciMLBase.AutoSpecialize}(sys, args...; kwargs...)
end

function ODEFunctionExpr{false}(sys::AbstractODESystem, args...; kwargs...)
    return ODEFunctionExpr{false, SciMLBase.FullSpecialize}(sys, args...; kwargs...)
end

"""
```julia
DAEFunctionExpr{iip}(sys::AbstractODESystem, dvs = unknowns(sys),
                     ps = parameters(sys);
                     version = nothing, tgrad = false,
                     jac = false,
                     sparse = false,
                     kwargs...) where {iip}
```

Create a Julia expression for an `ODEFunction` from the [`ODESystem`](@ref).
The arguments `dvs` and `ps` are used to set the order of the dependent
variable and parameter vectors, respectively.
"""
struct DAEFunctionExpr{iip} end

function DAEFunctionExpr{iip}(sys::AbstractODESystem, dvs = unknowns(sys),
        ps = parameters(sys), u0 = nothing;
        version = nothing, tgrad = false,
        jac = false, p = nothing,
        linenumbers = false,
        sparse = false, simplify = false,
        kwargs...) where {iip}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system before creating an `DAEFunctionExpr`")
    end
    f_oop, f_iip = generate_function(sys, dvs, ps; expression = Val{true},
        implicit_dae = true, kwargs...)
    fsym = gensym(:f)
    _f = :($fsym = $(GeneratedFunctionWrapper{(3, 4, is_split(sys))})($f_oop, $f_iip))
    ex = quote
        $_f
        ODEFunction{$iip}($fsym)
    end
    !linenumbers ? Base.remove_linenums!(ex) : ex
end

function DAEFunctionExpr(sys::AbstractODESystem, args...; kwargs...)
    DAEFunctionExpr{true}(sys, args...; kwargs...)
end

struct SymbolicTstops{F}
    fn::F
end

function (st::SymbolicTstops)(p, tspan)
    unique!(sort!(reduce(vcat, st.fn(p, tspan...))))
end

function SymbolicTstops(
        sys::AbstractSystem; eval_expression = false, eval_module = @__MODULE__)
    tstops = symbolic_tstops(sys)
    isempty(tstops) && return nothing
    t0 = gensym(:t0)
    t1 = gensym(:t1)
    tstops = map(tstops) do val
        if is_array_of_symbolics(val) || val isa AbstractArray
            collect(val)
        else
            term(:, t0, unwrap(val), t1; type = AbstractArray{Real})
        end
    end
    rps = reorder_parameters(sys)
    tstops, _ = build_function_wrapper(sys, tstops,
        rps...,
        t0,
        t1;
        expression = Val{true},
        p_start = 1, p_end = length(rps), add_observed = false, force_SA = true)
    tstops = eval_or_rgf(tstops; eval_expression, eval_module)
    tstops = GeneratedFunctionWrapper{(1, 3, is_split(sys))}(tstops, nothing)
    return SymbolicTstops(tstops)
end

"""
```julia
DiffEqBase.ODEProblem{iip}(sys::AbstractODESystem, u0map, tspan,
                           parammap = DiffEqBase.NullParameters();
                           allow_cost = false,
                           version = nothing, tgrad = false,
                           jac = false,
                           checkbounds = false, sparse = false,
                           simplify = false,
                           linenumbers = true, parallel = SerialForm(),
                           kwargs...) where {iip}
```

Generates an ODEProblem from an ODESystem and allows for automatically
symbolically calculating numerical enhancements.
"""
function DiffEqBase.ODEProblem(sys::AbstractODESystem, args...; kwargs...)
    ODEProblem{true}(sys, args...; kwargs...)
end

function DiffEqBase.ODEProblem(sys::AbstractODESystem,
        u0map::StaticArray,
        args...;
        kwargs...)
    ODEProblem{false, SciMLBase.FullSpecialize}(sys, u0map, args...; kwargs...)
end

function DiffEqBase.ODEProblem{true}(sys::AbstractODESystem, args...; kwargs...)
    ODEProblem{true, SciMLBase.AutoSpecialize}(sys, args...; kwargs...)
end

function DiffEqBase.ODEProblem{false}(sys::AbstractODESystem, args...; kwargs...)
    ODEProblem{false, SciMLBase.FullSpecialize}(sys, args...; kwargs...)
end

function DiffEqBase.ODEProblem{iip, specialize}(sys::AbstractODESystem, u0map = [],
        tspan = get_tspan(sys),
        parammap = DiffEqBase.NullParameters();
        allow_cost = false,
        callback = nothing,
        check_length = true,
        warn_initialize_determined = true,
        eval_expression = false,
        eval_module = @__MODULE__,
        kwargs...) where {iip, specialize}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system before creating an `ODEProblem`")
    end

    if !isnothing(get_constraintsystem(sys))
        error("An ODESystem with constraints cannot be used to construct a regular ODEProblem. 
              Consider a BVProblem instead.")
    end

    if !isempty(get_costs(sys)) && !allow_cost
        error("ODEProblem will not optimize solutions of ODESystems that have associated cost functions. 
              Solvers for optimal control problems are forthcoming. In order to bypass this error (e.g. 
              to check the cost of a regular solution), pass `allow_cost` = true into the constructor.")
    end

    f, u0, p = process_SciMLProblem(ODEFunction{iip, specialize}, sys, u0map, parammap;
        t = tspan !== nothing ? tspan[1] : tspan,
        check_length, warn_initialize_determined, eval_expression, eval_module, kwargs...)
    cbs = process_events(sys; callback, eval_expression, eval_module, kwargs...)

    kwargs = filter_kwargs(kwargs)
    pt = something(get_metadata(sys), StandardODEProblem())

    kwargs1 = (;)
    if cbs !== nothing
        kwargs1 = merge(kwargs1, (callback = cbs,))
    end

    tstops = SymbolicTstops(sys; eval_expression, eval_module)
    if tstops !== nothing
        kwargs1 = merge(kwargs1, (; tstops))
    end

    # Call `remake` so it runs initialization if it is trivial
    return remake(ODEProblem{iip}(f, u0, tspan, p, pt; kwargs1..., kwargs...))
end
get_callback(prob::ODEProblem) = prob.kwargs[:callback]

"""
```julia
SciMLBase.BVProblem{iip}(sys::AbstractODESystem, u0map, tspan,
                         parammap = DiffEqBase.NullParameters();
                         constraints = nothing, guesses = nothing,
                         version = nothing, tgrad = false,
                         jac = true, sparse = true,
                         simplify = false,
                         kwargs...) where {iip}
```

Create a boundary value problem from the [`ODESystem`](@ref). 

`u0map` is used to specify fixed initial values for the states. Every variable 
must have either an initial guess supplied using `guesses` or a fixed initial 
value specified using `u0map`.

Boundary value conditions are supplied to ODESystems
in the form of a ConstraintsSystem. These equations 
should specify values that state variables should
take at specific points, as in `x(0.5) ~ 1`). More general constraints that 
should hold over the entire solution, such as `x(t)^2 + y(t)^2`, should be 
specified as one of the equations used to build the `ODESystem`.

If an ODESystem without `constraints` is specified, it will be treated as an initial value problem. 

```julia
    @parameters g t_c = 0.5
    @variables x(..) y(t) λ(t)
    eqs = [D(D(x(t))) ~ λ * x(t)
           D(D(y)) ~ λ * y - g
           x(t)^2 + y^2 ~ 1]
    cstr = [x(0.5) ~ 1]
    @mtkbuild pend = ODESystem(eqs, t; constraints = cstrs)

    tspan = (0.0, 1.5)
    u0map = [x(t) => 0.6, y => 0.8]
    parammap = [g => 1]
    guesses = [λ => 1]

    bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; guesses, check_length = false)
```

If the `ODESystem` has algebraic equations, like `x(t)^2 + y(t)^2`, the resulting 
`BVProblem` must be solved using BVDAE solvers, such as Ascher.
"""
function SciMLBase.BVProblem(sys::AbstractODESystem, args...; kwargs...)
    BVProblem{true}(sys, args...; kwargs...)
end

function SciMLBase.BVProblem(sys::AbstractODESystem,
        u0map::StaticArray,
        args...;
        kwargs...)
    BVProblem{false, SciMLBase.FullSpecialize}(sys, u0map, args...; kwargs...)
end

function SciMLBase.BVProblem{true}(sys::AbstractODESystem, args...; kwargs...)
    BVProblem{true, SciMLBase.AutoSpecialize}(sys, args...; kwargs...)
end

function SciMLBase.BVProblem{false}(sys::AbstractODESystem, args...; kwargs...)
    BVProblem{false, SciMLBase.FullSpecialize}(sys, args...; kwargs...)
end

function SciMLBase.BVProblem{iip, specialize}(sys::AbstractODESystem, u0map = [],
        tspan = get_tspan(sys),
        parammap = DiffEqBase.NullParameters();
        guesses = Dict(),
        allow_cost = false,
        version = nothing, tgrad = false,
        callback = nothing,
        check_length = true,
        warn_initialize_determined = true,
        eval_expression = false,
        eval_module = @__MODULE__,
        cse = true,
        kwargs...) where {iip, specialize}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system before creating an `BVProblem`")
    end
    !isnothing(callback) && error("BVP solvers do not support callbacks.")

    if !isempty(get_costs(sys)) && !allow_cost
        error("BVProblem will not optimize solutions of ODESystems that have associated cost functions.
              Solvers for optimal control problems are forthcoming. In order to bypass this error (e.g.
              to check the cost of a regular solution), pass `allow_cost` = true into the constructor.")
    end

    has_alg_eqs(sys) &&
        error("The BVProblem constructor currently does not support ODESystems with algebraic equations.") # Remove this when the BVDAE solvers get updated, the codegen should work when it does.

    sts = unknowns(sys)
    ps = parameters(sys)
    constraintsys = get_constraintsystem(sys)

    if !isnothing(constraintsys)
        (length(constraints(constraintsys)) + length(u0map) > length(sts)) &&
            @warn "The BVProblem is overdetermined. The total number of conditions (# constraints + # fixed initial values given by u0map) exceeds the total number of states. The BVP solvers will default to doing a nonlinear least-squares optimization."
    end

    # ODESystems without algebraic equations should use both fixed values + guesses
    # for initialization.
    _u0map = has_alg_eqs(sys) ? u0map : merge(Dict(u0map), Dict(guesses))
    f, u0, p = process_SciMLProblem(ODEFunction{iip, specialize}, sys, _u0map, parammap;
        t = tspan !== nothing ? tspan[1] : tspan, guesses,
        check_length, warn_initialize_determined, eval_expression, eval_module, cse, kwargs...)

    stidxmap = Dict([v => i for (i, v) in enumerate(sts)])
    u0_idxs = has_alg_eqs(sys) ? collect(1:length(sts)) : [stidxmap[k] for (k, v) in u0map]

    fns = generate_function_bc(sys, u0, u0_idxs, tspan; cse)
    bc_oop, bc_iip = eval_or_rgf.(fns; eval_expression, eval_module)
    bc(sol, p, t) = bc_oop(sol, p, t)
    bc(resid, u, p, t) = bc_iip(resid, u, p, t)

    return BVProblem{iip}(f, bc, u0, tspan, p; kwargs...)
end

get_callback(prob::BVProblem) = error("BVP solvers do not support callbacks.")

"""
    generate_function_bc(sys::ODESystem, u0, u0_idxs, tspan)

    Given an ODESystem with constraints, generate the boundary condition function to pass to boundary value problem solvers.
    Expression uses the constraints and the provided initial conditions.
"""
function generate_function_bc(sys::ODESystem, u0, u0_idxs, tspan; kwargs...)
    iv = get_iv(sys)
    sts = unknowns(sys)
    ps = parameters(sys)
    np = length(ps)
    ns = length(sts)
    stidxmap = Dict([v => i for (i, v) in enumerate(sts)])
    pidxmap = Dict([v => i for (i, v) in enumerate(ps)])

    @variables sol(..)[1:ns]

    conssys = get_constraintsystem(sys)
    cons = Any[]
    if !isnothing(conssys)
        cons = [con.lhs - con.rhs for con in constraints(conssys)]

        for st in get_unknowns(conssys)
            x = operation(st)
            t = only(arguments(st))
            idx = stidxmap[x(iv)]

            cons = map(c -> Symbolics.substitute(c, Dict(x(t) => sol(t)[idx])), cons)
        end
    end

    init_conds = Any[]
    for i in u0_idxs
        expr = sol(tspan[1])[i] - u0[i]
        push!(init_conds, expr)
    end

    exprs = vcat(init_conds, cons)
    _p = reorder_parameters(sys, ps)

    build_function_wrapper(sys, exprs, sol, _p..., iv; output_type = Array, kwargs...)
end

"""
```julia
DiffEqBase.DAEProblem{iip}(sys::AbstractODESystem, du0map, u0map, tspan,
                           parammap = DiffEqBase.NullParameters();
                           version = nothing, tgrad = false,
                           jac = false,
                           checkbounds = false, sparse = false,
                           simplify = false,
                           linenumbers = true, parallel = SerialForm(),
                           kwargs...) where {iip}
```

Generates a DAEProblem from an ODESystem and allows for automatically
symbolically calculating numerical enhancements.

Note: Solvers for DAEProblems like DFBDF, DImplicitEuler, DABDF2 are 
generally slower than the ones for ODEProblems. We recommend trying 
ODEProblem and its solvers for your problem first.
"""
function DiffEqBase.DAEProblem(sys::AbstractODESystem, args...; kwargs...)
    DAEProblem{true}(sys, args...; kwargs...)
end

function DiffEqBase.DAEProblem{iip}(sys::AbstractODESystem, du0map, u0map, tspan,
        parammap = DiffEqBase.NullParameters();
        allow_cost = false,
        warn_initialize_determined = true,
        check_length = true, eval_expression = false, eval_module = @__MODULE__, kwargs...) where {iip}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system before creating a `DAEProblem`.")
    end

    if !isempty(get_costs(sys)) && !allow_cost
        error("DAEProblem will not optimize solutions of ODESystems that have associated cost functions. 
              Solvers for optimal control problems are forthcoming. In order to bypass this error (e.g. 
              to check the cost of a regular solution), pass `allow_cost` = true into the constructor.")
    end

    f, du0, u0, p = process_SciMLProblem(DAEFunction{iip}, sys, u0map, parammap;
        implicit_dae = true, du0map = du0map, check_length,
        t = tspan !== nothing ? tspan[1] : tspan,
        warn_initialize_determined, kwargs...)
    diffvars = collect_differential_variables(sys)
    sts = unknowns(sys)
    differential_vars = map(Base.Fix2(in, diffvars), sts)
    kwargs = filter_kwargs(kwargs)

    kwargs1 = (;)

    tstops = SymbolicTstops(sys; eval_expression, eval_module)
    if tstops !== nothing
        kwargs1 = merge(kwargs1, (; tstops))
    end

    # Call `remake` so it runs initialization if it is trivial
    return remake(DAEProblem{iip}(
        f, du0, u0, tspan, p; differential_vars = differential_vars,
        kwargs..., kwargs1...))
end

function generate_history(sys::AbstractODESystem, u0; expression = Val{false}, kwargs...)
    p = reorder_parameters(sys)
    build_function_wrapper(
        sys, u0, p..., get_iv(sys); expression, p_start = 1, p_end = length(p),
        similarto = typeof(u0), wrap_delays = false, kwargs...)
end

function DiffEqBase.DDEProblem(sys::AbstractODESystem, args...; kwargs...)
    DDEProblem{true}(sys, args...; kwargs...)
end
function DiffEqBase.DDEProblem{iip}(sys::AbstractODESystem, u0map = [],
        tspan = get_tspan(sys),
        parammap = DiffEqBase.NullParameters();
        callback = nothing,
        check_length = true,
        eval_expression = false,
        eval_module = @__MODULE__,
        u0_constructor = identity,
        cse = true,
        kwargs...) where {iip}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system before creating a `DDEProblem`")
    end
    f, u0, p = process_SciMLProblem(DDEFunction{iip}, sys, u0map, parammap;
        t = tspan !== nothing ? tspan[1] : tspan,
        symbolic_u0 = true, u0_constructor, cse,
        check_length, eval_expression, eval_module, kwargs...)
    h_gen = generate_history(sys, u0; expression = Val{true}, cse)
    h_oop, h_iip = eval_or_rgf.(h_gen; eval_expression, eval_module)
    h = h_oop
    u0 = float.(h(p, tspan[1]))
    if u0 !== nothing
        u0 = u0_constructor(u0)
    end

    cbs = process_events(sys; callback, eval_expression, eval_module, kwargs...)
    kwargs = filter_kwargs(kwargs)

    kwargs1 = (;)
    if cbs !== nothing
        kwargs1 = merge(kwargs1, (callback = cbs,))
    end
    # Call `remake` so it runs initialization if it is trivial
    return remake(DDEProblem{iip}(f, u0, h, tspan, p; kwargs1..., kwargs...))
end

function DiffEqBase.SDDEProblem(sys::AbstractODESystem, args...; kwargs...)
    SDDEProblem{true}(sys, args...; kwargs...)
end
function DiffEqBase.SDDEProblem{iip}(sys::AbstractODESystem, u0map = [],
        tspan = get_tspan(sys),
        parammap = DiffEqBase.NullParameters();
        callback = nothing,
        check_length = true,
        sparsenoise = nothing,
        eval_expression = false,
        eval_module = @__MODULE__,
        u0_constructor = identity,
        cse = true,
        kwargs...) where {iip}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system before creating a `SDDEProblem`")
    end
    f, u0, p = process_SciMLProblem(SDDEFunction{iip}, sys, u0map, parammap;
        t = tspan !== nothing ? tspan[1] : tspan,
        symbolic_u0 = true, eval_expression, eval_module, u0_constructor,
        check_length, cse, kwargs...)
    h_gen = generate_history(sys, u0; expression = Val{true}, cse)
    h_oop, h_iip = eval_or_rgf.(h_gen; eval_expression, eval_module)
    h = h_oop
    u0 = h(p, tspan[1])
    if u0 !== nothing
        u0 = u0_constructor(u0)
    end

    cbs = process_events(sys; callback, eval_expression, eval_module, kwargs...)
    kwargs = filter_kwargs(kwargs)

    kwargs1 = (;)
    if cbs !== nothing
        kwargs1 = merge(kwargs1, (callback = cbs,))
    end

    noiseeqs = get_noiseeqs(sys)
    sparsenoise === nothing && (sparsenoise = get(kwargs, :sparse, false))
    if noiseeqs isa AbstractVector
        noise_rate_prototype = nothing
    elseif sparsenoise
        I, J, V = findnz(SparseArrays.sparse(noiseeqs))
        noise_rate_prototype = SparseArrays.sparse(I, J, zero(eltype(u0)))
    else
        noise_rate_prototype = zeros(eltype(u0), size(noiseeqs))
    end
    # Call `remake` so it runs initialization if it is trivial
    return remake(SDDEProblem{iip}(f, f.g, u0, h, tspan, p;
        noise_rate_prototype =
        noise_rate_prototype, kwargs1..., kwargs...))
end

"""
```julia
ODEProblemExpr{iip}(sys::AbstractODESystem, u0map, tspan,
                    parammap = DiffEqBase.NullParameters();
                    version = nothing, tgrad = false,
                    jac = false,
                    checkbounds = false, sparse = false,
                    linenumbers = true, parallel = SerialForm(),
                    skipzeros = true, fillzeros = true,
                    simplify = false,
                    kwargs...) where {iip}
```

Generates a Julia expression for constructing an ODEProblem from an
ODESystem and allows for automatically symbolically calculating
numerical enhancements.
"""
struct ODEProblemExpr{iip} end

function ODEProblemExpr{iip}(sys::AbstractODESystem, u0map, tspan,
        parammap = DiffEqBase.NullParameters(); check_length = true,
        kwargs...) where {iip}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system before creating a `ODEProblemExpr`")
    end
    f, u0, p = process_SciMLProblem(
        ODEFunctionExpr{iip}, sys, u0map, parammap; check_length,
        t = tspan !== nothing ? tspan[1] : tspan,
        kwargs...)
    linenumbers = get(kwargs, :linenumbers, true)
    kwargs = filter_kwargs(kwargs)
    kwarg_params = gen_quoted_kwargs(kwargs)
    odep = Expr(:call, :ODEProblem, kwarg_params, :f, :u0, :tspan, :p)
    ex = quote
        f = $f
        u0 = $u0
        tspan = $tspan
        p = $p
        $odep
    end
    !linenumbers ? Base.remove_linenums!(ex) : ex
end

function ODEProblemExpr(sys::AbstractODESystem, args...; kwargs...)
    ODEProblemExpr{true}(sys, args...; kwargs...)
end

"""
```julia
DAEProblemExpr{iip}(sys::AbstractODESystem, u0map, tspan,
                    parammap = DiffEqBase.NullParameters();
                    version = nothing, tgrad = false,
                    jac = false,
                    checkbounds = false, sparse = false,
                    linenumbers = true, parallel = SerialForm(),
                    skipzeros = true, fillzeros = true,
                    simplify = false,
                    kwargs...) where {iip}
```

Generates a Julia expression for constructing a DAEProblem from an
ODESystem and allows for automatically symbolically calculating
numerical enhancements.
"""
struct DAEProblemExpr{iip} end

function DAEProblemExpr{iip}(sys::AbstractODESystem, du0map, u0map, tspan,
        parammap = DiffEqBase.NullParameters(); check_length = true,
        kwargs...) where {iip}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system before creating a `DAEProblemExpr`")
    end
    f, du0, u0, p = process_SciMLProblem(DAEFunctionExpr{iip}, sys, u0map, parammap;
        t = tspan !== nothing ? tspan[1] : tspan,
        implicit_dae = true, du0map = du0map, check_length,
        kwargs...)
    linenumbers = get(kwargs, :linenumbers, true)
    diffvars = collect_differential_variables(sys)
    sts = unknowns(sys)
    differential_vars = map(Base.Fix2(in, diffvars), sts)
    kwargs = filter_kwargs(kwargs)
    kwarg_params = gen_quoted_kwargs(kwargs)
    push!(kwarg_params.args, Expr(:kw, :differential_vars, :differential_vars))
    prob = Expr(:call, :(DAEProblem{$iip}), kwarg_params, :f, :du0, :u0, :tspan, :p)
    ex = quote
        f = $f
        u0 = $u0
        du0 = $du0
        tspan = $tspan
        p = $p
        differential_vars = $differential_vars
        $prob
    end
    !linenumbers ? Base.remove_linenums!(ex) : ex
end

function DAEProblemExpr(sys::AbstractODESystem, args...; kwargs...)
    DAEProblemExpr{true}(sys, args...; kwargs...)
end

"""
```julia
SciMLBase.SteadyStateProblem(sys::AbstractODESystem, u0map,
                             parammap = DiffEqBase.NullParameters();
                             version = nothing, tgrad = false,
                             jac = false,
                             checkbounds = false, sparse = false,
                             linenumbers = true, parallel = SerialForm(),
                             kwargs...) where {iip}
```

Generates an SteadyStateProblem from an ODESystem and allows for automatically
symbolically calculating numerical enhancements.
"""
function SciMLBase.SteadyStateProblem(sys::AbstractODESystem, args...; kwargs...)
    SteadyStateProblem{true}(sys, args...; kwargs...)
end

function DiffEqBase.SteadyStateProblem{iip}(sys::AbstractODESystem, u0map,
        parammap = SciMLBase.NullParameters();
        check_length = true, kwargs...) where {iip}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system before creating a `SteadyStateProblem`")
    end
    f, u0, p = process_SciMLProblem(ODEFunction{iip}, sys, u0map, parammap;
        steady_state = true,
        check_length, force_initialization_time_independent = true, kwargs...)
    kwargs = filter_kwargs(kwargs)
    SteadyStateProblem{iip}(f, u0, p; kwargs...)
end

"""
```julia
SciMLBase.SteadyStateProblemExpr(sys::AbstractODESystem, u0map,
                                 parammap = DiffEqBase.NullParameters();
                                 version = nothing, tgrad = false,
                                 jac = false,
                                 checkbounds = false, sparse = false,
                                 skipzeros = true, fillzeros = true,
                                 linenumbers = true, parallel = SerialForm(),
                                 kwargs...) where {iip}
```

Generates a Julia expression for building a SteadyStateProblem from
an ODESystem and allows for automatically symbolically calculating
numerical enhancements.
"""
struct SteadyStateProblemExpr{iip} end

function SteadyStateProblemExpr{iip}(sys::AbstractODESystem, u0map,
        parammap = SciMLBase.NullParameters();
        check_length = true,
        kwargs...) where {iip}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system before creating a `SteadyStateProblemExpr`")
    end
    f, u0, p = process_SciMLProblem(ODEFunctionExpr{iip}, sys, u0map, parammap;
        steady_state = true,
        check_length, kwargs...)
    linenumbers = get(kwargs, :linenumbers, true)
    kwargs = filter_kwargs(kwargs)
    kwarg_params = gen_quoted_kwargs(kwargs)
    prob = Expr(:call, :SteadyStateProblem, kwarg_params, :f, :u0, :p)
    ex = quote
        f = $f
        u0 = $u0
        p = $p
        $prob
    end
    !linenumbers ? Base.remove_linenums!(ex) : ex
end

function SteadyStateProblemExpr(sys::AbstractODESystem, args...; kwargs...)
    SteadyStateProblemExpr{true}(sys, args...; kwargs...)
end

function _match_eqs(eqs1, eqs2)
    eqpairs = Pair[]
    for (i, eq) in enumerate(eqs1)
        for (j, eq2) in enumerate(eqs2)
            if isequal(eq, eq2)
                push!(eqpairs, i => j)
                break
            end
        end
    end
    eqpairs
end

function isisomorphic(sys1::AbstractODESystem, sys2::AbstractODESystem)
    sys1 = flatten(sys1)
    sys2 = flatten(sys2)

    iv2 = only(independent_variables(sys2))
    sys1 = convert_system(ODESystem, sys1, iv2)
    s1, s2 = unknowns(sys1), unknowns(sys2)
    p1, p2 = parameters(sys1), parameters(sys2)

    (length(s1) != length(s2)) || (length(p1) != length(p2)) && return false

    eqs1 = equations(sys1)
    eqs2 = equations(sys2)

    pps = permutations(p2)
    psts = permutations(s2)
    orig = [p1; s1]
    perms = ([x; y] for x in pps for y in psts)

    for perm in perms
        rules = Dict(orig .=> perm)
        neweqs1 = substitute(eqs1, rules)
        eqpairs = _match_eqs(neweqs1, eqs2)
        if length(eqpairs) == length(eqs1)
            return true
        end
    end
    return false
end

function flatten_equations(eqs)
    mapreduce(vcat, eqs; init = Equation[]) do eq
        islhsarr = eq.lhs isa AbstractArray || Symbolics.isarraysymbolic(eq.lhs)
        isrhsarr = eq.rhs isa AbstractArray || Symbolics.isarraysymbolic(eq.rhs)
        if islhsarr || isrhsarr
            islhsarr && isrhsarr ||
                error("LHS ($(eq.lhs)) and RHS ($(eq.rhs)) must either both be array expressions or both scalar")
            size(eq.lhs) == size(eq.rhs) ||
                error("Size of LHS ($(eq.lhs)) and RHS ($(eq.rhs)) must match: got $(size(eq.lhs)) and $(size(eq.rhs))")
            return vec(collect(eq.lhs) .~ collect(eq.rhs))
        else
            eq
        end
    end
end

struct InitializationProblem{iip, specialization} end

"""
```julia
InitializationProblem{iip}(sys::AbstractODESystem, t, u0map,
                           parammap = DiffEqBase.NullParameters();
                           version = nothing, tgrad = false,
                           jac = false,
                           checkbounds = false, sparse = false,
                           simplify = false,
                           linenumbers = true, parallel = SerialForm(),
                           initialization_eqs = [],
                           fully_determined = false,
                           kwargs...) where {iip}
```

Generates a NonlinearProblem or NonlinearLeastSquaresProblem from an ODESystem
which represents the initialization, i.e. the calculation of the consistent
initial conditions for the given DAE.
"""
function InitializationProblem(sys::AbstractSystem, args...; kwargs...)
    InitializationProblem{true}(sys, args...; kwargs...)
end

function InitializationProblem(sys::AbstractSystem, t,
        u0map::StaticArray,
        args...;
        kwargs...)
    InitializationProblem{false, SciMLBase.FullSpecialize}(
        sys, t, u0map, args...; kwargs...)
end

function InitializationProblem{true}(sys::AbstractSystem, args...; kwargs...)
    InitializationProblem{true, SciMLBase.AutoSpecialize}(sys, args...; kwargs...)
end

function InitializationProblem{false}(sys::AbstractSystem, args...; kwargs...)
    InitializationProblem{false, SciMLBase.FullSpecialize}(sys, args...; kwargs...)
end

const INCOMPLETE_INITIALIZATION_MESSAGE = """
                                Initialization incomplete. Not all of the state variables of the
                                DAE system can be determined by the initialization. Missing
                                variables:
                                """

struct IncompleteInitializationError <: Exception
    uninit::Any
end

function Base.showerror(io::IO, e::IncompleteInitializationError)
    println(io, INCOMPLETE_INITIALIZATION_MESSAGE)
    println(io, e.uninit)
end

function InitializationProblem{iip, specialize}(sys::AbstractSystem,
        t, u0map = [],
        parammap = DiffEqBase.NullParameters();
        guesses = [],
        check_length = true,
        warn_initialize_determined = true,
        initialization_eqs = [],
        fully_determined = nothing,
        check_units = true,
        use_scc = true,
        allow_incomplete = false,
        force_time_independent = false,
        algebraic_only = false,
        kwargs...) where {iip, specialize}
    if !iscomplete(sys)
        error("A completed system is required. Call `complete` or `structural_simplify` on the system before creating an `ODEProblem`")
    end
    if isempty(u0map) && get_initializesystem(sys) !== nothing
        isys = get_initializesystem(sys; initialization_eqs, check_units)
        simplify_system = false
    elseif isempty(u0map) && get_initializesystem(sys) === nothing
        isys = generate_initializesystem(
            sys; initialization_eqs, check_units, pmap = parammap,
            guesses, extra_metadata = (; use_scc), algebraic_only)
        simplify_system = true
    else
        isys = generate_initializesystem(
            sys; u0map, initialization_eqs, check_units,
            pmap = parammap, guesses, extra_metadata = (; use_scc), algebraic_only)
        simplify_system = true
    end

    # useful for `SteadyStateProblem` since `f` has to be autonomous and the
    # initialization should be too
    if force_time_independent
        idx = findfirst(isequal(get_iv(sys)), get_ps(isys))
        idx === nothing || deleteat!(get_ps(isys), idx)
    end

    if simplify_system
        isys = structural_simplify(isys; fully_determined)
    end

    meta = get_metadata(isys)
    if meta isa InitializationSystemMetadata
        @set! isys.metadata.oop_reconstruct_u0_p = ReconstructInitializeprob(
            sys, isys)
    end

    ts = get_tearing_state(isys)
    unassigned_vars = StructuralTransformations.singular_check(ts)
    if warn_initialize_determined && !isempty(unassigned_vars)
        errmsg = """
        The initialization system is structurally singular. Guess values may \
        significantly affect the initial values of the ODE. The problematic variables \
        are $unassigned_vars.

        Note that the identification of problematic variables is a best-effort heuristic.
        """
        @warn errmsg
    end

    uninit = setdiff(unknowns(sys), [unknowns(isys); getfield.(observed(isys), :lhs)])

    # TODO: throw on uninitialized arrays
    filter!(x -> !(x isa Symbolics.Arr), uninit)
    if is_time_dependent(sys) && !isempty(uninit)
        allow_incomplete || throw(IncompleteInitializationError(uninit))
        # for incomplete initialization, we will add the missing variables as parameters.
        # they will be updated by `update_initializeprob!` and `initializeprobmap` will
        # use them to construct the new `u0`.
        newparams = map(toparam, uninit)
        append!(get_ps(isys), newparams)
        isys = complete(isys)
    end

    neqs = length(equations(isys))
    nunknown = length(unknowns(isys))

    if use_scc
        scc_message = "`SCCNonlinearProblem` can only be used for initialization of fully determined systems and hence will not be used here. "
    else
        scc_message = ""
    end

    if warn_initialize_determined && neqs > nunknown
        @warn "Initialization system is overdetermined. $neqs equations for $nunknown unknowns. Initialization will default to using least squares. $(scc_message)To suppress this warning pass warn_initialize_determined = false. To make this warning into an error, pass fully_determined = true"
    end
    if warn_initialize_determined && neqs < nunknown
        @warn "Initialization system is underdetermined. $neqs equations for $nunknown unknowns. Initialization will default to using least squares. $(scc_message)To suppress this warning pass warn_initialize_determined = false. To make this warning into an error, pass fully_determined = true"
    end

    parammap = recursive_unwrap(anydict(parammap))
    if t !== nothing
        parammap[get_iv(sys)] = t
    end
    filter!(kvp -> kvp[2] !== missing, parammap)

    u0map = to_varmap(u0map, unknowns(sys))
    if isempty(guesses)
        guesses = Dict()
    end

    filter_missing_values!(u0map)
    filter_missing_values!(parammap)
    u0map = merge(ModelingToolkit.guesses(sys), todict(guesses), u0map)

    fullmap = merge(u0map, parammap)
    u0T = Union{}
    for sym in unknowns(isys)
        val = fixpoint_sub(sym, fullmap)
        symbolic_type(val) == NotSymbolic() || continue
        u0T = promote_type(u0T, typeof(val))
    end
    for eq in observed(isys)
        # ignore HACK-ed observed equations
        symbolic_type(eq.lhs) == ArraySymbolic() && continue
        val = fixpoint_sub(eq.lhs, fullmap)
        symbolic_type(val) == NotSymbolic() || continue
        u0T = promote_type(u0T, typeof(val))
    end
    if u0T != Union{}
        u0T = eltype(u0T)
        u0map = Dict(k => if v === nothing
                         nothing
                     elseif symbolic_type(v) == NotSymbolic() && !is_array_of_symbolics(v)
                         v isa AbstractArray ? u0T.(v) : u0T(v)
                     else
                         v
                     end
        for (k, v) in u0map)
    end

    TProb = if neqs == nunknown && isempty(unassigned_vars)
        if use_scc && neqs > 0
            if is_split(isys)
                SCCNonlinearProblem
            else
                @warn "`SCCNonlinearProblem` can only be used with `split = true` systems. Simplify your `ODESystem` with `split = true` or pass `use_scc = false` to disable this warning"
                NonlinearProblem
            end
        else
            NonlinearProblem
        end
    else
        NonlinearLeastSquaresProblem
    end
    TProb(isys, u0map, parammap; kwargs...,
        build_initializeprob = false, is_initializeprob = true)
end
