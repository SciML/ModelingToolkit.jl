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

    # We need to remove explicit time dependence on the state because when we
    # have `u(t) * t` we want to have the tgrad to be `u(t)` instead of `u'(t) *
    # t + u(t)`.
    rhs = [detime_dvs(eq.rhs) for eq in full_equations(sys)]
    iv = get_iv(sys)
    xs = states(sys)
    rule = Dict(map((x, xt) -> xt => x, detime_dvs.(xs), xs))
    rhs = substitute.(rhs, Ref(rule))
    tgrad = [expand_derivatives(Differential(iv)(r), simplify) for r in rhs]
    reverse_rule = Dict(map((x, xt) -> x => xt, detime_dvs.(xs), xs))
    tgrad = Num.(substitute.(tgrad, Ref(reverse_rule)))
    get_tgrad(sys)[] = tgrad
    return tgrad
end

function calculate_jacobian(sys::AbstractODESystem;
                            sparse = false, simplify = false, dvs = states(sys))
    if isequal(dvs, states(sys))
        cache = get_jac(sys)[]
        if cache isa Tuple && cache[2] == (sparse, simplify)
            return cache[1]
        end
    end

    rhs = [eq.rhs - eq.lhs for eq in full_equations(sys)] #need du terms on rhs for differentiating wrt du

    iv = get_iv(sys)

    if sparse
        jac = sparsejacobian(rhs, dvs, simplify = simplify)
    else
        jac = jacobian(rhs, dvs, simplify = simplify)
    end

    if isequal(dvs, states(sys))
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

    iv = get_iv(sys)
    ctrls = controls(sys)

    if sparse
        jac = sparsejacobian(rhs, ctrls, simplify = simplify)
    else
        jac = jacobian(rhs, ctrls, simplify = simplify)
    end

    get_ctrl_jac(sys)[] = jac, (sparse, simplify) # cache Jacobian
    return jac
end

function generate_tgrad(sys::AbstractODESystem, dvs = states(sys), ps = parameters(sys);
                        simplify = false, kwargs...)
    tgrad = calculate_tgrad(sys, simplify = simplify)
    return build_function(tgrad, dvs, ps, get_iv(sys); kwargs...)
end

function generate_jacobian(sys::AbstractODESystem, dvs = states(sys), ps = parameters(sys);
                           simplify = false, sparse = false, kwargs...)
    jac = calculate_jacobian(sys; simplify = simplify, sparse = sparse)
    return build_function(jac, dvs, ps, get_iv(sys); kwargs...)
end

function generate_control_jacobian(sys::AbstractODESystem, dvs = states(sys),
                                   ps = parameters(sys);
                                   simplify = false, sparse = false, kwargs...)
    jac = calculate_control_jacobian(sys; simplify = simplify, sparse = sparse)
    return build_function(jac, dvs, ps, get_iv(sys); kwargs...)
end

function generate_dae_jacobian(sys::AbstractODESystem, dvs = states(sys),
                               ps = parameters(sys); simplify = false, sparse = false,
                               kwargs...)
    jac_u = calculate_jacobian(sys; simplify = simplify, sparse = sparse)
    derivatives = Differential(get_iv(sys)).(states(sys))
    jac_du = calculate_jacobian(sys; simplify = simplify, sparse = sparse,
                                dvs = derivatives)
    dvs = states(sys)
    @variables ˍ₋gamma
    jac = ˍ₋gamma * jac_du + jac_u
    return build_function(jac, derivatives, dvs, ps, ˍ₋gamma, get_iv(sys); kwargs...)
end

function generate_function(sys::AbstractODESystem, dvs = states(sys), ps = parameters(sys);
                           implicit_dae = false,
                           ddvs = implicit_dae ? map(Differential(get_iv(sys)), dvs) :
                                  nothing,
                           has_difference = false,
                           kwargs...)
    eqs = [eq for eq in equations(sys) if !isdifferenceeq(eq)]
    if !implicit_dae
        check_operator_variables(eqs, Differential)
        check_lhs(eqs, Differential, Set(dvs))
    end
    # substitute x(t) by just x
    rhss = implicit_dae ? [_iszero(eq.lhs) ? eq.rhs : eq.rhs - eq.lhs for eq in eqs] :
           [eq.rhs for eq in eqs]

    # TODO: add an optional check on the ordering of observed equations
    u = map(x -> time_varying_as_func(value(x), sys), dvs)
    p = map(x -> time_varying_as_func(value(x), sys), ps)
    t = get_iv(sys)

    pre, sol_states = get_substitutions_and_solved_states(sys,
                                                          no_postprocess = has_difference)

    if implicit_dae
        build_function(rhss, ddvs, u, p, t; postprocess_fbody = pre, states = sol_states,
                       kwargs...)
    else
        build_function(rhss, u, p, t; postprocess_fbody = pre, states = sol_states,
                       kwargs...)
    end
end

function generate_difference_cb(sys::ODESystem, dvs = states(sys), ps = parameters(sys);
                                kwargs...)
    eqs = equations(sys)
    check_operator_variables(eqs, Difference)

    var2eq = Dict(arguments(eq.lhs)[1] => eq for eq in eqs if isdifference(eq.lhs))

    u = map(x -> time_varying_as_func(value(x), sys), dvs)
    p = map(x -> time_varying_as_func(value(x), sys), ps)
    t = get_iv(sys)

    body = map(dvs) do v
        eq = get(var2eq, v, nothing)
        eq === nothing && return v
        d = operation(eq.lhs)
        d.update ? eq.rhs : eq.rhs + v
    end

    pre = get_postprocess_fbody(sys)
    f_oop, f_iip = build_function(body, u, p, t; expression = Val{false},
                                  postprocess_fbody = pre, kwargs...)

    cb_affect! = let f_oop = f_oop, f_iip = f_iip
        function cb_affect!(integ)
            if DiffEqBase.isinplace(integ.sol.prob)
                tmp, = DiffEqBase.get_tmp_cache(integ)
                f_iip(tmp, integ.u, integ.p, integ.t) # aliasing `integ.u` would be bad.
                copyto!(integ.u, tmp)
            else
                integ.u = f_oop(integ.u, integ.p, integ.t)
            end
            return nothing
        end
    end

    getdt(eq) = operation(eq.lhs).dt
    deqs = values(var2eq)
    dt = getdt(first(deqs))
    all(dt == getdt(eq) for eq in deqs) ||
        error("All difference variables should have same time steps.")

    PeriodicCallback(cb_affect!, first(dt))
end

function calculate_massmatrix(sys::AbstractODESystem; simplify = false)
    eqs = [eq for eq in full_equations(sys) if !isdifferenceeq(eq)]
    dvs = states(sys)
    M = zeros(length(eqs), length(eqs))
    state2idx = Dict(s => i for (i, s) in enumerate(dvs))
    for (i, eq) in enumerate(eqs)
        if eq.lhs isa Term && operation(eq.lhs) isa Differential
            st = var_from_nested_derivative(eq.lhs)[1]
            j = state2idx[st]
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
                      [dv for dv in states(sys)])
end

function jacobian_dae_sparsity(sys::AbstractODESystem)
    J1 = jacobian_sparsity([eq.rhs for eq in full_equations(sys)],
                           [dv for dv in states(sys)])
    derivatives = Differential(get_iv(sys)).(states(sys))
    J2 = jacobian_sparsity([eq.rhs for eq in full_equations(sys)],
                           [dv for dv in derivatives])
    J1 + J2
end

function isautonomous(sys::AbstractODESystem)
    tgrad = calculate_tgrad(sys; simplify = true)
    all(iszero, tgrad)
end

"""
```julia
function DiffEqBase.ODEFunction{iip}(sys::AbstractODESystem, dvs = states(sys),
                                     ps = parameters(sys);
                                     version = nothing, tgrad=false,
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

function DiffEqBase.ODEFunction{iip, specialize}(sys::AbstractODESystem, dvs = states(sys),
                                                 ps = parameters(sys), u0 = nothing;
                                                 version = nothing, tgrad = false,
                                                 jac = false, p = nothing,
                                                 t = nothing,
                                                 eval_expression = true,
                                                 sparse = false, simplify = false,
                                                 eval_module = @__MODULE__,
                                                 steady_state = false,
                                                 checkbounds = false,
                                                 sparsity = false,
                                                 kwargs...) where {iip, specialize}
    f_gen = generate_function(sys, dvs, ps; expression = Val{eval_expression},
                              expression_module = eval_module, checkbounds = checkbounds,
                              kwargs...)
    f_oop, f_iip = eval_expression ?
                   (@RuntimeGeneratedFunction(eval_module, ex) for ex in f_gen) : f_gen
    f(u, p, t) = f_oop(u, p, t)
    f(du, u, p, t) = f_iip(du, u, p, t)

    if specialize === SciMLBase.FunctionWrapperSpecialize && iip
        if u0 === nothing || p === nothing || t === nothing
            error("u0, p, and t must be specified for FunctionWrapperSpecialize on ODEFunction.")
        end
        f = SciMLBase.wrapfun_iip(f, (u0, u0, p, t))
    end

    if tgrad
        tgrad_gen = generate_tgrad(sys, dvs, ps;
                                   simplify = simplify,
                                   expression = Val{eval_expression},
                                   expression_module = eval_module,
                                   checkbounds = checkbounds, kwargs...)
        tgrad_oop, tgrad_iip = eval_expression ?
                               (@RuntimeGeneratedFunction(eval_module, ex) for ex in tgrad_gen) :
                               tgrad_gen
        _tgrad(u, p, t) = tgrad_oop(u, p, t)
        _tgrad(J, u, p, t) = tgrad_iip(J, u, p, t)
    else
        _tgrad = nothing
    end

    if jac
        jac_gen = generate_jacobian(sys, dvs, ps;
                                    simplify = simplify, sparse = sparse,
                                    expression = Val{eval_expression},
                                    expression_module = eval_module,
                                    checkbounds = checkbounds, kwargs...)
        jac_oop, jac_iip = eval_expression ?
                           (@RuntimeGeneratedFunction(eval_module, ex) for ex in jac_gen) :
                           jac_gen
        _jac(u, p, t) = jac_oop(u, p, t)
        _jac(J, u, p, t) = jac_iip(J, u, p, t)
    else
        _jac = nothing
    end

    M = calculate_massmatrix(sys)

    _M = if sparse && !(u0 === nothing || M === I)
        SparseArrays.sparse(M)
    elseif u0 === nothing || M === I
        M
    else
        ArrayInterfaceCore.restructure(u0 .* u0', M)
    end

    obs = observed(sys)
    observedfun = if steady_state
        let sys = sys, dict = Dict()
            function generated_observed(obsvar, u, p, t = Inf)
                obs = get!(dict, value(obsvar)) do
                    build_explicit_observed_function(sys, obsvar)
                end
                obs(u, p, t)
            end
        end
    else
        let sys = sys, dict = Dict()
            function generated_observed(obsvar, u, p, t)
                obs = get!(dict, value(obsvar)) do
                    build_explicit_observed_function(sys, obsvar; checkbounds = checkbounds)
                end
                obs(u, p, t)
            end
        end
    end

    jac_prototype = if sparse
        uElType = u0 === nothing ? Float64 : eltype(u0)
        if jac
            similar(calculate_jacobian(sys, sparse = sparse), uElType)
        else
            similar(jacobian_sparsity(sys), uElType)
        end
    else
        nothing
    end
    ODEFunction{iip, specialize}(f,
                                 sys = sys,
                                 jac = _jac === nothing ? nothing : _jac,
                                 tgrad = _tgrad === nothing ? nothing : _tgrad,
                                 mass_matrix = _M,
                                 jac_prototype = jac_prototype,
                                 syms = Symbol.(states(sys)),
                                 indepsym = Symbol(get_iv(sys)),
                                 observed = observedfun,
                                 sparsity = sparsity ? jacobian_sparsity(sys) : nothing)
end

"""
```julia
function DiffEqBase.DAEFunction{iip}(sys::AbstractODESystem, dvs = states(sys),
                                     ps = parameters(sys);
                                     version = nothing, tgrad=false,
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

function DiffEqBase.DAEFunction{iip}(sys::AbstractODESystem, dvs = states(sys),
                                     ps = parameters(sys), u0 = nothing;
                                     ddvs = map(diff2term ∘ Differential(get_iv(sys)), dvs),
                                     version = nothing, p = nothing,
                                     jac = false,
                                     eval_expression = true,
                                     sparse = false, simplify = false,
                                     eval_module = @__MODULE__,
                                     checkbounds = false,
                                     kwargs...) where {iip}
    f_gen = generate_function(sys, dvs, ps; implicit_dae = true,
                              expression = Val{eval_expression},
                              expression_module = eval_module, checkbounds = checkbounds,
                              kwargs...)
    f_oop, f_iip = eval_expression ?
                   (@RuntimeGeneratedFunction(eval_module, ex) for ex in f_gen) : f_gen
    f(du, u, p, t) = f_oop(du, u, p, t)
    f(out, du, u, p, t) = f_iip(out, du, u, p, t)

    if jac
        jac_gen = generate_dae_jacobian(sys, dvs, ps;
                                        simplify = simplify, sparse = sparse,
                                        expression = Val{eval_expression},
                                        expression_module = eval_module,
                                        checkbounds = checkbounds, kwargs...)
        jac_oop, jac_iip = eval_expression ?
                           (@RuntimeGeneratedFunction(eval_module, ex) for ex in jac_gen) :
                           jac_gen
        _jac(du, u, p, ˍ₋gamma, t) = jac_oop(du, u, p, ˍ₋gamma, t)

        _jac(J, du, u, p, ˍ₋gamma, t) = jac_iip(J, du, u, p, ˍ₋gamma, t)
    else
        _jac = nothing
    end

    obs = observed(sys)
    observedfun = let sys = sys, dict = Dict()
        function generated_observed(obsvar, u, p, t)
            obs = get!(dict, value(obsvar)) do
                build_explicit_observed_function(sys, obsvar; checkbounds = checkbounds)
            end
            obs(u, p, t)
        end
    end

    jac_prototype = if sparse
        uElType = u0 === nothing ? Float64 : eltype(u0)
        if jac
            J1 = calculate_jacobian(sys, sparse = sparse)
            derivatives = Differential(get_iv(sys)).(states(sys))
            J2 = calculate_jacobian(sys; sparse = sparse, dvs = derivatives)
            similar(J1 + J2, uElType)
        else
            similar(jacobian_dae_sparsity(sys), uElType)
        end
    else
        nothing
    end

    DAEFunction{iip}(f,
                     sys = sys,
                     jac = _jac === nothing ? nothing : _jac,
                     syms = Symbol.(dvs),
                     jac_prototype = jac_prototype,
                     # missing fields in `DAEFunction`
                     #indepsym = Symbol(get_iv(sys)),
                     observed = observedfun)
end

"""
```julia
function ODEFunctionExpr{iip}(sys::AbstractODESystem, dvs = states(sys),
                                     ps = parameters(sys);
                                     version = nothing, tgrad=false,
                                     jac = false,
                                     sparse = false,
                                     kwargs...) where {iip}
```

Create a Julia expression for an `ODEFunction` from the [`ODESystem`](@ref).
The arguments `dvs` and `ps` are used to set the order of the dependent
variable and parameter vectors, respectively.
"""
struct ODEFunctionExpr{iip} end

struct ODEFunctionClosure{O, I} <: Function
    f_oop::O
    f_iip::I
end
(f::ODEFunctionClosure)(u, p, t) = f.f_oop(u, p, t)
(f::ODEFunctionClosure)(du, u, p, t) = f.f_iip(du, u, p, t)

function ODEFunctionExpr{iip}(sys::AbstractODESystem, dvs = states(sys),
                              ps = parameters(sys), u0 = nothing;
                              version = nothing, tgrad = false,
                              jac = false, p = nothing,
                              linenumbers = false,
                              sparse = false, simplify = false,
                              steady_state = false,
                              kwargs...) where {iip}
    f_oop, f_iip = generate_function(sys, dvs, ps; expression = Val{true}, kwargs...)

    dict = Dict()

    fsym = gensym(:f)
    _f = :($fsym = $ODEFunctionClosure($f_oop, $f_iip))
    tgradsym = gensym(:tgrad)
    if tgrad
        tgrad_oop, tgrad_iip = generate_tgrad(sys, dvs, ps;
                                              simplify = simplify,
                                              expression = Val{true}, kwargs...)
        _tgrad = :($tgradsym = $ODEFunctionClosure($tgrad_oop, $tgrad_iip))
    else
        _tgrad = :($tgradsym = nothing)
    end

    jacsym = gensym(:jac)
    if jac
        jac_oop, jac_iip = generate_jacobian(sys, dvs, ps;
                                             sparse = sparse, simplify = simplify,
                                             expression = Val{true}, kwargs...)
        _jac = :($jacsym = $ODEFunctionClosure($jac_oop, $jac_iip))
    else
        _jac = :($jacsym = nothing)
    end

    M = calculate_massmatrix(sys)

    _M = if sparse && !(u0 === nothing || M === I)
        SparseArrays.sparse(M)
    elseif u0 === nothing || M === I
        M
    else
        ArrayInterfaceCore.restructure(u0 .* u0', M)
    end

    jp_expr = sparse ? :(similar($(get_jac(sys)[]), Float64)) : :nothing
    ex = quote
        $_f
        $_tgrad
        $_jac
        M = $_M
        ODEFunction{$iip}($fsym,
                          jac = $jacsym,
                          tgrad = $tgradsym,
                          mass_matrix = M,
                          jac_prototype = $jp_expr,
                          syms = $(Symbol.(states(sys))),
                          indepsym = $(QuoteNode(Symbol(get_iv(sys)))))
    end
    !linenumbers ? striplines(ex) : ex
end

function process_DEProblem(constructor, sys::AbstractODESystem, u0map, parammap;
                           implicit_dae = false, du0map = nothing,
                           version = nothing, tgrad = false,
                           jac = false,
                           checkbounds = false, sparse = false,
                           simplify = false,
                           linenumbers = true, parallel = nothing,
                           eval_expression = true,
                           use_union = false,
                           kwargs...)
    eqs = equations(sys)
    dvs = states(sys)
    ps = parameters(sys)
    iv = get_iv(sys)

    defs = defaults(sys)
    defs = mergedefaults(defs, parammap, ps)
    defs = mergedefaults(defs, u0map, dvs)

    u0 = varmap_to_vars(u0map, dvs; defaults = defs, tofloat = true)
    p = varmap_to_vars(parammap, ps; defaults = defs, tofloat = !use_union, use_union)
    p = p === nothing ? SciMLBase.NullParameters() : p

    if implicit_dae && du0map !== nothing
        ddvs = map(Differential(iv), dvs)
        defs = mergedefaults(defs, du0map, ddvs)
        du0 = varmap_to_vars(du0map, ddvs; defaults = defs, toterm = identity,
                             tofloat = true)
    else
        du0 = nothing
        ddvs = nothing
    end

    check_eqs_u0(eqs, dvs, u0; kwargs...)

    f = constructor(sys, dvs, ps, u0; ddvs = ddvs, tgrad = tgrad, jac = jac,
                    checkbounds = checkbounds, p = p,
                    linenumbers = linenumbers, parallel = parallel, simplify = simplify,
                    sparse = sparse, eval_expression = eval_expression, kwargs...)
    implicit_dae ? (f, du0, u0, p) : (f, u0, p)
end

function ODEFunctionExpr(sys::AbstractODESystem, args...; kwargs...)
    ODEFunctionExpr{true}(sys, args...; kwargs...)
end

"""
```julia
function DAEFunctionExpr{iip}(sys::AbstractODESystem, dvs = states(sys),
                                     ps = parameters(sys);
                                     version = nothing, tgrad=false,
                                     jac = false,
                                     sparse = false,
                                     kwargs...) where {iip}
```

Create a Julia expression for an `ODEFunction` from the [`ODESystem`](@ref).
The arguments `dvs` and `ps` are used to set the order of the dependent
variable and parameter vectors, respectively.
"""
struct DAEFunctionExpr{iip} end

struct DAEFunctionClosure{O, I} <: Function
    f_oop::O
    f_iip::I
end
(f::DAEFunctionClosure)(du, u, p, t) = f.f_oop(du, u, p, t)
(f::DAEFunctionClosure)(out, du, u, p, t) = f.f_iip(out, du, u, p, t)

function DAEFunctionExpr{iip}(sys::AbstractODESystem, dvs = states(sys),
                              ps = parameters(sys), u0 = nothing;
                              version = nothing, tgrad = false,
                              jac = false, p = nothing,
                              linenumbers = false,
                              sparse = false, simplify = false,
                              kwargs...) where {iip}
    f_oop, f_iip = generate_function(sys, dvs, ps; expression = Val{true},
                                     implicit_dae = true, kwargs...)
    fsym = gensym(:f)
    _f = :($fsym = $DAEFunctionClosure($f_oop, $f_iip))
    ex = quote
        $_f
        ODEFunction{$iip}($fsym)
    end
    !linenumbers ? striplines(ex) : ex
end

function DAEFunctionExpr(sys::AbstractODESystem, args...; kwargs...)
    DAEFunctionExpr{true}(sys, args...; kwargs...)
end

"""
```julia
function DiffEqBase.ODEProblem{iip}(sys::AbstractODESystem,u0map,tspan,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false,
                                    checkbounds = false, sparse = false,
                                    simplify=false,
                                    linenumbers = true, parallel=nothing,
                                    kwargs...) where iip
```

Generates an ODEProblem from an ODESystem and allows for automatically
symbolically calculating numerical enhancements.
"""
function DiffEqBase.ODEProblem(sys::AbstractODESystem, args...; kwargs...)
    ODEProblem{true}(sys, args...; kwargs...)
end

function DiffEqBase.ODEProblem{true}(sys::AbstractODESystem, args...; kwargs...)
    ODEProblem{true, SciMLBase.AutoSpecialize}(sys, args...; kwargs...)
end

function DiffEqBase.ODEProblem{false}(sys::AbstractODESystem, args...; kwargs...)
    ODEProblem{false, SciMLBase.FullSpecialize}(sys, args...; kwargs...)
end

function DiffEqBase.ODEProblem{iip, specialize}(sys::AbstractODESystem, u0map, tspan,
                                                parammap = DiffEqBase.NullParameters();
                                                callback = nothing,
                                                check_length = true,
                                                kwargs...) where {iip, specialize}
    has_difference = any(isdifferenceeq, equations(sys))
    f, u0, p = process_DEProblem(ODEFunction{iip, specialize}, sys, u0map, parammap;
                                 t = tspan !== nothing ? tspan[1] : tspan,
                                 has_difference = has_difference,
                                 check_length, kwargs...)
    cbs = process_events(sys; callback, has_difference, kwargs...)
    kwargs = filter_kwargs(kwargs)
    pt = something(get_metadata(sys), StandardODEProblem())

    if cbs === nothing
        ODEProblem{iip}(f, u0, tspan, p, pt; kwargs...)
    else
        ODEProblem{iip}(f, u0, tspan, p, pt; callback = cbs, kwargs...)
    end
end
get_callback(prob::ODEProblem) = prob.kwargs[:callback]

"""
```julia
function DiffEqBase.DAEProblem{iip}(sys::AbstractODESystem,du0map,u0map,tspan,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false,
                                    checkbounds = false, sparse = false,
                                    simplify=false,
                                    linenumbers = true, parallel=nothing,
                                    kwargs...) where iip
```

Generates a DAEProblem from an ODESystem and allows for automatically
symbolically calculating numerical enhancements.
"""
function DiffEqBase.DAEProblem(sys::AbstractODESystem, args...; kwargs...)
    DAEProblem{true}(sys, args...; kwargs...)
end

function DiffEqBase.DAEProblem{iip}(sys::AbstractODESystem, du0map, u0map, tspan,
                                    parammap = DiffEqBase.NullParameters();
                                    check_length = true, kwargs...) where {iip}
    has_difference = any(isdifferenceeq, equations(sys))
    f, du0, u0, p = process_DEProblem(DAEFunction{iip}, sys, u0map, parammap;
                                      implicit_dae = true, du0map = du0map,
                                      has_difference = has_difference, check_length,
                                      kwargs...)
    diffvars = collect_differential_variables(sys)
    sts = states(sys)
    differential_vars = map(Base.Fix2(in, diffvars), sts)
    kwargs = filter_kwargs(kwargs)

    if has_difference
        DAEProblem{iip}(f, du0, u0, tspan, p;
                        difference_cb = generate_difference_cb(sys; kwargs...),
                        differential_vars = differential_vars, kwargs...)
    else
        DAEProblem{iip}(f, du0, u0, tspan, p; differential_vars = differential_vars,
                        kwargs...)
    end
end

"""
```julia
function ODEProblemExpr{iip}(sys::AbstractODESystem,u0map,tspan,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false,
                                    checkbounds = false, sparse = false,
                                    linenumbers = true, parallel=nothing,
                                    skipzeros=true, fillzeros=true,
                                    simplify=false,
                                    kwargs...) where iip
```

Generates a Julia expression for constructing an ODEProblem from an
ODESystem and allows for automatically symbolically calculating
numerical enhancements.
"""
struct ODEProblemExpr{iip} end

function ODEProblemExpr{iip}(sys::AbstractODESystem, u0map, tspan,
                             parammap = DiffEqBase.NullParameters(); check_length = true,
                             kwargs...) where {iip}
    f, u0, p = process_DEProblem(ODEFunctionExpr{iip}, sys, u0map, parammap; check_length,
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
    !linenumbers ? striplines(ex) : ex
end

function ODEProblemExpr(sys::AbstractODESystem, args...; kwargs...)
    ODEProblemExpr{true}(sys, args...; kwargs...)
end

"""
```julia
function DAEProblemExpr{iip}(sys::AbstractODESystem,u0map,tspan,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false,
                                    checkbounds = false, sparse = false,
                                    linenumbers = true, parallel=nothing,
                                    skipzeros=true, fillzeros=true,
                                    simplify=false,
                                    kwargs...) where iip
```

Generates a Julia expression for constructing a DAEProblem from an
ODESystem and allows for automatically symbolically calculating
numerical enhancements.
"""
struct DAEProblemExpr{iip} end

function DAEProblemExpr{iip}(sys::AbstractODESystem, du0map, u0map, tspan,
                             parammap = DiffEqBase.NullParameters(); check_length = true,
                             kwargs...) where {iip}
    f, du0, u0, p = process_DEProblem(DAEFunctionExpr{iip}, sys, u0map, parammap;
                                      implicit_dae = true, du0map = du0map, check_length,
                                      kwargs...)
    linenumbers = get(kwargs, :linenumbers, true)
    diffvars = collect_differential_variables(sys)
    sts = states(sys)
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
    !linenumbers ? striplines(ex) : ex
end

function DAEProblemExpr(sys::AbstractODESystem, args...; kwargs...)
    DAEProblemExpr{true}(sys, args...; kwargs...)
end

"""
```julia
function SciMLBase.SteadyStateProblem(sys::AbstractODESystem,u0map,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false,
                                    checkbounds = false, sparse = false,
                                    linenumbers = true, parallel=nothing,
                                    kwargs...) where iip
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
    f, u0, p = process_DEProblem(ODEFunction{iip}, sys, u0map, parammap;
                                 steady_state = true,
                                 check_length, kwargs...)
    kwargs = filter_kwargs(kwargs)
    SteadyStateProblem{iip}(f, u0, p; kwargs...)
end

"""
```julia
function SciMLBase.SteadyStateProblemExpr(sys::AbstractODESystem,u0map,
                                    parammap=DiffEqBase.NullParameters();
                                    version = nothing, tgrad=false,
                                    jac = false,
                                    checkbounds = false, sparse = false,
                                    skipzeros=true, fillzeros=true,
                                    linenumbers = true, parallel=nothing,
                                    kwargs...) where iip
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
    f, u0, p = process_DEProblem(ODEFunctionExpr{iip}, sys, u0map, parammap;
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
    !linenumbers ? striplines(ex) : ex
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
    s1, s2 = states(sys1), states(sys2)
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
