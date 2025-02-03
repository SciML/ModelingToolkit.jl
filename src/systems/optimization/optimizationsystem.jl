"""
$(TYPEDEF)

A scalar equation for optimization.

# Fields
$(FIELDS)

# Examples

```julia
@variables x y z
@parameters a b c

obj = a * (y - x) + x * (b - z) - y + x * y - c * z
cons = [x^2 + y^2 ≲ 1]
@named os = OptimizationSystem(obj, [x, y, z], [a, b, c]; constraints = cons)
```
"""
struct OptimizationSystem <: AbstractOptimizationSystem
    """
    A tag for the system. If two systems have the same tag, then they are
    structurally identical.
    """
    tag::UInt
    """Objective function of the system."""
    op::Any
    """Unknown variables."""
    unknowns::Array
    """Parameters."""
    ps::Vector
    """Array variables."""
    var_to_name::Any
    """Observed variables."""
    observed::Vector{Equation}
    """List of constraint equations of the system."""
    constraints::Vector{Union{Equation, Inequality}}
    """The name of the system."""
    name::Symbol
    """A description of the system."""
    description::String
    """The internal systems. These are required to have unique names."""
    systems::Vector{OptimizationSystem}
    """
    The default values to use when initial guess and/or
    parameters are not supplied in `OptimizationProblem`.
    """
    defaults::Dict
    """
    Metadata for the system, to be used by downstream packages.
    """
    metadata::Any
    """
    Metadata for MTK GUI.
    """
    gui_metadata::Union{Nothing, GUIMetadata}
    """
    If a model `sys` is complete, then `sys.x` no longer performs namespacing.
    """
    complete::Bool
    """
    Cached data for fast symbolic indexing.
    """
    index_cache::Union{Nothing, IndexCache}
    """
    The hierarchical parent system before simplification.
    """
    parent::Any
    isscheduled::Bool

    function OptimizationSystem(tag, op, unknowns, ps, var_to_name, observed,
            constraints, name, description, systems, defaults, metadata = nothing,
            gui_metadata = nothing, complete = false, index_cache = nothing, parent = nothing,
            isscheduled = false;
            checks::Union{Bool, Int} = true)
        if checks == true || (checks & CheckUnits) > 0
            u = __get_unit_type(unknowns, ps)
            unwrap(op) isa Symbolic && check_units(u, op)
            check_units(u, observed)
            check_units(u, constraints)
        end
        new(tag, op, unknowns, ps, var_to_name, observed,
            constraints, name, description, systems, defaults, metadata, gui_metadata, complete,
            index_cache, parent, isscheduled)
    end
end

equations(sys::AbstractOptimizationSystem) = objective(sys) # needed for Base.show

function OptimizationSystem(op, unknowns, ps;
        observed = [],
        constraints = [],
        default_u0 = Dict(),
        default_p = Dict(),
        defaults = _merge(Dict(default_u0), Dict(default_p)),
        name = nothing,
        description = "",
        systems = OptimizationSystem[],
        checks = true,
        metadata = nothing,
        gui_metadata = nothing)
    name === nothing &&
        throw(ArgumentError("The `name` keyword must be provided. Please consider using the `@named` macro"))
    constraints = value.(reduce(vcat, scalarize(constraints); init = []))
    unknowns′ = value.(reduce(vcat, scalarize(unknowns); init = []))
    ps′ = value.(ps)
    op′ = value(scalarize(op))

    irreducible_subs = Dict()
    for i in eachindex(unknowns′)
        var = unknowns′[i]
        if hasbounds(var)
            irreducible_subs[var] = irrvar = setirreducible(var, true)
            unknowns′[i] = irrvar
        end
    end
    op′ = substitute(op′, irreducible_subs)
    constraints = substitute.(constraints, (irreducible_subs,))

    if !(isempty(default_u0) && isempty(default_p))
        Base.depwarn(
            "`default_u0` and `default_p` are deprecated. Use `defaults` instead.",
            :OptimizationSystem, force = true)
    end
    sysnames = nameof.(systems)
    if length(unique(sysnames)) != length(sysnames)
        throw(ArgumentError("System names must be unique."))
    end
    defaults = todict(defaults)
    defaults = Dict(substitute(value(k), irreducible_subs) => substitute(
                        value(v), irreducible_subs)
    for (k, v) in pairs(defaults) if value(v) !== nothing)

    var_to_name = Dict()
    process_variables!(var_to_name, defaults, Dict(), unknowns′)
    process_variables!(var_to_name, defaults, Dict(), ps′)
    isempty(observed) || collect_var_to_name!(var_to_name, (eq.lhs for eq in observed))

    OptimizationSystem(Threads.atomic_add!(SYSTEM_COUNT, UInt(1)),
        op′, unknowns′, ps′, var_to_name,
        observed,
        constraints,
        name, description, systems, defaults, metadata, gui_metadata;
        checks = checks)
end

function OptimizationSystem(objective; constraints = [], kwargs...)
    allunknowns = OrderedSet()
    ps = OrderedSet()
    collect_vars!(allunknowns, ps, objective, nothing)
    for cons in constraints
        collect_vars!(allunknowns, ps, cons, nothing)
    end
    for ssys in get(kwargs, :systems, OptimizationSystem[])
        collect_scoped_vars!(allunknowns, ps, ssys, nothing)
    end
    new_ps = OrderedSet()
    for p in ps
        if iscall(p) && operation(p) === getindex
            par = arguments(p)[begin]
            if Symbolics.shape(Symbolics.unwrap(par)) !== Symbolics.Unknown() &&
               all(par[i] in ps for i in eachindex(par))
                push!(new_ps, par)
            else
                push!(new_ps, p)
            end
        else
            push!(new_ps, p)
        end
    end
    return OptimizationSystem(
        objective, collect(allunknowns), collect(new_ps); constraints, kwargs...)
end

function flatten(sys::OptimizationSystem)
    systems = get_systems(sys)
    isempty(systems) && return sys

    return OptimizationSystem(
        objective(sys),
        unknowns(sys),
        parameters(sys);
        observed = observed(sys),
        constraints = constraints(sys),
        defaults = defaults(sys),
        name = nameof(sys),
        metadata = get_metadata(sys),
        checks = false
    )
end

function calculate_gradient(sys::OptimizationSystem)
    expand_derivatives.(gradient(objective(sys), unknowns(sys)))
end

function generate_gradient(sys::OptimizationSystem, vs = unknowns(sys),
        ps = parameters(sys); kwargs...)
    grad = calculate_gradient(sys)
    p = reorder_parameters(sys, ps)
    return build_function_wrapper(sys, grad, vs, p...; kwargs...)
end

function calculate_hessian(sys::OptimizationSystem)
    expand_derivatives.(hessian(objective(sys), unknowns(sys)))
end

function generate_hessian(
        sys::OptimizationSystem, vs = unknowns(sys), ps = parameters(sys);
        sparse = false, kwargs...)
    if sparse
        hess = sparsehessian(objective(sys), unknowns(sys))
    else
        hess = calculate_hessian(sys)
    end
    p = reorder_parameters(sys, ps)
    return build_function_wrapper(sys, hess, vs, p...; kwargs...)
end

function generate_function(sys::OptimizationSystem, vs = unknowns(sys),
        ps = parameters(sys);
        kwargs...)
    eqs = objective(sys)
    p = reorder_parameters(sys, ps)
    return build_function_wrapper(sys, eqs, vs, p...; kwargs...)
end

function namespace_objective(sys::AbstractSystem)
    op = objective(sys)
    namespace_expr(op, sys)
end

function objective(sys)
    op = get_op(sys)
    systems = get_systems(sys)
    if isempty(systems)
        op
    else
        op + reduce(+, map(sys_ -> namespace_objective(sys_), systems))
    end
end

namespace_constraint(eq::Equation, sys) = namespace_equation(eq, sys)

namespace_constraint(ineq::Inequality, sys) = namespace_inequality(ineq, sys)

function namespace_inequality(ineq::Inequality, sys, n = nameof(sys))
    _lhs = namespace_expr(ineq.lhs, sys, n)
    _rhs = namespace_expr(ineq.rhs, sys, n)
    Inequality(_lhs,
        _rhs,
        ineq.relational_op)
end

function namespace_constraints(sys)
    cstrs = constraints(sys)
    isempty(cstrs) && return Vector{Union{Equation, Inequality}}(undef, 0)
    map(cstr -> namespace_constraint(cstr, sys), cstrs)
end

function constraints(sys)
    cs = get_constraints(sys)
    systems = get_systems(sys)
    isempty(systems) ? cs : [cs; reduce(vcat, namespace_constraints.(systems))]
end

hessian_sparsity(sys::OptimizationSystem) = hessian_sparsity(get_op(sys), unknowns(sys))

"""
```julia
DiffEqBase.OptimizationProblem{iip}(sys::OptimizationSystem, u0map,
                                    parammap = DiffEqBase.NullParameters();
                                    grad = false,
                                    hess = false, sparse = false,
                                    cons_j = false, cons_h = false,
                                    checkbounds = false,
                                    linenumbers = true, parallel = SerialForm(),
                                    kwargs...) where {iip}
```

Generates an OptimizationProblem from an OptimizationSystem and allows for automatically
symbolically calculating numerical enhancements.

Certain solvers require setting `cons_j`, `cons_h` to `true` for constrained-optimization problems.
"""
function DiffEqBase.OptimizationProblem(sys::OptimizationSystem, args...; kwargs...)
    DiffEqBase.OptimizationProblem{true}(sys::OptimizationSystem, args...; kwargs...)
end
function DiffEqBase.OptimizationProblem{iip}(sys::OptimizationSystem, u0map,
        parammap = DiffEqBase.NullParameters();
        lb = nothing, ub = nothing,
        grad = false,
        hess = false, sparse = false,
        cons_j = false, cons_h = false,
        cons_sparse = false, checkbounds = false,
        linenumbers = true, parallel = SerialForm(),
        eval_expression = false, eval_module = @__MODULE__,
        use_union = false,
        checks = true,
        kwargs...) where {iip}
    if !iscomplete(sys)
        error("A completed `OptimizationSystem` is required. Call `complete` or `structural_simplify` on the system before creating a `OptimizationProblem`")
    end
    if haskey(kwargs, :lcons) || haskey(kwargs, :ucons)
        Base.depwarn(
            "`lcons` and `ucons` are deprecated. Specify constraints directly instead.",
            :OptimizationProblem, force = true)
    end

    dvs = unknowns(sys)
    ps = parameters(sys)
    cstr = constraints(sys)

    if isnothing(lb) && isnothing(ub) # use the symbolically specified bounds
        lb = first.(getbounds.(dvs))
        ub = last.(getbounds.(dvs))
        isboolean = symtype.(unwrap.(dvs)) .<: Bool
        lb[isboolean] .= 0
        ub[isboolean] .= 1
    else # use the user supplied variable bounds
        xor(isnothing(lb), isnothing(ub)) &&
            throw(ArgumentError("Expected both `lb` and `ub` to be supplied"))
        !isnothing(lb) && length(lb) != length(dvs) &&
            throw(ArgumentError("Expected both `lb` to be of the same length as the vector of optimization variables"))
        !isnothing(ub) && length(ub) != length(dvs) &&
            throw(ArgumentError("Expected both `ub` to be of the same length as the vector of optimization variables"))
    end

    int = symtype.(unwrap.(dvs)) .<: Integer

    defs = defaults(sys)
    defs = mergedefaults(defs, parammap, ps)
    defs = mergedefaults(defs, u0map, dvs)

    u0 = varmap_to_vars(u0map, dvs; defaults = defs, tofloat = false)
    if parammap isa MTKParameters
        p = parammap
    elseif has_index_cache(sys) && get_index_cache(sys) !== nothing
        p = MTKParameters(sys, parammap, u0map)
    else
        p = varmap_to_vars(parammap, ps; defaults = defs, tofloat = false, use_union)
    end
    lb = varmap_to_vars(dvs .=> lb, dvs; defaults = defs, tofloat = false, use_union)
    ub = varmap_to_vars(dvs .=> ub, dvs; defaults = defs, tofloat = false, use_union)

    if !isnothing(lb) && all(lb .== -Inf) && !isnothing(ub) && all(ub .== Inf)
        lb = nothing
        ub = nothing
    end

    f = let _f = eval_or_rgf(
            generate_function(
                sys, checkbounds = checkbounds, linenumbers = linenumbers,
                expression = Val{true}, wrap_mtkparameters = false);
            eval_expression,
            eval_module)
        __f(u, p) = _f(u, p)
        __f(u, p::MTKParameters) = _f(u, p...)
        __f
    end
    obj_expr = subs_constants(objective(sys))

    if grad
        _grad = let (grad_oop, grad_iip) = eval_or_rgf.(
                generate_gradient(
                    sys, checkbounds = checkbounds,
                    linenumbers = linenumbers,
                    parallel = parallel, expression = Val{true},
                    wrap_mtkparameters = false);
                eval_expression,
                eval_module)
            _grad(u, p) = grad_oop(u, p)
            _grad(J, u, p) = (grad_iip(J, u, p); J)
            _grad(u, p::MTKParameters) = grad_oop(u, p...)
            _grad(J, u, p::MTKParameters) = (grad_iip(J, u, p...); J)
            _grad
        end
    else
        _grad = nothing
    end

    if hess
        _hess = let (hess_oop, hess_iip) = eval_or_rgf.(
                generate_hessian(
                    sys, checkbounds = checkbounds,
                    linenumbers = linenumbers,
                    sparse = sparse, parallel = parallel,
                    expression = Val{true}, wrap_mtkparameters = false);
                eval_expression,
                eval_module)
            _hess(u, p) = hess_oop(u, p)
            _hess(J, u, p) = (hess_iip(J, u, p); J)
            _hess(u, p::MTKParameters) = hess_oop(u, p...)
            _hess(J, u, p::MTKParameters) = (hess_iip(J, u, p...); J)
            _hess
        end
    else
        _hess = nothing
    end

    if sparse
        hess_prototype = hessian_sparsity(sys)
    else
        hess_prototype = nothing
    end

    observedfun = ObservedFunctionCache(sys; eval_expression, eval_module, checkbounds)

    if length(cstr) > 0
        @named cons_sys = ConstraintsSystem(cstr, dvs, ps; checks)
        cons_sys = complete(cons_sys)
        cons, lcons_, ucons_ = generate_function(cons_sys, checkbounds = checkbounds,
            linenumbers = linenumbers,
            expression = Val{true}; wrap_mtkparameters = false)
        cons = let (cons_oop, cons_iip) = eval_or_rgf.(cons; eval_expression, eval_module)
            _cons(u, p) = cons_oop(u, p)
            _cons(resid, u, p) = cons_iip(resid, u, p)
            _cons(u, p::MTKParameters) = cons_oop(u, p...)
            _cons(resid, u, p::MTKParameters) = cons_iip(resid, u, p...)
        end
        if cons_j
            _cons_j = let (cons_jac_oop, cons_jac_iip) = eval_or_rgf.(
                    generate_jacobian(cons_sys;
                        checkbounds = checkbounds,
                        linenumbers = linenumbers,
                        parallel = parallel, expression = Val{true},
                        sparse = cons_sparse, wrap_mtkparameters = false);
                    eval_expression,
                    eval_module)
                _cons_j(u, p) = cons_jac_oop(u, p)
                _cons_j(J, u, p) = (cons_jac_iip(J, u, p); J)
                _cons_j(u, p::MTKParameters) = cons_jac_oop(u, p...)
                _cons_j(J, u, p::MTKParameters) = (cons_jac_iip(J, u, p...); J)
                _cons_j
            end
        else
            _cons_j = nothing
        end
        if cons_h
            _cons_h = let (cons_hess_oop, cons_hess_iip) = eval_or_rgf.(
                    generate_hessian(
                        cons_sys, checkbounds = checkbounds,
                        linenumbers = linenumbers,
                        sparse = cons_sparse, parallel = parallel,
                        expression = Val{true}, wrap_mtkparameters = false);
                    eval_expression,
                    eval_module)
                _cons_h(u, p) = cons_hess_oop(u, p)
                _cons_h(J, u, p) = (cons_hess_iip(J, u, p); J)
                _cons_h(u, p::MTKParameters) = cons_hess_oop(u, p...)
                _cons_h(J, u, p::MTKParameters) = (cons_hess_iip(J, u, p...); J)
                _cons_h
            end
        else
            _cons_h = nothing
        end
        cons_expr = subs_constants(constraints(cons_sys))

        if !haskey(kwargs, :lcons) && !haskey(kwargs, :ucons) # use the symbolically specified bounds
            lcons = lcons_
            ucons = ucons_
        else # use the user supplied constraints bounds
            (haskey(kwargs, :lcons) ⊻ haskey(kwargs, :ucons)) &&
                throw(ArgumentError("Expected both `ucons` and `lcons` to be supplied"))
            haskey(kwargs, :lcons) && length(kwargs[:lcons]) != length(cstr) &&
                throw(ArgumentError("Expected `lcons` to be of the same length as the vector of constraints"))
            haskey(kwargs, :ucons) && length(kwargs[:ucons]) != length(cstr) &&
                throw(ArgumentError("Expected `ucons` to be of the same length as the vector of constraints"))
            lcons = haskey(kwargs, :lcons)
            ucons = haskey(kwargs, :ucons)
        end

        if cons_sparse
            cons_jac_prototype = jacobian_sparsity(cons_sys)
            cons_hess_prototype = hessian_sparsity(cons_sys)
        else
            cons_jac_prototype = nothing
            cons_hess_prototype = nothing
        end
        _f = DiffEqBase.OptimizationFunction{iip}(f,
            sys = sys,
            SciMLBase.NoAD();
            grad = _grad,
            hess = _hess,
            hess_prototype = hess_prototype,
            cons = cons,
            cons_j = _cons_j,
            cons_h = _cons_h,
            cons_jac_prototype = cons_jac_prototype,
            cons_hess_prototype = cons_hess_prototype,
            expr = obj_expr,
            cons_expr = cons_expr,
            observed = observedfun)
        OptimizationProblem{iip}(_f, u0, p; lb = lb, ub = ub, int = int,
            lcons = lcons, ucons = ucons, kwargs...)
    else
        _f = DiffEqBase.OptimizationFunction{iip}(f,
            sys = sys,
            SciMLBase.NoAD();
            grad = _grad,
            hess = _hess,
            hess_prototype = hess_prototype,
            expr = obj_expr,
            observed = observedfun)
        OptimizationProblem{iip}(_f, u0, p; lb = lb, ub = ub, int = int,
            kwargs...)
    end
end

"""
```julia
DiffEqBase.OptimizationProblemExpr{iip}(sys::OptimizationSystem,
                                        parammap = DiffEqBase.NullParameters();
                                        u0 = nothing,
                                        grad = false,
                                        hes = false, sparse = false,
                                        checkbounds = false,
                                        linenumbers = true, parallel = SerialForm(),
                                        kwargs...) where {iip}
```

Generates a Julia expression for an OptimizationProblem from an
OptimizationSystem and allows for automatically symbolically
calculating numerical enhancements.
"""
struct OptimizationProblemExpr{iip} end

function OptimizationProblemExpr(sys::OptimizationSystem, args...; kwargs...)
    OptimizationProblemExpr{true}(sys::OptimizationSystem, args...; kwargs...)
end

function OptimizationProblemExpr{iip}(sys::OptimizationSystem, u0map,
        parammap = DiffEqBase.NullParameters();
        lb = nothing, ub = nothing,
        grad = false,
        hess = false, sparse = false,
        cons_j = false, cons_h = false,
        checkbounds = false,
        linenumbers = false, parallel = SerialForm(),
        eval_expression = false, eval_module = @__MODULE__,
        use_union = false,
        kwargs...) where {iip}
    if !iscomplete(sys)
        error("A completed `OptimizationSystem` is required. Call `complete` or `structural_simplify` on the system before creating a `OptimizationProblemExpr`")
    end
    if haskey(kwargs, :lcons) || haskey(kwargs, :ucons)
        Base.depwarn(
            "`lcons` and `ucons` are deprecated. Specify constraints directly instead.",
            :OptimizationProblem, force = true)
    end

    dvs = unknowns(sys)
    ps = parameters(sys)
    cstr = constraints(sys)

    if isnothing(lb) && isnothing(ub) # use the symbolically specified bounds
        lb = first.(getbounds.(dvs))
        ub = last.(getbounds.(dvs))
        isboolean = symtype.(unwrap.(dvs)) .<: Bool
        lb[isboolean] .= 0
        ub[isboolean] .= 1
    else # use the user supplied variable bounds
        xor(isnothing(lb), isnothing(ub)) &&
            throw(ArgumentError("Expected both `lb` and `ub` to be supplied"))
        !isnothing(lb) && length(lb) != length(dvs) &&
            throw(ArgumentError("Expected `lb` to be of the same length as the vector of optimization variables"))
        !isnothing(ub) && length(ub) != length(dvs) &&
            throw(ArgumentError("Expected `ub` to be of the same length as the vector of optimization variables"))
    end

    int = symtype.(unwrap.(dvs)) .<: Integer

    defs = defaults(sys)
    defs = mergedefaults(defs, parammap, ps)
    defs = mergedefaults(defs, u0map, dvs)

    u0 = varmap_to_vars(u0map, dvs; defaults = defs, tofloat = false)
    if has_index_cache(sys) && get_index_cache(sys) !== nothing
        p = MTKParameters(sys, parammap, u0map)
    else
        p = varmap_to_vars(parammap, ps; defaults = defs, tofloat = false, use_union)
    end
    lb = varmap_to_vars(dvs .=> lb, dvs; defaults = defs, tofloat = false, use_union)
    ub = varmap_to_vars(dvs .=> ub, dvs; defaults = defs, tofloat = false, use_union)

    if !isnothing(lb) && all(lb .== -Inf) && !isnothing(ub) && all(ub .== Inf)
        lb = nothing
        ub = nothing
    end

    idx = iip ? 2 : 1
    f = generate_function(sys, checkbounds = checkbounds, linenumbers = linenumbers,
        expression = Val{true})
    if grad
        _grad = eval_or_rgf(
            generate_gradient(
                sys, checkbounds = checkbounds, linenumbers = linenumbers,
                parallel = parallel, expression = Val{true})[idx];
            eval_expression,
            eval_module)
    else
        _grad = :nothing
    end

    if hess
        _hess = eval_or_rgf(
            generate_hessian(sys, checkbounds = checkbounds, linenumbers = linenumbers,
                sparse = sparse, parallel = parallel,
                expression = Val{false})[idx];
            eval_expression,
            eval_module)
    else
        _hess = :nothing
    end

    if sparse
        hess_prototype = hessian_sparsity(sys)
    else
        hess_prototype = nothing
    end

    obj_expr = toexpr(subs_constants(objective(sys)))
    pairs_arr = if p isa SciMLBase.NullParameters
        [Symbol(_s) => Expr(:ref, :x, i) for (i, _s) in enumerate(dvs)]
    else
        vcat([Symbol(_s) => Expr(:ref, :x, i) for (i, _s) in enumerate(dvs)],
            [Symbol(_p) => p[i] for (i, _p) in enumerate(ps)])
    end
    rep_pars_vals!(obj_expr, pairs_arr)

    if length(cstr) > 0
        @named cons_sys = ConstraintsSystem(cstr, dvs, ps)
        cons, lcons_, ucons_ = generate_function(cons_sys, checkbounds = checkbounds,
            linenumbers = linenumbers,
            expression = Val{true})
        cons = eval_or_rgf(cons; eval_expression, eval_module)
        if cons_j
            _cons_j = eval_or_rgf(
                generate_jacobian(cons_sys; expression = Val{true}, sparse = sparse)[2];
                eval_expression, eval_module)
        else
            _cons_j = nothing
        end
        if cons_h
            _cons_h = eval_or_rgf(
                generate_hessian(cons_sys; expression = Val{true}, sparse = sparse)[2];
                eval_expression, eval_module)
        else
            _cons_h = nothing
        end

        cons_expr = toexpr.(subs_constants(constraints(cons_sys)))
        rep_pars_vals!.(cons_expr, Ref(pairs_arr))

        if !haskey(kwargs, :lcons) && !haskey(kwargs, :ucons) # use the symbolically specified bounds
            lcons = lcons_
            ucons = ucons_
        else # use the user supplied constraints bounds
            (haskey(kwargs, :lcons) ⊻ haskey(kwargs, :ucons)) &&
                throw(ArgumentError("Expected both `ucons` and `lcons` to be supplied"))
            haskey(kwargs, :lcons) && length(kwargs[:lcons]) != length(cstr) &&
                throw(ArgumentError("Expected `lcons` to be of the same length as the vector of constraints"))
            haskey(kwargs, :ucons) && length(kwargs[:ucons]) != length(cstr) &&
                throw(ArgumentError("Expected `ucons` to be of the same length as the vector of constraints"))
            lcons = haskey(kwargs, :lcons)
            ucons = haskey(kwargs, :ucons)
        end

        if sparse
            cons_jac_prototype = jacobian_sparsity(cons_sys)
            cons_hess_prototype = hessian_sparsity(cons_sys)
        else
            cons_jac_prototype = nothing
            cons_hess_prototype = nothing
        end

        quote
            f = $f
            p = $p
            u0 = $u0
            grad = $_grad
            hess = $_hess
            lb = $lb
            ub = $ub
            int = $int
            cons = $cons[1]
            lcons = $lcons
            ucons = $ucons
            cons_j = $_cons_j
            cons_h = $_cons_h
            _f = OptimizationFunction{iip}(f, SciMLBase.NoAD();
                grad = grad,
                hess = hess,
                hess_prototype = hess_prototype,
                cons = cons,
                cons_j = cons_j,
                cons_h = cons_h,
                cons_jac_prototype = cons_jac_prototype,
                cons_hess_prototype = cons_hess_prototype,
                expr = obj_expr,
                cons_expr = cons_expr)
            OptimizationProblem{$iip}(
                _f, u0, p; lb = lb, ub = ub, int = int, lcons = lcons,
                ucons = ucons, kwargs...)
        end
    else
        quote
            f = $f
            p = $p
            u0 = $u0
            grad = $_grad
            hess = $_hess
            lb = $lb
            ub = $ub
            int = $int
            _f = OptimizationFunction{iip}(f, SciMLBase.NoAD();
                grad = grad,
                hess = hess,
                hess_prototype = hess_prototype,
                expr = obj_expr)
            OptimizationProblem{$iip}(_f, u0, p; lb = lb, ub = ub, int = int, kwargs...)
        end
    end
end

function structural_simplify(sys::OptimizationSystem; split = true, kwargs...)
    sys = flatten(sys)
    cons = constraints(sys)
    econs = Equation[]
    icons = similar(cons, 0)
    for e in cons
        if e isa Equation
            push!(econs, e)
        else
            push!(icons, e)
        end
    end
    nlsys = NonlinearSystem(econs, unknowns(sys), parameters(sys); name = :___tmp_nlsystem)
    snlsys = structural_simplify(nlsys; fully_determined = false, kwargs...)
    obs = observed(snlsys)
    subs = Dict(eq.lhs => eq.rhs for eq in observed(snlsys))
    seqs = equations(snlsys)
    cons_simplified = similar(cons, length(icons) + length(seqs))
    for (i, eq) in enumerate(Iterators.flatten((seqs, icons)))
        cons_simplified[i] = fixpoint_sub(eq, subs)
    end
    newsts = setdiff(unknowns(sys), keys(subs))
    @set! sys.constraints = cons_simplified
    @set! sys.observed = [observed(sys); obs]
    neweqs = fixpoint_sub.(equations(sys), (subs,))
    @set! sys.op = length(neweqs) == 1 ? first(neweqs) : neweqs
    @set! sys.unknowns = newsts
    sys = complete(sys; split)
    return sys
end
