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
    states::Vector
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
    The hierarchical parent system before simplification.
    """
    parent::Any

    function OptimizationSystem(tag, op, states, ps, var_to_name, observed,
            constraints, name, systems, defaults, metadata = nothing,
            gui_metadata = nothing, complete = false, parent = nothing;
            checks::Union{Bool, Int} = true)
        if checks == true || (checks & CheckUnits) > 0
            unwrap(op) isa Symbolic && check_units(op)
            check_units(observed)
            all_dimensionless([states; ps]) || check_units(constraints)
        end
        new(tag, op, states, ps, var_to_name, observed,
            constraints, name, systems, defaults, metadata, gui_metadata, complete,
            parent)
    end
end

equations(sys::AbstractOptimizationSystem) = objective(sys) # needed for Base.show

function OptimizationSystem(op, states, ps;
        observed = [],
        constraints = [],
        default_u0 = Dict(),
        default_p = Dict(),
        defaults = _merge(Dict(default_u0), Dict(default_p)),
        name = nothing,
        systems = OptimizationSystem[],
        checks = true,
        metadata = nothing,
        gui_metadata = nothing)
    name === nothing &&
        throw(ArgumentError("The `name` keyword must be provided. Please consider using the `@named` macro"))
    constraints = value.(scalarize(constraints))
    states′ = value.(scalarize(states))
    ps′ = value.(scalarize(ps))
    op′ = value(scalarize(op))

    if !(isempty(default_u0) && isempty(default_p))
        Base.depwarn("`default_u0` and `default_p` are deprecated. Use `defaults` instead.",
            :OptimizationSystem, force = true)
    end
    sysnames = nameof.(systems)
    if length(unique(sysnames)) != length(sysnames)
        throw(ArgumentError("System names must be unique."))
    end
    defaults = todict(defaults)
    defaults = Dict(value(k) => value(v) for (k, v) in pairs(defaults))

    var_to_name = Dict()
    process_variables!(var_to_name, defaults, states′)
    process_variables!(var_to_name, defaults, ps′)
    isempty(observed) || collect_var_to_name!(var_to_name, (eq.lhs for eq in observed))

    OptimizationSystem(Threads.atomic_add!(SYSTEM_COUNT, UInt(1)),
        op′, states′, ps′, var_to_name,
        observed,
        constraints,
        name, systems, defaults, metadata, gui_metadata;
        checks = checks)
end

function calculate_gradient(sys::OptimizationSystem)
    expand_derivatives.(gradient(objective(sys), states(sys)))
end

function generate_gradient(sys::OptimizationSystem, vs = states(sys), ps = parameters(sys);
        kwargs...)
    grad = calculate_gradient(sys)
    pre = get_preprocess_constants(grad)
    return build_function(grad, vs, ps; postprocess_fbody = pre,
        kwargs...)
end

function calculate_hessian(sys::OptimizationSystem)
    expand_derivatives.(hessian(objective(sys), states(sys)))
end

function generate_hessian(sys::OptimizationSystem, vs = states(sys), ps = parameters(sys);
        sparse = false, kwargs...)
    if sparse
        hess = sparsehessian(objective(sys), states(sys))
    else
        hess = calculate_hessian(sys)
    end
    pre = get_preprocess_constants(hess)
    return build_function(hess, vs, ps; postprocess_fbody = pre,
        kwargs...)
end

function generate_function(sys::OptimizationSystem, vs = states(sys), ps = parameters(sys);
        kwargs...)
    eqs = subs_constants(objective(sys))
    return build_function(eqs, vs, ps;
        kwargs...)
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

hessian_sparsity(sys::OptimizationSystem) = hessian_sparsity(get_op(sys), states(sys))

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
        use_union = false,
        kwargs...) where {iip}
    if haskey(kwargs, :lcons) || haskey(kwargs, :ucons)
        Base.depwarn("`lcons` and `ucons` are deprecated. Specify constraints directly instead.",
            :OptimizationProblem, force = true)
    end

    dvs = states(sys)
    ps = parameters(sys)
    cstr = constraints(sys)

    if isnothing(lb) && isnothing(ub) # use the symbolically specified bounds
        lb = first.(getbounds.(dvs))
        ub = last.(getbounds.(dvs))
        lb[isbinaryvar.(dvs)] .= 0
        ub[isbinaryvar.(dvs)] .= 1
    else # use the user supplied variable bounds
        xor(isnothing(lb), isnothing(ub)) &&
            throw(ArgumentError("Expected both `lb` and `ub` to be supplied"))
        !isnothing(lb) && length(lb) != length(dvs) &&
            throw(ArgumentError("Expected both `lb` to be of the same length as the vector of optimization variables"))
        !isnothing(ub) && length(ub) != length(dvs) &&
            throw(ArgumentError("Expected both `ub` to be of the same length as the vector of optimization variables"))
    end

    int = isintegervar.(dvs) .| isbinaryvar.(dvs)

    defs = defaults(sys)
    defs = mergedefaults(defs, parammap, ps)
    defs = mergedefaults(defs, u0map, dvs)

    u0 = varmap_to_vars(u0map, dvs; defaults = defs, tofloat = false)
    p = varmap_to_vars(parammap, ps; defaults = defs, tofloat = false, use_union)
    lb = varmap_to_vars(dvs .=> lb, dvs; defaults = defs, tofloat = false, use_union)
    ub = varmap_to_vars(dvs .=> ub, dvs; defaults = defs, tofloat = false, use_union)

    if !isnothing(lb) && all(lb .== -Inf) && !isnothing(ub) && all(ub .== Inf)
        lb = nothing
        ub = nothing
    end

    f = generate_function(sys, checkbounds = checkbounds, linenumbers = linenumbers,
        expression = Val{false})

    obj_expr = subs_constants(objective(sys))

    if grad
        grad_oop, grad_iip = generate_gradient(sys, checkbounds = checkbounds,
            linenumbers = linenumbers,
            parallel = parallel, expression = Val{false})
        _grad(u, p) = grad_oop(u, p)
        _grad(J, u, p) = (grad_iip(J, u, p); J)
    else
        _grad = nothing
    end

    if hess
        hess_oop, hess_iip = generate_hessian(sys, checkbounds = checkbounds,
            linenumbers = linenumbers,
            sparse = sparse, parallel = parallel,
            expression = Val{false})
        _hess(u, p) = hess_oop(u, p)
        _hess(J, u, p) = (hess_iip(J, u, p); J)
    else
        _hess = nothing
    end

    if sparse
        hess_prototype = hessian_sparsity(sys)
    else
        hess_prototype = nothing
    end

    observedfun = let sys = sys, dict = Dict()
        function generated_observed(obsvar, args...)
            obs = get!(dict, value(obsvar)) do
                build_explicit_observed_function(sys, obsvar)
            end
            if args === ()
                let obs = obs
                    (u, p) -> obs(u, p)
                end
            else
                obs(args...)
            end
        end
    end

    if length(cstr) > 0
        @named cons_sys = ConstraintsSystem(cstr, dvs, ps)
        cons, lcons_, ucons_ = generate_function(cons_sys, checkbounds = checkbounds,
            linenumbers = linenumbers,
            expression = Val{false})
        if cons_j
            _cons_j = generate_jacobian(cons_sys; expression = Val{false},
                sparse = cons_sparse)[2]
        else
            _cons_j = nothing
        end
        if cons_h
            _cons_h = generate_hessian(cons_sys; expression = Val{false},
                sparse = cons_sparse)[2]
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
            syms = Symbol.(states(sys)),
            paramsyms = Symbol.(parameters(sys)),
            cons = cons[2],
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
            syms = Symbol.(states(sys)),
            paramsyms = Symbol.(parameters(sys)),
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
        use_union = false,
        kwargs...) where {iip}
    if haskey(kwargs, :lcons) || haskey(kwargs, :ucons)
        Base.depwarn("`lcons` and `ucons` are deprecated. Specify constraints directly instead.",
            :OptimizationProblem, force = true)
    end

    dvs = states(sys)
    ps = parameters(sys)
    cstr = constraints(sys)

    if isnothing(lb) && isnothing(ub) # use the symbolically specified bounds
        lb = first.(getbounds.(dvs))
        ub = last.(getbounds.(dvs))
        lb[isbinaryvar.(dvs)] .= 0
        ub[isbinaryvar.(dvs)] .= 1
    else # use the user supplied variable bounds
        xor(isnothing(lb), isnothing(ub)) &&
            throw(ArgumentError("Expected both `lb` and `ub` to be supplied"))
        !isnothing(lb) && length(lb) != length(dvs) &&
            throw(ArgumentError("Expected `lb` to be of the same length as the vector of optimization variables"))
        !isnothing(ub) && length(ub) != length(dvs) &&
            throw(ArgumentError("Expected `ub` to be of the same length as the vector of optimization variables"))
    end

    int = isintegervar.(dvs) .| isbinaryvar.(dvs)

    defs = defaults(sys)
    defs = mergedefaults(defs, parammap, ps)
    defs = mergedefaults(defs, u0map, dvs)

    u0 = varmap_to_vars(u0map, dvs; defaults = defs, tofloat = false)
    p = varmap_to_vars(parammap, ps; defaults = defs, tofloat = false, use_union)
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
        _grad = generate_gradient(sys, checkbounds = checkbounds, linenumbers = linenumbers,
            parallel = parallel, expression = Val{false})[idx]
    else
        _grad = :nothing
    end

    if hess
        _hess = generate_hessian(sys, checkbounds = checkbounds, linenumbers = linenumbers,
            sparse = sparse, parallel = parallel,
            expression = Val{false})[idx]
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
            expression = Val{false})
        if cons_j
            _cons_j = generate_jacobian(cons_sys; expression = Val{false}, sparse = sparse)[2]
        else
            _cons_j = nothing
        end
        if cons_h
            _cons_h = generate_hessian(cons_sys; expression = Val{false}, sparse = sparse)[2]
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
            syms = $(Symbol.(states(sys)))
            paramsyms = $(Symbol.(parameters(sys)))
            _f = OptimizationFunction{iip}(f, SciMLBase.NoAD();
                grad = grad,
                hess = hess,
                syms = syms,
                paramsyms = paramsyms,
                hess_prototype = hess_prototype,
                cons = cons,
                cons_j = cons_j,
                cons_h = cons_h,
                cons_jac_prototype = cons_jac_prototype,
                cons_hess_prototype = cons_hess_prototype,
                expr = obj_expr,
                cons_expr = cons_expr)
            OptimizationProblem{$iip}(_f, u0, p; lb = lb, ub = ub, int = int, lcons = lcons,
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
            syms = $(Symbol.(states(sys)))
            paramsyms = $(Symbol.(parameters(sys)))
            _f = OptimizationFunction{iip}(f, SciMLBase.NoAD();
                grad = grad,
                hess = hess,
                syms = syms,
                paramsyms = paramsyms,
                hess_prototype = hess_prototype,
                expr = obj_expr)
            OptimizationProblem{$iip}(_f, u0, p; lb = lb, ub = ub, int = int, kwargs...)
        end
    end
end

function structural_simplify(sys::OptimizationSystem; kwargs...)
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
    nlsys = NonlinearSystem(econs, states(sys), parameters(sys); name = :___tmp_nlsystem)
    snlsys = structural_simplify(nlsys; fully_determined = false, kwargs...)
    obs = observed(snlsys)
    subs = Dict(eq.lhs => eq.rhs for eq in observed(snlsys))
    seqs = equations(snlsys)
    cons_simplified = similar(cons, length(icons) + length(seqs))
    for (i, eq) in enumerate(Iterators.flatten((seqs, icons)))
        cons_simplified[i] = fixpoint_sub(eq, subs)
    end
    newsts = setdiff(states(sys), keys(subs))
    @set! sys.constraints = cons_simplified
    @set! sys.observed = [observed(sys); obs]
    neweqs = fixpoint_sub.(equations(sys), (subs,))
    @set! sys.op = length(neweqs) == 1 ? first(neweqs) : neweqs
    @set! sys.states = newsts
    return sys
end
