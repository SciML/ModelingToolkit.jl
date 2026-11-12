"""$(function_docstring(OptimizationFunction, false, [:jac, :grad, :hess, :cons_h, :cons_j, :adtype]; extra_kwargs = WEIGHTS_KWARGS))"""
function SciMLBase.OptimizationFunction(sys::System, args...; kwargs...)
    return OptimizationFunction{true}(sys, args...; kwargs...)
end

function SciMLBase.OptimizationFunction{iip}(
        sys::System, adtype::ADTypes.AbstractADType; kwargs...
    ) where {iip}
    return OptimizationFunction{iip}(sys; adtype, kwargs...)
end

function SciMLBase.OptimizationFunction{iip}(
        sys::System;
        u0 = nothing, p = nothing, t = nothing, grad = false, hess = false,
        sparse = false, cons_j = false, cons_h = false,
        cons_sparse = false, adtype::ADTypes.AbstractADType = SciMLBase.NoAD(),
        linenumbers = true, eval_expression = false,
        eval_module = @__MODULE__,
        simplify = false, check_compatibility = true, checkbounds = false,
        expression = Val{false}, optimize = nothing, weights = nothing,
        compiler_options::CompilerOptions = CompilerOptions(), kwargs...
    ) where {iip}
    opts = SciMLFunctionOptions(;
        u0, p, t, sparse, simplify, expression, check_compatibility,
        eval_expression, eval_module, compiler_options, checkbounds, optimize, kwargs...,
    )
    return OptimizationFunction{iip}(
        sys, opts; grad, hess, cons_j, cons_h, cons_sparse, weights, adtype
    )
end

"""
    SciMLBase.OptimizationFunction{iip}(sys::System, opts::SciMLFunctionOptions; kwargs...)

Public entry point that builds an `OptimizationFunction` directly from a pre-assembled
`SciMLFunctionOptions`, bypassing the `kwargs...` wrapper above.
"""
function SciMLBase.OptimizationFunction{iip}(
        sys::System, opts::SciMLFunctionOptions{E};
        grad::Bool = false, hess::Bool = false, cons_j::Bool = false, cons_h::Bool = false,
        cons_sparse::Bool = false, weights = nothing,
        adtype::ADTypes.AbstractADType = SciMLBase.NoAD()
    ) where {iip, E}
    check_complete(sys, OptimizationFunction)
    opts.check_compatibility && check_compatible_system(OptimizationFunction, sys)

    if weights !== nothing
        sys = system_with_cost_weights(sys, weights)
    end

    (; u0, p, sparse, simplify) = opts
    codegen_opts = opts.codegen

    f = generate_cost(sys, codegen_opts)

    if grad
        _grad = generate_cost_gradient(sys, codegen_opts)
    else
        _grad = nothing
    end
    if hess
        _hess,
            hess_prototype = generate_cost_hessian(
            sys, codegen_opts; sparse, simplify, return_sparsity = true
        )
    else
        _hess = hess_prototype = nothing
        if sparse
            hess_prototype = cost_hessian_sparsity(sys)
        end
    end
    constraint_fields = generate_constraint_fields(
        sys, codegen_opts; cons_j, cons_h, cons_sparse, simplify
    )

    obj_expr = Code.toexpr(expand(cost(sys)))

    observedfun = ObservedFunctionCache(sys, codegen_opts)

    args = (; f, ad = adtype)
    kwargs = (;
        sys = sys,
        grad = _grad,
        hess = _hess,
        hess_prototype = hess_prototype,
        constraint_fields...,
        expr = obj_expr,
        observed = observedfun,
    )

    return maybe_codegen_scimlfn(Val{E}, OptimizationFunction{iip}, args; kwargs...)
end

"""$(function_docstring(MultiObjectiveOptimizationFunction, false, [:jac, :hess, :cons_h, :cons_j, :adtype]; extra_body = "The generated objective is vector-valued: it returns [`costs`](@ref) elementwise rather than the `consolidate`d scalar of `OptimizationFunction`."))"""
function SciMLBase.MultiObjectiveOptimizationFunction(sys::System, args...; kwargs...)
    return MultiObjectiveOptimizationFunction{true}(sys, args...; kwargs...)
end

function SciMLBase.MultiObjectiveOptimizationFunction{iip}(
        sys::System, adtype::ADTypes.AbstractADType; kwargs...
    ) where {iip}
    return MultiObjectiveOptimizationFunction{iip}(sys; adtype, kwargs...)
end

function SciMLBase.MultiObjectiveOptimizationFunction{iip}(
        sys::System;
        u0 = nothing, p = nothing, t = nothing, jac = false, hess = false,
        sparse = false, cons_j = false, cons_h = false,
        cons_sparse = false, adtype::ADTypes.AbstractADType = SciMLBase.NoAD(),
        linenumbers = true, eval_expression = false,
        eval_module = @__MODULE__,
        simplify = false, check_compatibility = true, checkbounds = false,
        expression = Val{false}, optimize = nothing, weights = nothing,
        compiler_options::CompilerOptions = CompilerOptions(), kwargs...
    ) where {iip}
    opts = SciMLFunctionOptions(;
        u0, p, t, jac, sparse, simplify, expression, check_compatibility,
        eval_expression, eval_module, compiler_options, checkbounds, optimize, kwargs...,
    )
    return MultiObjectiveOptimizationFunction{iip}(
        sys, opts; hess, cons_j, cons_h, cons_sparse, weights, adtype
    )
end

"""
    SciMLBase.MultiObjectiveOptimizationFunction{iip}(sys::System, opts::SciMLFunctionOptions; kwargs...)

Public entry point that builds a `MultiObjectiveOptimizationFunction` directly from a
pre-assembled `SciMLFunctionOptions`, bypassing the `kwargs...` wrapper above.
"""
function SciMLBase.MultiObjectiveOptimizationFunction{iip}(
        sys::System, opts::SciMLFunctionOptions{E};
        hess::Bool = false, cons_j::Bool = false, cons_h::Bool = false,
        cons_sparse::Bool = false, weights = nothing,
        adtype::ADTypes.AbstractADType = SciMLBase.NoAD()
    ) where {iip, E}
    check_complete(sys, MultiObjectiveOptimizationFunction)
    opts.check_compatibility &&
        check_compatible_system(MultiObjectiveOptimizationFunction, sys)
    weights !== nothing && throw(
        ArgumentError(
            "`weights` scalarizes the costs of `sys`; a multiobjective objective keeps \
            them separate. Use `OptimizationFunction` for a weighted-sum objective."
        )
    )

    (; u0, p, jac, sparse, simplify) = opts
    codegen_opts = opts.codegen

    f = generate_multiobjective_cost(sys, codegen_opts)

    if jac
        _jac = generate_multiobjective_jacobian(sys, codegen_opts)
    else
        _jac = nothing
    end
    if hess
        _hess,
            hess_prototype = generate_multiobjective_hessian(
            sys, codegen_opts; sparse, simplify, return_sparsity = true
        )
    else
        _hess = hess_prototype = nothing
        if sparse
            _, hess_prototype = calculate_multiobjective_hessian(
                sys; sparse = true, return_sparsity = true
            )
        end
    end
    constraint_fields = generate_constraint_fields(
        sys, codegen_opts; cons_j, cons_h, cons_sparse, simplify
    )

    obj_expr = Code.toexpr.(expand.(costs(sys)))

    observedfun = ObservedFunctionCache(sys, codegen_opts)

    args = (; f, ad = adtype)
    kwargs = (;
        sys = sys,
        jac = _jac,
        hess = _hess,
        hess_prototype = hess_prototype,
        constraint_fields...,
        expr = obj_expr,
        observed = observedfun,
    )

    return maybe_codegen_scimlfn(Val{E}, MultiObjectiveOptimizationFunction{iip}, args; kwargs...)
end

"""
    LazyConstraintExprs(sys::System, len::Int)

The `cons_expr` of an `OptimizationFunction` built from `sys`: an `AbstractVector{Expr}`
with one expression per row of the constraint function, the `Code.toexpr` of the expanded
rows of [`canonical_constraints`](@ref). Scalarizing and expanding large array-valued
constraints is expensive and only expression-graph consumers need it, so the expressions
are built on the first access to an element and cached; solvers that only call the
generated functions never pay for them.
"""
mutable struct LazyConstraintExprs <: AbstractVector{Expr}
    const sys::System
    const len::Int
    exprs::Union{Nothing, Vector{Expr}}
    const lock::ReentrantLock
end

function LazyConstraintExprs(sys::System, len::Int)
    return LazyConstraintExprs(sys, len, nothing, ReentrantLock())
end

Base.size(c::LazyConstraintExprs) = (c.len,)
Base.IndexStyle(::Type{LazyConstraintExprs}) = IndexLinear()
Base.getindex(c::LazyConstraintExprs, i::Int) = materialize_constraint_exprs(c)[i]

function materialize_constraint_exprs(c::LazyConstraintExprs)
    return @lock c.lock begin
        exprs = c.exprs
        if exprs === nothing
            exprs = Expr[Code.toexpr(expand(row)) for row in canonical_constraints(c.sys)]
            c.exprs = exprs
        end
        exprs
    end::Vector{Expr}
end

"""
    $(TYPEDSIGNATURES)

The constraint-related fields of an `OptimizationFunction` or
`MultiObjectiveOptimizationFunction` built from `sys`: the constraint function `cons`,
its symbolic form `cons_expr` (a `LazyConstraintExprs`), and, when requested, the
constraint jacobian `cons_j` and hessian `cons_h` with their sparsity prototypes. Every
field is `nothing` when `sys` has no constraints.
"""
function generate_constraint_fields(
        sys::System, codegen_opts::GeneratedFunctionOptions;
        cons_j::Bool, cons_h::Bool, cons_sparse::Bool, simplify::Bool
    )
    cstr = constraints(sys)
    if isempty(cstr)
        return (;
            cons = nothing, cons_j = nothing, cons_jac_prototype = nothing,
            cons_h = nothing, cons_hess_prototype = nothing, cons_expr = nothing,
        )
    end
    cons = generate_cons(sys, codegen_opts)
    _cons_j = cons_jac_prototype = nothing
    if cons_j
        _cons_j,
            cons_jac_prototype = generate_constraint_jacobian(
            sys, codegen_opts; simplify, sparse = cons_sparse, return_sparsity = true
        )
    end
    _cons_h = cons_hess_prototype = nothing
    if cons_h
        _cons_h,
            cons_hess_prototype = generate_constraint_hessian(
            sys, codegen_opts; simplify, sparse = cons_sparse, return_sparsity = true
        )
    end
    cons_expr = LazyConstraintExprs(sys, sum(constraint_length, cstr))
    return (;
        cons, cons_j = _cons_j, cons_jac_prototype, cons_h = _cons_h,
        cons_hess_prototype, cons_expr,
    )
end

"""$(problem_docstring(SciMLBase.OptimizationProblem, OptimizationFunction, false; init = false, extra_kwargs = WEIGHTS_KWARGS * MULTIOBJECTIVE_KWARGS * ADTYPE_PROBLEM_KWARGS))"""
function SciMLBase.OptimizationProblem(sys::System, args...; kwargs...)
    return OptimizationProblem{true}(sys, args...; kwargs...)
end

function SciMLBase.OptimizationProblem{iip}(
        sys::System, op; lb = nothing,
        ub = nothing, check_compatibility = true, expression = Val{false},
        multiobjective::Bool = false,
        kwargs...
    ) where {iip}
    check_complete(sys, OptimizationProblem)
    check_compatibility && check_compatible_system(OptimizationProblem, sys)

    f, u0,
        p = process_SciMLProblem(
        (multiobjective ? MultiObjectiveOptimizationFunction{iip} : OptimizationFunction{iip}),
        sys, op;
        check_compatibility, tofloat = false, check_length = false, expression, kwargs...
    )

    dvs = flat_unknowns(sys)
    int = symtype.(unwrap.(dvs)) .<: Integer
    if lb === nothing && ub === nothing
        lb = first.(getbounds.(dvs))
        ub = last.(getbounds.(dvs))
        isboolean = symtype.(unwrap.(dvs)) .<: Bool
        lb[isboolean] .= 0
        ub[isboolean] .= 1
    else
        xor(isnothing(lb), isnothing(ub)) &&
            throw(ArgumentError("Expected both `lb` and `ub` to be supplied"))
        !isnothing(lb) && length(lb) != length(dvs) &&
            throw(ArgumentError("Expected both `lb` to be of the same length as the vector of optimization variables"))
        !isnothing(ub) && length(ub) != length(dvs) &&
            throw(ArgumentError("Expected both `ub` to be of the same length as the vector of optimization variables"))
    end

    op = build_operating_point(sys, op)
    lbmap = as_atomic_dict_with_defaults(Dict{SymbolicT, SymbolicT}(dvs .=> lb), COMMON_NOTHING)
    left_merge!(lbmap, op)
    lb = varmap_to_vars(lbmap, dvs; tofloat = false)
    ubmap = as_atomic_dict_with_defaults(Dict{SymbolicT, SymbolicT}(dvs .=> ub), COMMON_NOTHING)
    left_merge!(ubmap, op)
    ub = varmap_to_vars(ubmap, dvs; tofloat = false)

    if !isnothing(lb) && all(lb .== -Inf) && !isnothing(ub) && all(ub .== Inf)
        lb = nothing
        ub = nothing
    end

    cstr = constraints(sys)
    if isempty(cstr)
        lcons = ucons = nothing
    else
        lcons = Float64[]
        for c in cstr
            append!(lcons, Iterators.repeated(c isa Equation ? 0.0 : -Inf, constraint_length(c)))
        end
        ucons = zeros(length(lcons))
    end

    kwargs = process_kwargs(sys; kwargs...)
    ptype = getmetadata(sys, ProblemTypeCtx, nothing)
    kwargs = (; lb, ub, int, lcons, ucons, problem_type = ptype, kwargs...)
    args = (; f, u0, p)
    return maybe_codegen_scimlproblem(expression, OptimizationProblem{iip}, args; kwargs...)
end

function check_compatible_system(
        T::Union{
            Type{OptimizationFunction}, Type{MultiObjectiveOptimizationFunction},
            Type{OptimizationProblem},
        }, sys::System
    )
    check_time_independent(sys, T)
    check_not_dde(sys)
    check_has_cost(sys, T)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
    return check_no_equations(sys, T)
end

"""
    $(TYPEDSIGNATURES)

Return a `consolidate` function computing `sum(weights .* costs)` plus the sum of the
consolidated `subcosts` of all subsystems. See [`system_with_cost_weights`](@ref).
"""
function weighted_consolidate(weights)
    ws = unwrap.(weights)
    return function (costs, subcosts)
        return _sum_costs(SymbolicT[w * c for (w, c) in zip(ws, costs)]) +
            _sum_costs(subcosts)
    end
end

"""
    $(TYPEDSIGNATURES)

Return a copy of `sys` whose `consolidate` function computes the `weights`-weighted sum of
the top-level costs of `sys` instead of deferring to the system's own `consolidate`. The
costs of subsystems are still consolidated recursively by their own `consolidate` and
added to the result.

`weights` must have one entry per top-level cost of `sys` (that is,
`length(weights) == length(get_costs(sys))`). Entries may be real numbers or symbolic
parameters of `sys`. Symbolic weights must be declared as parameters of `sys` so that
they are discoverable in the parameter object (`prob.p`) and can be updated via `remake`
between solves.
"""
function system_with_cost_weights(sys::System, weights)
    weights isa AbstractVector || throw(
        ArgumentError(
            "`weights` must be a vector with one entry per cost of `sys`."
        )
    )
    cs = get_costs(sys)
    length(weights) == length(cs) || throw(
        ArgumentError(
            "Expected `weights` to have one entry per cost of `sys`, but got \
            $(length(weights)) weights for $(length(cs)) costs."
        )
    )
    for w in weights
        # `Num <: Number`, so `unwrap` before checking for numeric weights.
        w = unwrap(w)
        w isa Number && continue
        symbolic_type(w) === NotSymbolic() && throw(
            ArgumentError(
                "Entries of `weights` must be real numbers or symbolic parameters of \
                `sys`; got `$w`."
            )
        )
        is_parameter(sys, w) || throw(
            ArgumentError(
                "Symbolic weight `$w` is not a parameter of `sys`. Declare it via \
                `@parameters` in the system so that it can be provided and updated \
                through the parameter object."
            )
        )
    end
    @set! sys.consolidate = weighted_consolidate(weights)
    return sys
end

"""
    constraints_to_penalties(sys::System; weights = 1.0)

Return a new [`System`](@ref) in which every constraint of `sys` - including those of
its subsystems - is removed from `constraints` and appended to `costs` as a weighted
quadratic penalty term. This is the "classical PINN" formulation of a constrained
problem: all constraint residuals are folded into the objective, producing an
unconstrained system that can be solved by optimizers which do not accept explicit
`cons`/`lcons`/`ucons` constraints.

- `Equation` constraints `l ~ r` contribute `weights[i] * (l - r)^2`.
- `Inequality` constraints contribute `weights[i] * max(residual, 0)^2`, where
  `residual` is the constraint rewritten in canonical `residual ≲ 0` form via
  `Symbolics.canonical_form`, so only violations of the constraint are penalized.

`weights` may be a scalar applied to every constraint, or a vector with one entry per
element of `constraints(sys)`, which orders a system's own constraints before its
subsystems'. Symbolic weights (e.g. penalty parameters that should be tunable through
the problem's parameter object) that are not already parameters or unknowns of the
system are automatically added to its parameters.

Penalty terms are appended to `costs(sys)` and are therefore combined with the rest of
the objective through the system's `consolidate` function. Note that a finite `weights`
makes the constraint satisfaction soft: increasing the magnitude of the weights enforces
the constraints more tightly at the cost of a stiffer objective.

# Example

```julia
using ModelingToolkitBase
@variables x
@named sys = OptimizationSystem((x - 2)^2, [x], []; constraints = [x ≲ 1])
pen_sys = constraints_to_penalties(complete(sys); weights = 1.0e3)
```
"""
function constraints_to_penalties(sys::System; weights = 1.0)
    if weights isa Union{AbstractVector, Tuple} &&
            length(weights) != length(constraints(sys))
        throw(
            ArgumentError(
                """
                Expected `weights` to be a scalar or have one entry per constraint of the \
                system (including subsystem constraints). Got $(length(weights)) weights \
                for $(length(constraints(sys))) constraints.
                """
            )
        )
    end
    return _constraints_to_penalties(sys, weights)
end

function _constraints_to_penalties(sys::System, weights)
    cstrs = get_constraints(sys)
    own_weights = weights isa Union{AbstractVector, Tuple} ?
        weights[1:length(cstrs)] : Iterators.repeated(weights, length(cstrs))
    penalties = SymbolicT[]
    new_ps = SymbolicT[]
    for (cstr, w) in zip(cstrs, own_weights)
        w = unwrap(w)
        res = Symbolics.canonical_form(cstr).lhs
        pen = cstr isa Equation ? w * res^2 : w * max(res, 0)^2
        push!(penalties, value(pen))
        _weight_parameters!(new_ps, sys, w)
    end
    @set! sys.costs = [get_costs(sys); penalties]
    @set! sys.constraints = Union{Equation, Inequality}[]

    subsystems = get_systems(sys)
    if !isempty(subsystems)
        newsystems = System[]
        offset = length(cstrs)
        for ssys in subsystems
            subweights = weights
            if weights isa Union{AbstractVector, Tuple}
                nsub = length(constraints(ssys))
                subweights = weights[(offset + 1):(offset + nsub)]
                offset += nsub
            end
            push!(newsystems, _constraints_to_penalties(ssys, subweights))
        end
        @set! sys.systems = newsystems
    end

    if !isempty(new_ps)
        @set! sys.ps = [get_ps(sys); new_ps]
        if has_index_cache(sys) && get_index_cache(sys) !== nothing
            # the `IndexCache` constructor reads `is_parameter` and family, so it must be
            # cleared before rebuilding to avoid seeing the stale cache.
            @set! sys.index_cache = nothing
            @set! sys.index_cache = IndexCache(sys)
        end
    end
    return sys
end

# Collect symbolic variables in a penalty weight that are not already parameters or
# unknowns of `sys`, so they can be added to the system's parameters.
function _weight_parameters!(new_ps::Vector{SymbolicT}, sys::System, w)
    symbolic_type(w) === NotSymbolic() && return new_ps
    for v in get_variables(w)
        symbolic_type(v) === NotSymbolic() && continue
        is_parameter(sys, v) && continue
        any(Base.Fix2(isequal, v), get_unknowns(sys)) && continue
        any(Base.Fix2(isequal, v), new_ps) && continue
        push!(new_ps, v)
    end
    return new_ps
end
