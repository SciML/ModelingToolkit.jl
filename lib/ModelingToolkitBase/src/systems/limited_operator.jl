# The `limited(actual, limiter)` operator: symbolic declaration of SPICE-style
# iterate limiting (the predictor/corrector Newton-Raphson (PCNR) formulation of
# Aadithya, Keiter & Mei), lowered to an augmented nonlinear system plus a
# `NonlinearFunction.postcondition` hook. Mirrors the architecture of the Modelica
# `homotopy` operator in `systems/homotopy_operator.jl`.
@register_symbolic limited(actual, limiter)

"""
    limited(actual, limiter)

Declare that the scalar `Real` expression `actual` is a *limited quantity* of a nonlinear
solve: a quantity whose Newton updates must be clipped to a trusted move per iteration,
in the manner of SPICE junction-voltage limiting. `limiter` is a scalar expression in the
reserved placeholder variables [`limitnew`](@ref) and [`limitold`](@ref) — the proposed
and the previously accepted value of the quantity — plus any parameters of the system,
evaluating to the corrected (limited) value.

Following the predictor/corrector Newton-Raphson (PCNR) method of Aadithya, Keiter & Mei,
`mtkcompile` on a time-independent system lowers every `limited` node by introducing an
auxiliary irreducible unknown `limited_k` for the quantity, replacing the node with it,
appending the consistency equation `limited_k ~ actual`, and recording the limiter. The
`SciMLBase.NonlinearFunction` built from the compiled system then carries a generated
`postcondition` hook that applies each limiter to its auxiliary unknown at every iterate a
solver accepts — the corrector phase — while the solver's ordinary Newton step on the
augmented system is the predictor. Residuals and Jacobians are evaluated at the corrected
iterates, so the PCNR consistency property holds. Solving such problems requires a solver
that supports `postcondition` (e.g. `NewtonRaphson` and the other native NonlinearSolve.jl
methods).

For time-dependent systems the operator is stripped to `actual` during `mtkcompile`, so
components carrying limiters compile unchanged for transient simulation; wherever the
operator is evaluated numerically it is simply `actual`.

# Arguments

  - `actual`: the expression being limited. It is what the operator means everywhere
    outside a limited nonlinear solve, and the consistency equation ties the auxiliary
    unknown to it.
  - `limiter`: the correction rule, written in terms of `limitnew`, `limitold`, and
    parameters. It must satisfy `limiter == limitnew` whenever `limitnew == limitold`
    (no proposed movement means no correction) so solutions are fixed points.

# Example

A diode's junction voltage with the classic SPICE3 `pnjlim` limiting rule, registered as
an opaque function so its branches keep Julia short-circuit semantics:

```julia
function pnjlim(vnew, vold, vt, vcrit)
    if vnew > vcrit && abs(vnew - vold) > 2vt
        if vold > 0
            arg = 1 + (vnew - vold) / vt
            vnew = arg > 0 ? vold + vt * log(arg) : vcrit
        else
            vnew = vt * log(vnew / vt)
        end
    end
    return vnew
end
@register_symbolic pnjlim(vnew, vold, vt, vcrit)

@variables v
@parameters Vs R Is Vt vcrit
eqs = [0 ~ (v - Vs) / R + Is * (exp(limited(v, pnjlim(limitnew, limitold, Vt, vcrit)) / Vt) - 1)]
```

**Reference:** K. V. Aadithya, E. R. Keiter, T. Mei, *Predictor/Corrector Newton-Raphson
(PCNR): A Simple, Flexible, Scalable, Modular, and Consistent Replacement for Limiting
in Circuit Simulation*, Scientific Computing in Electrical Engineering, 2020.
"""
limited(actual::Real, limiter::Real) = actual

# The limiter is a solver hint, not part of the residual: Jacobians treat the operator
# as `actual` (matching the runtime fallback and the postcondition hook's contract that
# derivatives do not chain through the corrector). `@register_derivative` requires the
# rule to return a symbolic result, hence the `SConst` wrapping.
Symbolics.@register_derivative limited(actual, limiter) 1 Symbolics.SConst(1)
Symbolics.@register_derivative limited(actual, limiter) 2 Symbolics.SConst(0)

# Reserved placeholder variables for limiter expressions. Fixed sentinel names, NOT
# `gensym`, for the same precompile-safety reasons as `HOMOTOPY_LAMBDA`; the `ₘₜₖ`
# suffix makes user collisions practically impossible. Declared as parameters so that
# variable discovery classifies them as (removable) parameters and the time-dependent
# `System` constructor's is-a-function-of-`t` validation does not reject them.
const LIMIT_NEW = unwrap(only(@parameters __limitnew_ₘₜₖ))
const LIMIT_OLD = unwrap(only(@parameters __limitold_ₘₜₖ))

"""
    limitnew

Reserved placeholder for the *proposed* value of a limited quantity inside the `limiter`
argument of [`limited`](@ref).
"""
const limitnew = wrap(LIMIT_NEW)

"""
    limitold

Reserved placeholder for the *previously accepted* value of a limited quantity inside the
`limiter` argument of [`limited`](@ref).
"""
const limitold = wrap(LIMIT_OLD)

"""
    has_limited(expr)

Return `true` iff `expr` contains at least one `limited(...)` node.
"""
function has_limited(expr)
    x = unwrap(expr)
    return _has_limited(x)
end

function _has_limited(x)
    if !iscall(x)
        return false
    end
    if operation(x) === limited
        return true
    end
    return any(_has_limited, arguments(x))
end

"""
    has_limited_in_equations(eqs)

Return `true` iff any equation in `eqs` (lhs or rhs) contains a `limited(...)` node.
"""
function has_limited_in_equations(eqs)
    for eq in eqs
        if has_limited(eq.lhs) || has_limited(eq.rhs)
            return true
        end
    end
    return false
end

"""
    has_any_limited(sys)

Return `true` iff `sys` contains a `limited(...)` node anywhere in its equations or its
observed equations.
"""
function has_any_limited(sys)
    has_limited_in_equations(equations(sys)) && return true
    return any(eq -> has_limited(eq.lhs) || has_limited(eq.rhs), observed(sys))
end

"""
    _strip_limited(x)

Recursively replace every `limited(actual, limiter)` node in the unwrapped expression `x`
with `actual`, discarding the limiter (and with it any `limitnew`/`limitold` sentinels).
"""
function _strip_limited(x)
    if !iscall(x)
        return x
    end
    op = operation(x)
    args = arguments(x)
    if op === limited
        return _strip_limited(args[1])
    end
    new_args = map(_strip_limited, args)
    return maketerm(typeof(x), op, new_args, metadata(x))
end

# The `limitnew`/`limitold` sentinels appear inside limiter arguments, so `System`'s
# variable discovery collects them as unknowns — namespaced (`diode₊__limitnew_ₘₜₖ`)
# when the `limited` call lives in a subsystem. Sentinels are therefore recognized by
# their name suffix, and limiter expressions are canonicalized back to the toplevel
# sentinels during lowering.
function _sentinel_kind(x)
    x = unwrap(x)
    iscall(x) && return :none
    hasname(x) || return :none
    s = string(getname(x))
    endswith(s, string(getname(LIMIT_NEW))) && return :new
    endswith(s, string(getname(LIMIT_OLD))) && return :old
    return :none
end

function _remove_sentinels(vars)
    return filter(v -> _sentinel_kind(v) === :none, vars)
end

function _canonicalize_sentinels(x)
    kind = _sentinel_kind(x)
    kind === :new && return LIMIT_NEW
    kind === :old && return LIMIT_OLD
    iscall(x) || return x
    new_args = map(_canonicalize_sentinels, arguments(x))
    return maketerm(typeof(x), operation(x), new_args, metadata(x))
end

"""
    strip_limited_system(sys::AbstractSystem)

Return a copy of `sys` with every `limited(actual, limiter)` node in the equations and
observed equations replaced by `actual`. Used for time-dependent compilation, where
iterate limiting does not apply.
"""
function strip_limited_system(sys::AbstractSystem)
    rew(x) = _strip_limited(unwrap(x))
    new_eqs = [Equation(rew(eq.lhs), rew(eq.rhs)) for eq in equations(sys)]
    new_obs = [Equation(rew(eq.lhs), rew(eq.rhs)) for eq in observed(sys)]
    return ConstructionBase.setproperties(
        sys; eqs = new_eqs, observed = new_obs,
        unknowns = _remove_sentinels(unknowns(sys)),
        ps = _remove_sentinels(get_ps(sys))
    )
end

"""
    $(TYPEDSIGNATURES)

Metadata key under which [`lower_limited`](@ref) records the limiter registry: a
`Vector{Pair{SymbolicT, SymbolicT}}` mapping each auxiliary limited unknown to its
limiter expression (in terms of `LIMIT_NEW`/`LIMIT_OLD` and parameters).
"""
struct LimitedCtx end

# Replace every node in `x` that appears as a key of `repl` (an ordered
# node => variable map keyed by `isequal`) with its replacement.
function _replace_limited_nodes(x, repl)
    if !iscall(x)
        return x
    end
    for (node, var) in repl
        if isequal(x, node)
            return var
        end
    end
    new_args = map(Base.Fix2(_replace_limited_nodes, repl), arguments(x))
    return maketerm(typeof(x), operation(x), new_args, metadata(x))
end

function _collect_limited_nodes!(nodes, x)
    if !iscall(x)
        return nodes
    end
    if operation(x) === limited
        if !any(Base.Fix1(isequal, x), nodes)
            push!(nodes, x)
        end
        # `limited` nodes may not be nested; the actual/limiter arguments are searched
        # so nesting is caught by the guard in `lower_limited`.
    end
    foreach(arg -> _collect_limited_nodes!(nodes, arg), arguments(x))
    return nodes
end

"""
    lower_limited(sys::AbstractSystem)

PCNR lowering pass for time-independent systems, run on the flattened system during
`mtkcompile` before structural simplification. For each unique `limited(actual, limiter)`
node in the equations:

 1. introduce an auxiliary irreducible unknown `limited_k` (irreducible so structural
    simplification cannot eliminate the augmentation),
 2. replace the node with `limited_k` everywhere,
 3. append the consistency equation `limited_k ~ actual`,
 4. add the symbolic guess `limited_k => actual` so initial values are consistent,
 5. record `limited_k => limiter` in the [`LimitedCtx`](@ref) metadata for
    `NonlinearFunction` construction to compile into a `postcondition` hook.

Returns the augmented system.
"""
function lower_limited(sys::AbstractSystem)
    eqs = equations(sys)
    nodes = SymbolicT[]
    for eq in eqs
        _collect_limited_nodes!(nodes, unwrap(eq.lhs))
        _collect_limited_nodes!(nodes, unwrap(eq.rhs))
    end
    isempty(nodes) && return sys

    # Guard against collisions with the auxiliary variable names about to be created.
    # (The sentinels themselves are legitimately present: variable discovery collects
    # them from the limiter arguments; they are removed from the unknowns below.)
    reserved = Set{Symbol}()
    for k in eachindex(nodes)
        push!(reserved, Symbol(:limited_, k))
    end
    for v in Iterators.flatten((unknowns(sys), parameters(sys; initial_parameters = true)))
        if hasname(v) && getname(v) in reserved
            throw(
                ArgumentError(
                    "the symbol `$(getname(v))` collides with a name reserved by the " *
                        "`limited` operator lowering; rename the variable."
                )
            )
        end
    end

    repl = Pair{SymbolicT, SymbolicT}[]
    specs = Pair{SymbolicT, SymbolicT}[]
    new_vars = SymbolicT[]
    new_eqs = Equation[]
    guessmap = Pair{SymbolicT, SymbolicT}[]
    for (k, node) in enumerate(nodes)
        args = arguments(node)
        actual, limiter = unwrap(args[1]), unwrap(args[2])
        if _has_limited(actual) || _has_limited(limiter)
            throw(ArgumentError("`limited` operators may not be nested; found $(node)."))
        end
        name = Symbol(:limited_, k)
        var = unwrap(only(@variables $name [irreducible = true]))
        push!(new_vars, var)
        push!(repl, node => var)
        push!(specs, var => _canonicalize_sentinels(limiter))
        push!(new_eqs, wrap(var) ~ wrap(actual))
        push!(guessmap, var => actual)
    end

    rew(x) = _replace_limited_nodes(unwrap(x), repl)
    lowered_eqs = [Equation(rew(eq.lhs), rew(eq.rhs)) for eq in eqs]

    # The consistency value `actual` is written both as a default (consumed by plain
    # `u0` construction) and as a guess (consumed by initialization machinery), so the
    # auxiliary unknowns never require user-provided initial values in either pipeline.
    newguesses = copy(get_guesses(sys))
    newics = copy(get_initial_conditions(sys))
    for (var, guess) in guessmap
        newguesses[var] = guess
        newics[var] = guess
    end
    # Keep the constructor-built caches coherent for the injected auxiliary variables:
    # `var_to_name` resolves names to the metadata-carrying instances, and structural
    # simplification reads irreducibility from the `irreducibles` field (which the
    # `System` constructor collects from variable metadata), not from per-variable
    # metadata alone.
    new_var_to_name = copy(get_var_to_name(sys))
    new_irreducibles = copy(get_irreducibles(sys))
    for var in new_vars
        new_var_to_name[getname(var)] = var
        push!(new_irreducibles, var)
    end
    newsys = ConstructionBase.setproperties(
        sys;
        eqs = [lowered_eqs; new_eqs],
        unknowns = [_remove_sentinels(unknowns(sys)); new_vars],
        ps = _remove_sentinels(get_ps(sys)),
        guesses = newguesses,
        initial_conditions = newics,
        var_to_name = new_var_to_name,
        irreducibles = new_irreducibles,
    )
    return setmetadata(newsys, LimitedCtx, specs)
end

"""
    apply_limited_lowering(sys::AbstractSystem)

The `mtkcompile` entry point for the [`limited`](@ref) operator, called on the flattened
system before structural simplification: time-independent systems are lowered via
[`lower_limited`](@ref) (PCNR augmentation), time-dependent systems are stripped via
[`strip_limited_system`](@ref). Systems without `limited` nodes pass through untouched.
"""
function apply_limited_lowering(sys::AbstractSystem)
    has_any_limited(sys) || return sys
    if is_time_dependent(sys)
        return strip_limited_system(sys)
    else
        return lower_limited(sys)
    end
end

"""
    apply_limited_lowering(sys::AbstractSystem, source_info)

Variant for compilation pipelines that track per-equation
[`EquationSourceInformation`](@ref): returns `(sys, source_info)` with the source
information padded for the consistency equations [`lower_limited`](@ref) appends
(unknown-source entries, not connection equations).
"""
function apply_limited_lowering(sys::AbstractSystem, source_info)
    n_before = length(equations(sys))
    sys = apply_limited_lowering(sys)
    if source_info !== nothing
        n_added = length(equations(sys)) - n_before
        if n_added > 0
            source_info = EquationSourceInformation(
                vcat(source_info.eqs_source, [Symbol[] for _ in 1:n_added]),
                vcat(source_info.is_connection_equation, falses(n_added))
            )
        end
    end
    return sys, source_info
end

"""
    generate_limited_postcondition(sys::AbstractSystem, iip::Bool; kwargs...)

Compile the limiter registry recorded by [`lower_limited`](@ref) into a
`NonlinearFunction.postcondition` hook `H(u_proposed, u_prev, p)` (mutating `u_proposed`
when `iip`). Returns `nothing` when the system carries no limiters. Each limiter is
compiled with [`generate_custom_function`](@ref) as a scalar function of
`(limitnew, limitold)` and the system's parameters; the hook applies it to the auxiliary
limited unknowns at their indices in `unknowns(sys)`.
"""
function generate_limited_postcondition(
        sys::AbstractSystem, iip::Bool;
        eval_expression = false, eval_module = @__MODULE__, kwargs...
    )
    specs = getmetadata(sys, LimitedCtx, nothing)
    specs === nothing && return nothing
    dvs = unknowns(sys)
    udvs = unwrap.(dvs)
    ps = parameters(sys)
    lims = map(specs) do (var, lexpr)
        idx = findfirst(isequal(var), dvs)
        if idx === nothing
            error(
                "the limited auxiliary unknown `$(var)` is not an unknown of the " *
                    "compiled system; this should not happen since it is marked " *
                    "irreducible. Please open an issue."
            )
        end
        # Limiters are functions of the proposed/previous value and (possibly bound)
        # parameters only; referencing other unknowns is not supported since the hook
        # sees the limited entries, not the whole state.
        bad = [
            v for v in unwrap.(Symbolics.get_variables(lexpr))
                if any(Base.Fix1(isequal, v), udvs)
        ]
        if !isempty(bad)
            throw(
                ArgumentError(
                    "the limiter expression `$(lexpr)` for `$(var)` may only reference " *
                        "`limitnew`, `limitold`, and parameters; it references the " *
                        "unknowns $(bad)."
                )
            )
        end
        fn = generate_custom_function(
            sys, lexpr, [LIMIT_NEW, LIMIT_OLD], ps;
            expression = Val{false}, eval_expression, eval_module
        )
        (idx, fn)
    end
    lims = Tuple(lims)
    if iip
        return let lims = lims
            function limited_postcondition_iip(up, uprev, p)
                for (idx, fn) in lims
                    up[idx] = fn((up[idx], uprev[idx]), p)
                end
                return nothing
            end
        end
    else
        return let lims = lims
            function limited_postcondition_oop(up, uprev, p)
                for (idx, fn) in lims
                    val = fn((up[idx], uprev[idx]), p)
                    up = _limited_setindex(up, val, idx)
                end
                return up
            end
        end
    end
end

function _limited_setindex(u, val, idx)
    if ArrayInterface.ismutable(u)
        u2 = copy(u)
        u2[idx] = val
        return u2
    else
        return Base.setindex(u, val, idx)
    end
end
