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
    IsLimitedNode()

Predicate matching `limited(actual, limiter)` call nodes. Used both as a `query`
predicate and as an `is_atomic` for `search_variables!`, which then collects the
`limited` nodes themselves rather than the variables inside them.
"""
struct IsLimitedNode end

(::IsLimitedNode)(x) = iscall(x) && operation(x) === limited

"""
    has_limited(expr)

Return `true` iff `expr` contains at least one `limited(...)` node.
"""
function has_limited(expr)
    return SU.query(IsLimitedNode(), unwrap(expr))
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
    collect_limited_nodes!(nodes, ir, exprs)

Collect every distinct `limited(...)` node reachable from `exprs` into the ordered set
`nodes`, in a deterministic order (`limited_k` numbering must be stable across processes
for precompilation). Backed by `IRStructure`, so shared subexpressions of large flattened
models are visited once rather than once per occurrence.
"""
function collect_limited_nodes!(nodes::OrderedSet{SymbolicT}, ir::IRStructure, exprs)
    buffer = SU.IRStructureSearchBuffer(ir, nodes)
    for expr in exprs
        SU.search_variables!(buffer, unwrap(expr); is_atomic = IsLimitedNode())
    end
    return nodes
end

# Rewrite `exprs` (equations) by substituting the given node => replacement rules through
# an `IRSubstituter`, which caches by IR index instead of rebuilding shared subtrees.
function _substitute_equations(ir::IRStructure, rules::Dict{SymbolicT, SymbolicT}, eqs)
    subber = SU.IRSubstituter{false}(ir, rules)
    return [Equation(subber(unwrap(eq.lhs)), subber(unwrap(eq.rhs))) for eq in eqs]
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

struct IsSentinel end

(::IsSentinel)(x) = _sentinel_kind(x) !== :none

# Map namespaced sentinel occurrences back to the toplevel `LIMIT_NEW`/`LIMIT_OLD`, so a
# limiter written inside a subsystem compiles against the same two placeholders. Built as
# a substitution rule set (the sentinels are leaves) and applied with the IR substituter.
function _canonicalize_sentinels(ir::IRStructure, x)
    x = unwrap(x)
    sentinels = OrderedSet{SymbolicT}()
    SU.search_variables!(
        SU.IRStructureSearchBuffer(ir, sentinels), x; is_atomic = IsSentinel()
    )
    isempty(sentinels) && return x
    rules = Dict{SymbolicT, SymbolicT}()
    for s in sentinels
        rules[s] = _sentinel_kind(s) === :new ? LIMIT_NEW : LIMIT_OLD
    end
    return SU.IRSubstituter{false}(ir, rules)(x)
end

"""
    strip_limited_system(sys::AbstractSystem)

Return a copy of `sys` with every `limited(actual, limiter)` node in the equations and
observed equations replaced by `actual`. Used for time-dependent compilation, where
iterate limiting does not apply.
"""
function strip_limited_system(sys::AbstractSystem)
    ir = get_irstructure(sys)
    eqs = equations(sys)
    obs = observed(sys)
    nodes = collect_limited_nodes!(
        OrderedSet{SymbolicT}(), ir,
        Iterators.flatten(((eq.lhs, eq.rhs) for eq in Iterators.flatten((eqs, obs))))
    )
    isempty(nodes) && return sys
    # `limited(actual, limiter) -> actual`; the substituter applies the rules bottom-up,
    # so a (rejected elsewhere, but harmless here) nested annotation still resolves.
    rules = Dict{SymbolicT, SymbolicT}(node => unwrap(arguments(node)[1]) for node in nodes)
    return ConstructionBase.setproperties(
        sys;
        eqs = _substitute_equations(ir, rules, eqs),
        observed = _substitute_equations(ir, rules, obs),
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
    ir = get_irstructure(sys)
    nodes = collect_limited_nodes!(
        OrderedSet{SymbolicT}(), ir,
        Iterators.flatten(((eq.lhs, eq.rhs) for eq in eqs))
    )
    isempty(nodes) && return sys

    rules = Dict{SymbolicT, SymbolicT}()
    specs = Pair{SymbolicT, SymbolicT}[]
    new_vars = SymbolicT[]
    new_eqs = Equation[]
    guessmap = Pair{SymbolicT, SymbolicT}[]
    for (k, node) in enumerate(nodes)
        args = arguments(node)
        actual, limiter = unwrap(args[1]), unwrap(args[2])
        if has_limited(actual) || has_limited(limiter)
            throw(ArgumentError("`limited` operators may not be nested; found $(node)."))
        end
        # `#` makes the generated name unwritable as a Julia identifier, so it cannot
        # collide with a user symbol and no collision guard is needed. The index makes
        # it a pure function of the (deterministically ordered) node list, keeping
        # codegen reproducible across processes.
        name = Symbol("#limited_", k)
        var = unwrap(only(@variables $name [irreducible = true]))
        push!(new_vars, var)
        rules[node] = var
        push!(specs, var => _canonicalize_sentinels(ir, limiter))
        push!(new_eqs, wrap(var) ~ wrap(actual))
        push!(guessmap, var => actual)
    end

    lowered_eqs = _substitute_equations(ir, rules, eqs)

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
    dvset = Set{SymbolicT}(unwrap.(dvs))
    ps = parameters(sys)
    ir = get_irstructure(sys)
    varbuffer = Set{SymbolicT}()
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
        empty!(varbuffer)
        SU.search_variables!(SU.IRStructureSearchBuffer(ir, varbuffer), unwrap(lexpr))
        bad = filter(in(dvset), varbuffer)
        if !isempty(bad)
            throw(
                ArgumentError(
                    "the limiter expression `$(lexpr)` for `$(var)` may only reference " *
                        "`limitnew`, `limitold`, and parameters; it references the " *
                        "unknowns $(collect(bad))."
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
    return iip ? LimitedPostcondition{true, typeof(lims)}(lims) :
        LimitedPostcondition{false, typeof(lims)}(lims)
end

"""
    LimitedPostcondition{iip}(limiters)

The `NonlinearFunction.postcondition` callable generated by
[`generate_limited_postcondition`](@ref). `limiters` is a tuple of
`(index into unknowns, compiled limiter)` pairs; calling it applies each limiter to its
entry of the proposed iterate. A struct rather than a closure so every captured field is
explicit.
"""
struct LimitedPostcondition{iip, L}
    limiters::L
end

function (h::LimitedPostcondition{true})(up, uprev, p)
    for (idx, fn) in h.limiters
        up[idx] = fn((up[idx], uprev[idx]), p)
    end
    return nothing
end

function (h::LimitedPostcondition{false})(up, uprev, p)
    for (idx, fn) in h.limiters
        val = fn((up[idx], uprev[idx]), p)
        up = _limited_setindex(up, val, idx)
    end
    return up
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
