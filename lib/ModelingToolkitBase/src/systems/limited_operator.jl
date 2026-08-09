# The `limited(actual, limiter)` operator: symbolic declaration of SPICE-style
# iterate limiting (the predictor/corrector Newton-Raphson (PCNR) formulation of
# Aadithya, Keiter & Mei), lowered to an augmented nonlinear system plus a
# `postcondition` corrector attached to the problem's solver options. Mirrors the
# architecture of the Modelica
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
problem built from the compiled system then carries a generated `postcondition` solver
option — forwarded to `solve`/`init` like any other keyword — that applies each limiter to
its auxiliary unknown at every iterate a solver accepts — the corrector phase — while the solver's ordinary Newton step on the
augmented system is the predictor. Residuals and Jacobians are evaluated at the corrected
iterates, so the PCNR consistency property holds. Solving such problems requires a solver
that supports `postcondition` (e.g. `NewtonRaphson` and the other native NonlinearSolve.jl
methods).

For time-dependent systems the operator is stripped to `actual` during `mtkcompile`, so
components carrying limiters compile unchanged for transient simulation; wherever the
operator is evaluated numerically it is simply `actual`. The limiters are remembered
though: an `ODEProblem` built with `nlstep = true` re-attaches them to the nonlinear system
of implicit stage equations, so an implicit solver limits its stage Newton iterates. There
`limitold` is the previous *Newton iterate of the stage solve*, not the previous time step
— which is the SPICE reading of limiting, since limiting damps the iteration rather than
the trajectory. This requires each limited quantity to be an affine function of a single
stage unknown.

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
`limiter` argument of [`limited`](@ref) — the previous iterate of the nonlinear solve, which
in an `nlstep` transient solve is the previous Newton iterate of the current implicit stage.
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
    $(TYPEDSIGNATURES)

Metadata key under which [`strip_limited_system`](@ref) records the limited quantities of a
time-dependent system: a `Vector{Pair{SymbolicT, SymbolicT}}` mapping each limited
expression `actual` to its limiter (in terms of `LIMIT_NEW`/`LIMIT_OLD` and parameters).

The nodes themselves are gone from the equations, so a transient compilation and the
`ODEFunction` it generates are exactly what they would be without the annotation. The
registry is what lets [`attach_stage_limiters`](@ref) put the limiting back for the
implicit stage solves of an `nlstep` problem.
"""
struct StageLimitedCtx end

"""
    strip_limited_system(sys::AbstractSystem)

Return a copy of `sys` with every `limited(actual, limiter)` node in the equations and
observed equations replaced by `actual`, recording the stripped quantities in the
[`StageLimitedCtx`](@ref) metadata. Used for time-dependent compilation, where the
right-hand side itself carries no limiting: a transient solve limits (if at all) inside
the nonlinear solve of an implicit step, which is a separate system built later.
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
    specs = Pair{SymbolicT, SymbolicT}[
        unwrap(arguments(node)[1]) => _canonicalize_sentinels(ir, unwrap(arguments(node)[2]))
            for node in nodes
    ]
    newsys = ConstructionBase.setproperties(
        sys;
        eqs = _substitute_equations(ir, rules, eqs),
        observed = _substitute_equations(ir, rules, obs),
        unknowns = _remove_sentinels(unknowns(sys)),
        ps = _remove_sentinels(get_ps(sys))
    )
    return setmetadata(newsys, StageLimitedCtx, specs)
end

"""
    $(TYPEDSIGNATURES)

Metadata key under which [`lower_limited`](@ref) records the limiter registry: a
`Vector{Pair{SymbolicT, SymbolicT}}` mapping each auxiliary limited unknown to its
limiter expression (in terms of `LIMIT_NEW`/`LIMIT_OLD` and parameters).

A registry entry is a correction, not a property of the unknown, so the same unknown may
appear more than once: [`generate_limited_postcondition`](@ref) applies the entries in
order, each seeing the previous one's output as `limitnew`, so several correctors on one
unknown compose rather than overwrite.
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
    $(TYPEDSIGNATURES)

The clamping limiters implied by the `bounds` metadata of the unknowns of `sys`: a
`Vector{Pair{SymbolicT, SymbolicT}}` sending each unknown with a finite bound to
`max(lo, min(limitnew, hi))`, with only the finite sides present.

These are ordinary limiters of the *physical* quantity, which is the point: handed to
[`attach_stage_limiters`](@ref) they are conjugated into stage coordinates like any other
limiter, turning the static box on the state into the per-stage box on the increment the
stage solver actually iterates on.

Unknowns whose bounds are not scalar are skipped, since a corrector acts on a single entry
of the iterate.
"""
function bounds_limiter_specs(sys::AbstractSystem)
    specs = Pair{SymbolicT, SymbolicT}[]
    for v in unknowns(sys)
        v = unwrap(v)
        lim = _clamp_limiter(getbounds(v)...)
        lim === nothing && continue
        push!(specs, v => lim)
    end
    return specs
end

# `min`/`max` rather than `clamp` so that a one-sided bound contributes a single comparison
# and an infinite side contributes nothing to the corrector at all.
function _clamp_limiter(lo, hi)
    lo, hi = unwrap(lo), unwrap(hi)
    haslo, hashi = _is_active_bound(lo), _is_active_bound(hi)
    haslo || hashi || return nothing
    expr = limitnew
    hashi && (expr = min(expr, wrap(hi)))
    haslo && (expr = max(expr, wrap(lo)))
    return unwrap(expr)
end

# Symbolic bounds are resolved against the parameters at call time, so they are always
# active; a numeric bound is active only when it actually constrains anything.
_is_active_bound(b) = b isa SymbolicT || (b isa Number && isfinite(b))

"""
    attach_stage_limiters(sys::AbstractSystem, stagesys::AbstractSystem, subrules::AbstractDict; limit_bounds = true)

Re-attach the limiters that [`strip_limited_system`](@ref) recorded for the time-dependent
system `sys` to `stagesys`, the time-independent system of implicit stage equations built
from it for `nlstep` problems (`M z = outer_tmp + γ₁ f(γ₂ z + inner_tmp, p, c)`).
`subrules` is the substitution that took the unknowns of `sys` to the stage coordinates of
`stagesys`, i.e. `v => γ₂ * v + inner_tmp[i]` plus `t => c`.

The stage system's unknowns must stay a subset of the ODE unknowns — `SciMLBase.ODENLStepData`
maps the two by index — so this cannot use the auxiliary-unknown augmentation
[`lower_limited`](@ref) applies to a standalone nonlinear system. Instead, each limited
quantity `q` must resolve to an affine function `q = a * z + b` of a *single* stage unknown
`z` (with `a` and `b` free of unknowns), which is what makes the correction on `q` a
correction on `z`: the recorded limiter `L` is conjugated by that affine map into

    φ(znew, zold) = (L(a * znew + b, a * zold + b) - b) / a

and registered against `z` in the [`LimitedCtx`](@ref) metadata, so the ordinary
[`merge_limited_postcondition`](@ref) path compiles it into the stage `NonlinearProblem`'s
`postcondition`. `φ` inherits the fixed-point property of `L`, and since the corrector runs
before the residual is evaluated at each stage iterate, the converged stage solution solves
the unmodified stage equations. `a` and `b` are expressions in the stage parameters, so the
compiled corrector reads them live and is correct for every stage, which a static bound
vector on the stage problem could never be.

Unless `limit_bounds = false`, the `bounds` metadata of the unknowns of `sys` contributes
clamping correctors the same way (see [`bounds_limiter_specs`](@ref)). They are registered
after the model's own limiters, so a quantity carrying both is corrected as
`clamp ∘ limiter` and the corrected iterate always lands in range. Clamping an intermediate
Newton iterate cannot change the root the stage solve converges to, only the path it takes
there, and it keeps `log`/`sqrt`/`exp` right-hand sides from being evaluated outside their
domain — so it is a safety property and applies by default.

Returns `stagesys` unchanged when there is nothing to attach. A *declared* limited quantity
that does not reduce to a single stage unknown affinely is an `ArgumentError`: the user
asked for limiting and cannot have it. The bounds-derived clamps are a best-effort safety
net rather than a request, so one that cannot be conjugated is dropped instead — a model
that integrates today must not stop integrating because a state carries a box.
"""
function attach_stage_limiters(
        sys::AbstractSystem, stagesys::AbstractSystem, subrules::AbstractDict;
        limit_bounds::Bool = true
    )
    specs = getmetadata(sys, StageLimitedCtx, nothing)
    todo = Tuple{SymbolicT, SymbolicT, Symbol}[]
    if specs !== nothing
        append!(todo, ((actual, limiter, :limited) for (actual, limiter) in specs))
    end
    if limit_bounds
        append!(todo, ((v, lim, :bounds) for (v, lim) in bounds_limiter_specs(sys)))
    end
    isempty(todo) && return stagesys
    dvset = Set{SymbolicT}(unwrap.(unknowns(stagesys)))
    stage_specs = Pair{SymbolicT, SymbolicT}[]
    for (actual, limiter, origin) in todo
        if has_limited(actual) || has_limited(limiter)
            throw(
                ArgumentError("`limited` operators may not be nested; found `$(actual)`.")
            )
        end
        what = _limiter_origin(actual, origin)
        # `actual` is written in the variables the system had before simplification, so
        # resolve it through the variables each of the two compilations eliminated, and
        # through the stage substitution in between.
        q = unwrap(substitute_observed(sys, unwrap(actual)))
        q = unwrap(substitute(q, subrules))
        q = unwrap(substitute_observed(stagesys, q))
        zs = unique(filter(in(dvset), unwrap.(get_variables(q))))
        # A limited quantity that no stage unknown feeds is constant throughout the stage
        # solve — the solver never iterates on it, so there is nothing to limit.
        isempty(zs) && continue
        if length(zs) != 1
            origin === :bounds && continue
            throw(
                ArgumentError(
                    "`nlstep = true` requires every limited quantity to be a function of " *
                        "exactly one unknown of the implicit stage system, since the stage " *
                        "unknowns are the only iterate the corrector can act on. The " *
                        "corrector from $(what) applies to `$(q)`, which depends on the " *
                        "stage unknowns $(zs). Limit a single state instead, or build the " *
                        "problem with `nlstep = false`."
                )
            )
        end
        z = only(zs)
        a = unwrap(Symbolics.derivative(q, z))
        if any(isequal(z), unwrap.(get_variables(a)))
            origin === :bounds && continue
            throw(
                ArgumentError(
                    "`nlstep = true` requires every limited quantity to be an affine " *
                        "function of its stage unknown, so that limiting the quantity is " *
                        "limiting the unknown. The corrector from $(what) applies to " *
                        "`$(q)`, which is nonlinear in `$(z)`."
                )
            )
        end
        b = unwrap(substitute(q, Dict(z => 0)))
        push!(stage_specs, z => _conjugate_limiter(limiter, a, b))
    end
    isempty(stage_specs) && return stagesys
    return setmetadata(stagesys, LimitedCtx, stage_specs)
end

function _limiter_origin(actual, origin::Symbol)
    origin === :bounds && return "the `bounds` metadata of `$(actual)`"
    return "the `limited` quantity `$(actual)`"
end

# `L` limits the physical quantity `a * z + b`; the stage solver only ever sees `z`, so the
# corrector it needs is `L` conjugated by that affine map. `Substituter` rewrites each node
# at most once, so the sentinels reintroduced by the replacements are not rewritten again.
function _conjugate_limiter(limiter, a, b)
    rules = Dict{SymbolicT, SymbolicT}(
        LIMIT_NEW => unwrap(a * limitnew + b),
        LIMIT_OLD => unwrap(a * limitold + b)
    )
    subber = SU.Substituter{false}(rules, SU.default_substitute_filter)
    return unwrap((wrap(subber(unwrap(limiter))) - b) / a)
end

"""
    merge_limited_postcondition(sys, iip, kwargs; expression, kwargs...)

Attach the corrector compiled by [`generate_limited_postcondition`](@ref) to the problem
keywords as `postcondition`, so it is forwarded to `solve`/`init` like any other solver
option. Returns `kwargs` unchanged when `sys` carries no `limited` quantities, and leaves
a user-supplied `postcondition` in place (an explicitly passed corrector wins over the
model-declared one, matching the precedence of solver options generally).
"""
function merge_limited_postcondition(
        sys::AbstractSystem, iip::Bool, kwargs;
        expression = Val{false}, eval_expression = false, eval_module = @__MODULE__
    )
    getmetadata(sys, LimitedCtx, nothing) === nothing && return kwargs
    haskey(kwargs, :postcondition) && return kwargs
    if expression == Val{true}
        throw(
            ArgumentError(
                "`expression = Val{true}` is not supported for systems with `limited` " *
                    "quantities; the generated `postcondition` corrector is a runtime " *
                    "closure."
            )
        )
    end
    post = generate_limited_postcondition(sys, iip; eval_expression, eval_module)
    post === nothing && return kwargs
    return merge(NamedTuple(kwargs), (; postcondition = post))
end

"""
    generate_limited_postcondition(sys::AbstractSystem, iip::Bool; kwargs...)

Compile the limiter registry recorded by [`lower_limited`](@ref) into a
`postcondition` corrector `H(u_proposed, u_prev, p, cache)` (mutating `u_proposed` when
`iip`; the solver cache is unused by generated limiters). Returns `nothing` when the system carries no limiters. Each limiter is
compiled with [`generate_custom_function`](@ref) as a scalar function of
`(limitnew, limitold)` and the system's parameters; the hook applies it to the auxiliary
limited unknowns at their indices in `unknowns(sys)`, in registry order and each on the
running value, so several correctors registered against one unknown compose.
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
entry of the proposed iterate. The solver-cache argument is accepted for the corrector
contract and ignored, since limiters are functions of the iterate and parameters only. A struct rather than a closure so every captured field is
explicit.
"""
struct LimitedPostcondition{iip, L}
    limiters::L
end

function (h::LimitedPostcondition{true})(up, uprev, p, cache)
    for (idx, fn) in h.limiters
        up[idx] = fn((up[idx], uprev[idx]), p)
    end
    return nothing
end

function (h::LimitedPostcondition{false})(up, uprev, p, cache)
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
