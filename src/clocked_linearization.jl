"""
    $(TYPEDEF)

The linearization of one clock partition of a hybrid system, as returned by
[`linearize_clocked`](@ref).

A partition is either the continuous part of the system, in which case `clock` is a
`ContinuousClock` and the matrices describe

```math
\\begin{aligned}
ẋ &= Ax + Bu \\\\
y &= Cx + Du
\\end{aligned}
```

or one clocked partition, in which case `clock` is a `PeriodicClock` and the matrices
describe the difference equation

```math
\\begin{aligned}
x_{k+1} &= Ax_k + Bu_k \\\\
y_k &= Cx_k + Du_k
\\end{aligned}
```

at the sample interval [`sampletime`](@ref) of that clock.

The inputs and outputs of a partition are its user-requested signals followed by the signals
that cross a clock boundary. `input_groups` and `output_groups` say which is which: the pair
`0 => range` is the user-requested signals in the order they were given, and `j => range` is
the group of signals that comes from, respectively goes to, partition `j`. The same variable
names the signal on both sides of a boundary, so the partitions can be connected by name.

If the partition could not be linearized, `sys`, `A`, `B`, `C`, `D` and `extras` are
`nothing` while the signal lists are still populated, so that a model of that partition
supplied by other means can be connected in its place.

# Fields

$(TYPEDFIELDS)
"""
struct ClockPartition{S, TA, TB, TC, TD, E}
    """
    The clock of this partition. `SciMLBase.ContinuousClock()` for the continuous partition.
    """
    clock::SciMLBase.AbstractClock
    """
    The simplified system this partition was linearized from.
    """
    sys::S
    """
    The dynamics matrix.
    """
    A::TA
    """
    The input matrix.
    """
    B::TB
    """
    The output matrix.
    """
    C::TC
    """
    The feedthrough matrix.
    """
    D::TD
    """
    Input variables, in the column order of `B` and `D`.
    """
    inputs::Vector{SymbolicT}
    """
    Output variables, in the row order of `C` and `D`.
    """
    outputs::Vector{SymbolicT}
    """
    Index ranges into `inputs`. `0 => range` is the user-requested inputs, `j => range` the
    signals entering from partition `j`.
    """
    input_groups::Vector{Pair{Int, UnitRange{Int}}}
    """
    Index ranges into `outputs`. `0 => range` is the user-requested outputs, `j => range`
    the signals leaving for partition `j`.
    """
    output_groups::Vector{Pair{Int, UnitRange{Int}}}
    """
    The operating point the partition was linearized at, `(; x, p, t)`.
    """
    extras::E
end

"""
    $(TYPEDSIGNATURES)

The sample interval of a [`ClockPartition`](@ref), or `nothing` for the continuous
partition.
"""
function sampletime(p::ClockPartition)
    clock = p.clock
    return clock isa SciMLBase.PeriodicClock ? clock.dt : nothing
end

"""
    $(TYPEDSIGNATURES)

Whether `p` is the continuous partition of the system.
"""
is_continuous_partition(p::ClockPartition) = p.clock isa SciMLBase.ContinuousClock

"""
    $(TYPEDSIGNATURES)

Whether `p` was linearized. A partition carrying a clock that is not supported is returned
with its signal lists but without matrices.
"""
is_linearized(p::ClockPartition) = p.A !== nothing

MTKBase.inputs(p::ClockPartition) = p.inputs
MTKBase.outputs(p::ClockPartition) = p.outputs

"""
    $(TYPEDSIGNATURES)

The input variables of `p` that come from partition `j`, or the user-requested inputs for
`j = 0`. Returns an empty vector if there is no such group.
"""
function input_group(p::ClockPartition, j::Int)
    idx = findfirst(x -> first(x) == j, p.input_groups)
    return idx === nothing ? SymbolicT[] : p.inputs[last(p.input_groups[idx])]
end

"""
    $(TYPEDSIGNATURES)

The output variables of `p` that go to partition `j`, or the user-requested outputs for
`j = 0`. Returns an empty vector if there is no such group.
"""
function output_group(p::ClockPartition, j::Int)
    idx = findfirst(x -> first(x) == j, p.output_groups)
    return idx === nothing ? SymbolicT[] : p.outputs[last(p.output_groups[idx])]
end

function Base.show(io::IO, ::MIME"text/plain", p::ClockPartition)
    printstyled(io, "ClockPartition"; bold = true, color = :blue)
    println(io, " on ", p.clock)
    if is_linearized(p)
        println(io, "  state dimension ", size(p.A, 1), ", ", length(p.inputs),
            " inputs, ", length(p.outputs), " outputs")
    else
        println(io, "  not linearized")
    end
    for (label, groups, vars) in (
        ("inputs", p.input_groups, p.inputs), ("outputs", p.output_groups, p.outputs))
        println(io, "  ", label, ":")
        for (j, range) in groups
            isempty(range) && continue
            origin = j == 0 ? "requested" : (label == "inputs" ? "from partition $j" :
                                             "to partition $j")
            println(io, "    ", origin, ": ", join(string.(vars[range]), ", "))
        end
    end
    return nothing
end

"""
Is `τ` a term that carries a signal across a clock boundary? A `Shift` is excluded: it
refers to a past value within the same partition, not to another partition.
"""
_is_clock_crossing(τ) = isoperator(τ, Union{Sample, Hold})

"""
Keep the entries of `op` that name something in `psys`, recording the keys used in `seen`.
An entry naming a parameter is kept for every partition that has it; an entry naming a
variable is kept for the partition that variable belongs to. An entry naming an observed
variable is kept only for a
continuous partition, where the initialization problem can solve for it; in a clocked
partition the state is named by its own past value, for instance `x(k-1)`.
"""
function _partition_op(
        psys::AbstractSystem, op::AbstractDict, belongs, seen::Set; allow_observed::Bool
    )
    out = anydict()
    for (k, v) in op
        key = unwrap(k)
        keep = is_parameter(psys, key) ||
               (belongs(key) &&
                (SymbolicIndexingInterface.is_variable(psys, key) ||
                 (allow_observed && SymbolicIndexingInterface.is_observed(psys, key))))
        keep || continue
        out[key] = v
        push!(seen, key)
    end
    return out
end

"""
The variable an operator expression is built from: `Shift(t, -1)(x)`, `Sample(clk)(x)` and
`Hold()(x)` all name the signal `x`.
"""
function _base_variable(v)
    while isoperator(v, Union{Shift, Sample, Hold})
        v = only(SU.arguments(v))
    end
    return v
end

"""
    partitions = linearize_clocked(sys::AbstractSystem, inputs, outputs; op, t = 0.0, kwargs...)

Linearize every clock partition of `sys` separately and return them as a vector of
[`ClockPartition`](@ref)s. The continuous partition, when the system has one, is last.

`linearize` builds a single continuous linearization and has no representation for a clocked
subsystem, so a model whose controller is a synchronous program linearizes to the plant with
the controller's contribution missing. This function instead splits the model the way the
compiler does, linearizes each partition in its own time domain, and reports the signals that
cross each boundary, leaving the assembly of the loop to the caller. Every partition carries
those crossing signals as inputs and outputs in addition to the requested ones, so the
partitions can be reconnected, for example with `ControlSystems.feedback` or
`RobustAndOptimalControl.connect`.

The continuous partition is returned in continuous time and every clocked partition in
discrete time at its own sample interval. Two routes to a closed-loop model are available and
they are not equivalent. Discretizing the continuous partition at the clock interval (`c2d`
with a zero-order hold) and connecting it to the clocked partitions is exact for a single
periodic clock, because `Hold` is a zero-order hold and `Sample` is ideal sampling.
Converting the clocked partitions to continuous time (`d2c`) instead keeps the resolution of
the fast dynamics, but it is an approximation whose reliability rests on the signals crossing
each boundary having little content above that boundary's Nyquist frequency — a property of
the nonlinear model that no linearization will reveal, and one worth measuring separately.

# Arguments

- `sys`: The system to linearize, **not** simplified. It may contain any number of periodic
  clocks.
- `inputs`: A variable, analysis point, or collection of either, treated as the inputs of the
  linearized model. Each one is assigned to the partition it belongs to.
- `outputs`: The same, for the outputs.

# Keywords

- `op`: The operating point, as for [`linearize`](@ref). A `LinearizationOpPoint` at a single
  time point is accepted. Entries are distributed over the partitions they belong to; an
  entry that names nothing in any partition produces a warning.
- `t = 0.0`: The value of the independent variable at the operating point.
- `allow_input_derivatives = false`: See [`linearize`](@ref).
- `initialize`: Whether to solve an initialization problem for the operating point. Defaults
  to `true` for the continuous partition and `false` for the clocked ones, whose operating
  point is given rather than implied.
- `loop_openings`, `system_modifier`: As for [`linearize`](@ref), when analysis points are
  used.
- `warn_unhandled = true`: Whether to warn about clocks and events whose contribution is not
  in the result.
- `kwargs`: Forwarded to [`linearization_function`](@ref), `autodiff` among them.
  `AutoFiniteDiff()` is the choice for a partition whose update law is not differentiable by
  `ForwardDiff`, such as one that calls an optimization solver.

# Returns

- `partitions::Vector{ClockPartition}`, one per clock, with the continuous partition last.

# Examples

```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

dt = 0.1
k = ShiftIndex(Clock(dt))
@variables x(t) y(t) u(t) yd(t) ud(t) r(t)
@parameters kp

@named sys = System(
    [
        yd ~ Sample(dt)(y)
        ud ~ ud(k - 1) + kp * (r - yd)
        u ~ Hold(ud)
        D(x) ~ -x + u
        y ~ x
    ], t)

partitions = linearize_clocked(sys, [r], [y]; op = Dict(x => 0.0, ud => 0.0, kp => 2.0))
```

See also [`linearize`](@ref).
"""
function linearize_clocked(
        sys::AbstractSystem,
        inputs::Union{Symbol, Vector{Symbol}, AnalysisPoint, Vector{AnalysisPoint}},
        outputs; loop_openings = [], system_modifier = identity, kwargs...
    )
    sys, input_vars, output_vars,
        loop_opening_params = linearization_ap_transform(
        sys, inputs, outputs, loop_openings
    )
    return linearize_clocked(
        system_modifier(sys), input_vars, output_vars; loop_opening_params, kwargs...
    )
end

function linearize_clocked(
        sys::AbstractSystem, inputs, outputs;
        op = Dict{SymbolicT, SymbolicT}(),
        t = 0.0,
        allow_input_derivatives = false,
        initialize = nothing,
        warn_unhandled = true,
        loop_opening_params = SymbolicT[],
        kwargs...
    )
    inputs isa AbstractVector || (inputs = [inputs])
    outputs isa AbstractVector || (outputs = [outputs])
    inputs = SymbolicT[unwrap(x) for x in inputs]
    outputs = SymbolicT[unwrap(x) for x in outputs]

    if op isa LinearizationOpPoint
        op.t isa AbstractVector &&
            throw(ArgumentError("`linearize_clocked` accepts a `LinearizationOpPoint` at a single time point only."))
        op = _build_op_from_solution(op)
    end
    op = anydict(op)

    ets = expand_connections(sys)
    ci = MTKTearing.infer_clocks!(MTKTearing.ClockInference(TearingState(ets)))
    tss, clocked_inputs, continuous_id,
        id_to_clock = MTKTearing.split_system(deepcopy(ci))
    npartitions = length(tss)

    domain_of = Dict{SymbolicT, Any}()
    for (v, d) in zip(ci.ts.fullvars, ci.var_domain)
        domain_of[v] = d
    end
    partition_of_clock = Dict{Any, Int}(d => i for (i, d) in enumerate(id_to_clock))
    function partition_of(v)
        d = get(domain_of, v, nothing)
        d === nothing && (d = get(ci.expression_clocks, v, nothing))
        d === nothing &&
            error("The variable $v provided to `linearize_clocked` was not found in the system, so it could not be assigned to a clock partition.")
        return partition_of_clock[d]
    end

    # Signals that cross a boundary: `(source, term, variable)` entering a partition, and
    # `(destination, variable)` leaving one. The variable names the signal on both sides.
    crossing_in = [Tuple{Int, SymbolicT, SymbolicT}[] for _ in 1:npartitions]
    crossing_out = [Tuple{Int, SymbolicT}[] for _ in 1:npartitions]
    for q in 1:npartitions, τ in clocked_inputs[q]
        _is_clock_crossing(τ) || continue
        v = only(SU.arguments(τ))
        p = partition_of(v)
        p == q && continue
        push!(crossing_in[q], (p, τ, v))
        (q, v) in crossing_out[p] || push!(crossing_out[p], (q, v))
    end

    user_inputs = [SymbolicT[] for _ in 1:npartitions]
    user_outputs = [SymbolicT[] for _ in 1:npartitions]
    for v in inputs
        push!(user_inputs[partition_of(v)], v)
    end
    for v in outputs
        push!(user_outputs[partition_of(v)], v)
    end

    if warn_unhandled
        _warn_unhandled_events(ets)
    end

    seen_op_keys = Set{Any}()
    partitions = ClockPartition[]
    iv = get_iv(ets)
    for i in 1:npartitions
        clock = id_to_clock[i]
        ts = deepcopy(tss[i])
        fullvars = Set{SymbolicT}(ts.fullvars)
        in_vars, in_groups = _partition_signals(user_inputs[i], crossing_in[i])
        out_vars, out_groups = _partition_signals(user_outputs[i], crossing_out[i])

        if !(clock isa SciMLBase.PeriodicClock || clock isa SciMLBase.ContinuousClock)
            warn_unhandled && @warn "The partition on $clock is not a periodic clock and was not linearized. Its contribution is not included in the returned partitions."
            push!(
                partitions,
                ClockPartition(
                    clock, nothing, nothing, nothing, nothing, nothing,
                    in_vars, out_vars, in_groups, out_groups, nothing
                )
            )
            continue
        end

        # Name every boundary signal the way the model does. A crossing signal arrives as the
        # `Sample`/`Hold` term that carries it, and `mark_discrete` has shifted a clocked
        # partition's own variables forward one step; both are operator expressions, which
        # the rest of the compiler does not accept as inputs. Rewriting them to the plain
        # variable leaves an ordinary system whose inputs, outputs and operating point are
        # addressed by the names the model uses.
        discrete = !(clock isa SciMLBase.ContinuousClock)
        ts = _rename_boundary_signals(ts, fullvars, crossing_in[i], iv, discrete)
        fullvars = Set{SymbolicT}(ts.fullvars)
        in_terms = SymbolicT[_partition_spelling(fullvars, v, iv) for v in in_vars]
        out_terms = SymbolicT[_partition_spelling(fullvars, v, iv) for v in out_vars]
        # A user-requested input of a clocked partition is still spelled shifted. It cannot
        # be registered as a system input, since the initialization machinery rejects an
        # operator expression there, so it becomes a plain parameter and is renamed back
        # once the partition is compiled.
        shifted_io = Dict{SymbolicT, SymbolicT}(
            τ => v
            for (v, τ) in Iterators.flatten((zip(in_vars, in_terms), zip(out_vars, out_terms)))
                if !isequal(v, τ)
        )
        plain_inputs = SymbolicT[τ for τ in in_terms if !haskey(shifted_io, τ)]
        shifted_inputs = SymbolicT[τ for τ in in_terms if haskey(shifted_io, τ)]

        psys = _mtkcompile!(
            ts; inputs = OrderedSet{SymbolicT}(plain_inputs),
            discrete_inputs = OrderedSet{SymbolicT}(shifted_inputs),
            outputs = OrderedSet{SymbolicT}(out_terms)
        )
        psys = _unshift_partition_names(psys, shifted_io)
        # A partition with no state variable is a static map, and its matrices are the same
        # in either time domain. `DiscreteProblem` cannot represent a system without
        # unknowns, so such a partition goes through the continuous path; only its reported
        # clock says which rate it runs at.
        stateless = isempty(unknowns(psys))
        discrete && !stateless && (@set! psys.is_discrete = true)
        psys = complete(psys)
        belongs = function (key)
            base = _base_variable(key)
            domain = get(domain_of, base, nothing)
            return domain !== nothing && partition_of_clock[domain] == i
        end
        partition_op = _partition_op(
            psys, op, belongs, seen_op_keys; allow_observed = !discrete)
        linfun, _ = linearization_function(
            psys, in_vars, out_vars;
            already_simplified = true,
            problem_constructor = discrete && !stateless ? DiscreteProblem : ODEProblem,
            # A clocked partition's operating point is given rather than solved for, and a
            # static map has no initialization problem to speak of: its would-be
            # initialization system consists of parameter equations alone.
            problem_kwargs = discrete || stateless ? (; build_initializeprob = false) : (;),
            initialize = initialize === nothing ? !discrete : initialize,
            op = partition_op, t, warn_empty_op = false,
            loop_opening_params = filter(Base.Fix1(is_parameter, psys), loop_opening_params),
            kwargs...
        )
        u0 = state_values(linfun)
        if discrete && u0 !== nothing
            # A clocked partition's state is given, not solved for. Writing it here rather
            # than relying on the problem's initialization keeps the operating point exactly
            # what was asked for; the state is named by its own past value, `x(k-1)`.
            u0 = copy(u0)
            for (key, value) in partition_op
                idx = SymbolicIndexingInterface.variable_index(psys, key)
                idx === nothing || (u0[idx] = value)
            end
        end
        linres = linfun(u0, parameter_values(linfun), t)
        matrices = _linearization_matrices(
            linres, in_vars; allow_input_derivatives)
        push!(
            partitions,
            ClockPartition(
                clock, psys, matrices.A, matrices.B, matrices.C, matrices.D,
                in_vars, out_vars, in_groups, out_groups,
                (; x = linres.x, p = linres.p, t = linres.t)
            )
        )
    end

    unused = setdiff(Set(unwrap.(keys(op))), seen_op_keys)
    if warn_unhandled && !isempty(unused)
        names = join(string.(collect(unused)), ", ")
        @warn "The operating point contains entries that name nothing in any partition and were ignored: $names."
    end
    return partitions
end

"""
How a model-level variable is spelled inside the tearing state of a clock partition.
`mark_discrete` shifts every variable of a clocked partition forward one step, so a variable
that the model calls `v` is called `Shift(t, 1)(v)` there. Which spelling applies is decided
by looking, rather than by assuming, because a variable that only ever appears as a past
value shifts back onto itself.
"""
function _partition_spelling(fullvars::Set{SymbolicT}, v::SymbolicT, iv)
    v in fullvars && return v
    shifted = MTKBase.simplify_shifts(Shift(iv, 1)(v))
    shifted in fullvars && return shifted
    error("The variable $v provided to `linearize_clocked` was not found in its clock partition.")
end

"""
Undo, in the compiled partition, the forward shift that `mark_discrete` applied to the
signals named in `subs`.

Reassembly writes a clocked partition's equations in unshifted form but leaves the parameter
list holding the shifted spelling that an input was declared with, so a variable promoted to
an input ends up referenced by the equations under a name the system does not have. This
renames those parameters and outputs back, leaving the partition consistent and addressable
by the names the model uses.
"""
function _unshift_partition_names(psys::AbstractSystem, subs::AbstractDict)
    isempty(subs) && return psys
    @set! psys.ps = SymbolicT[get(subs, p, p) for p in get_ps(psys)]
    binds = copy(parent(bindings(psys)))
    for (shifted, plain) in subs
        haskey(binds, shifted) || continue
        binds[plain] = binds[shifted]
        delete!(binds, shifted)
    end
    @set! psys.bindings = ROSymmapT(binds)
    @set! psys.outputs = OrderedSet{SymbolicT}(
        SymbolicT[get(subs, o, o) for o in get_outputs(psys)])
    return psys
end

"""
Rewrite the boundary signals of a partition's tearing state to the names the model uses.

A signal arriving from another partition appears as the `Sample` or `Hold` term that carries
it, an operator expression that the compiler does not accept as an input. Substituting it for
the variable it names, in place and before compilation, leaves an ordinary system.

Inside a clocked partition every variable is shifted forward one step by `mark_discrete`, and
reassembly shifts them all back again, so there the substitute must be the shifted variable;
a plain one would end up a step out of step with the rest of the partition. A name is only
taken when it is free.
"""
function _rename_boundary_signals(ts, fullvars::Set{SymbolicT}, crossing, iv, shifted::Bool)
    subs = Dict{SymbolicT, SymbolicT}()
    for (_, term, v) in crossing
        target = shifted ? MTKBase.simplify_shifts(Shift(iv, 1)(v)) : v
        target in fullvars || (subs[term] = target)
    end
    isempty(subs) && return ts

    for (i, v) in enumerate(ts.fullvars)
        ts.fullvars[i] = get(subs, v, v)
    end
    sys = ts.sys
    @set! sys.eqs = Equation[substitute(eq, subs) for eq in get_eqs(sys)]
    @set! sys.unknowns = SymbolicT[get(subs, v, v) for v in get_unknowns(sys)]
    @set! sys.initialization_eqs = Equation[substitute(eq, subs)
                                            for eq in get_initialization_eqs(sys)]
    ts.sys = sys
    ts.original_eqs = Equation[substitute(eq, subs) for eq in ts.original_eqs]
    return ts
end

"""
Order a partition's signals as the user-requested ones followed by one group per neighbouring
partition, and return them with the index range of each group.
"""
function _partition_signals(user, crossing)
    vars = copy(user)
    groups = [0 => 1:length(vars)]
    for j in sort!(unique(first.(crossing)))
        selected = filter(x -> first(x) == j, crossing)
        start = length(vars) + 1
        append!(vars, last.(selected))
        push!(groups, j => start:length(vars))
    end
    return vars, groups
end

"""
Warn about events whose contribution `linearize_clocked` does not include.
"""
function _warn_unhandled_events(sys::AbstractSystem)
    ces = continuous_events(sys)
    isempty(ces) ||
        @warn "The model contains $(length(ces)) continuous event(s). Their contribution is not included in the returned partitions."
    des = discrete_events(sys)
    isempty(des) ||
        @warn "The model contains $(length(des)) discrete event(s). Their contribution is not included in the returned partitions."
    return nothing
end
