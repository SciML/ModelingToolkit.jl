"""
    $(TYPEDEF)

Linearization of a single clock partition of a hybrid system, as returned in the
`partitions` field of a [`HybridLinearization`](@ref).

For the continuous partition the matrices describe

```math
\\begin{aligned}
ẋ &= Ax + Bu \\\\
y &= Cx + Du
\\end{aligned}
```

For a discrete partition with sample interval `Ts`, the state consists of the values of the
discrete variables at the previous tick, and the matrices describe the update performed at a
tick,

```math
\\begin{aligned}
x(k) &= A x(k-1) + B u(k) \\\\
y(k) &= C x(k-1) + D u(k)
\\end{aligned}
```

whose transfer function is `C(zI - A)⁻¹B + D`. The variables in `unknowns` are the history
variables (`xₜ₋₁` denotes the value of `x` at the previous tick).

The inputs of a partition consist of the user-specified inputs that belong to the partition
(first, in the order given by the user) followed by the signals that enter the partition
across a clock boundary, grouped by the partition they originate from. The outputs consist of
the user-specified outputs that belong to the partition (first) followed by the signals that
leave the partition across a clock boundary. Boundary inputs are the operator terms of the
model, such as `Sample(dt)(y)` or `Hold(u)`; boundary outputs are the variables these
operators are applied to.

# Fields

$(TYPEDFIELDS)
"""
struct ClockPartitionLinearization
    """The state matrix."""
    A::Matrix{Float64}
    """The input matrix."""
    B::Matrix{Float64}
    """The output matrix."""
    C::Matrix{Float64}
    """The feedthrough matrix."""
    D::Matrix{Float64}
    """The state variables, in the order of the rows and columns of `A`."""
    unknowns::Vector{SymbolicT}
    """The input variables, in the order of the columns of `B` and `D`."""
    inputs::Vector{SymbolicT}
    """The output variables, in the order of the rows of `C` and `D`."""
    outputs::Vector{SymbolicT}
    """The number of user-specified inputs; `inputs[1:nu_user]` is the user input group."""
    nu_user::Int
    """The number of user-specified outputs; `outputs[1:ny_user]` is the user output group."""
    ny_user::Int
    """The clock of the partition, `ContinuousClock()` or a `PeriodicClock`."""
    clock::TimeDomain
    """The sample interval of a periodic clock, `nothing` for the continuous partition or for a clock without a fixed interval."""
    Ts::Union{Nothing, Float64}
    """The state at the operating point, in the order of `unknowns`."""
    x0::Vector{Float64}
    """
    The compiled system of the partition. Its unknowns and parameters carry the operating
    point. A boundary term such as `Hold(u)` is represented in this system by a plain
    parameter named after the signal and the operator, `u_hold`.
    """
    sys::System
end

"""
    $(TYPEDEF)

The result of [`linearize_hybrid`](@ref): one [`ClockPartitionLinearization`](@ref) per clock
partition of the model together with the connections between them.

Each connection `(; from, output, to, input)` states that `partitions[from].outputs[output]`
drives `partitions[to].inputs[input]` across a clock boundary. These are the index pairs
required to assemble the sampled-data loop with the advanced interface of
`ControlSystemsBase.feedback`, or the signal names required by
`RobustAndOptimalControl.connect`, once each partition has been converted to a common time
domain (for example with `c2d` or `d2c`). One output may drive several inputs.

# Fields

$(TYPEDFIELDS)
"""
struct HybridLinearization
    """The linearized partitions. The continuous partition, if any, is first."""
    partitions::Vector{ClockPartitionLinearization}
    """The connections across clock boundaries."""
    connections::Vector{@NamedTuple{from::Int, output::Int, to::Int, input::Int}}
    """The index of the continuous partition in `partitions`, `nothing` for a purely discrete model."""
    continuous_index::Union{Nothing, Int}
end

"""
    clock_boundary(hl::HybridLinearization, from::Integer, to::Integer)

The signals crossing the clock boundary from partition `from` to partition `to` of `hl`, as
a named tuple `(; from, to, outputs, inputs)`. `outputs[i]` is the index into
`hl.partitions[from].outputs` of the `i`-th signal leaving partition `from`, and `inputs[i]`
is the index into `hl.partitions[to].inputs` of the entry it drives. With
`hl.partitions[from]` as the first system and `hl.partitions[to]` as the second, the keyword
arguments `Y1 = outputs, U2 = inputs` of the advanced interface of
`ControlSystemsBase.feedback` connect these signals. Both index vectors are empty if no
signal crosses from `from` to `to`.
"""
function clock_boundary(hl::HybridLinearization, from::Integer, to::Integer)
    outputs = Int[]
    inputs = Int[]
    for c in hl.connections
        c.from == from && c.to == to || continue
        push!(outputs, c.output)
        push!(inputs, c.input)
    end
    return (; from = Int(from), to = Int(to), outputs, inputs)
end

function Base.show(io::IO, ::MIME"text/plain", hl::HybridLinearization)
    n = length(hl.partitions)
    printstyled(io, "HybridLinearization"; bold = true, color = :blue)
    println(io, " with ", n, " clock partition", n == 1 ? "" : "s")
    for (i, p) in enumerate(hl.partitions)
        print(io, "  [", i, "] ")
        if p.clock isa SciMLBase.ContinuousClock
            print(io, "continuous")
        elseif p.Ts === nothing
            print(io, p.clock)
        else
            print(io, "discrete, Ts = ", p.Ts)
            iszero(p.clock.phase) || print(io, ", phase = ", p.clock.phase)
        end
        println(
            io, ": state dimension ", size(p.A, 1), ", ",
            length(p.inputs), " input", length(p.inputs) == 1 ? "" : "s",
            " (", p.nu_user, " user), ",
            length(p.outputs), " output", length(p.outputs) == 1 ? "" : "s",
            " (", p.ny_user, " user)"
        )
    end
    isempty(hl.connections) && return nothing
    println(io, "Connections:")
    for c in hl.connections
        from = hl.partitions[c.from]
        to = hl.partitions[c.to]
        println(
            io, "  [", c.from, "].outputs[", c.output, "] ", from.outputs[c.output],
            " => [", c.to, "].inputs[", c.input, "] ", to.inputs[c.input]
        )
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", p::ClockPartitionLinearization)
    printstyled(io, "ClockPartitionLinearization"; bold = true, color = :blue)
    if p.clock isa SciMLBase.ContinuousClock
        println(io, " of the continuous partition")
    else
        println(io, " of the partition on clock ", p.clock)
    end
    println(io, "  unknowns: ", p.unknowns)
    println(io, "  inputs:   ", p.inputs, " (", p.nu_user, " user)")
    println(io, "  outputs:  ", p.outputs, " (", p.ny_user, " user)")
    return nothing
end

"""
    $(TYPEDSIGNATURES)

The clock partitions that the operators of a boundary term connect. Returns the variable that
a boundary term such as `Sample(dt)(y)`, `Hold(u)` or `Sample(dt)(Hold(u))` is ultimately
applied to. For operators with several arguments, such as SynchToolkit's `Latest`, the first
argument is the signal that crosses the boundary.
"""
function _boundary_source(term::SymbolicT)
    src = arguments(term)[1]
    while iscall(src) && operation(src) isa Union{Sample, Hold}
        src = arguments(src)[1]
    end
    return src
end

"""
    $(TYPEDSIGNATURES)

Strip all `Shift`s from `var`. Discrete partitions store their variables shifted forward by
one step, so this recovers the variable of the model.
"""
function _strip_shifts(var::SymbolicT)
    while iscall(var) && operation(var) isa Shift
        var = arguments(var)[1]
    end
    return var
end

_is_boundary_operator(var::SymbolicT) = iscall(var) && operation(var) isa Operator &&
    !(operation(var) isa Union{Shift, Differential})

"""
    $(TYPEDSIGNATURES)

The array variable that the element `var` belongs to, or `var` itself if it is not an element
of an array variable.
"""
_array_root(var::SymbolicT) = iscall(var) && operation(var) === getindex ? arguments(var)[1] : var

"""
    $(TYPEDSIGNATURES)

Whether `var` occurs in no equation of the tearing state `ts`. This is the case for the
`Hold(w)` term of a nested clock change `Sample(clk)(Hold(w))`, which clock inference assigns
to the continuous partition although only the discrete partition reads it.
"""
function _is_dangling(ts::TearingState, var::SymbolicT)
    idx = findfirst(isequal(var), ts.fullvars)
    idx === nothing && return true
    return isempty(𝑑neighbors(ts.structure.graph, idx))
end

"""
    $(TYPEDSIGNATURES)

The set of model variables owned by the partition `ts`, that is, its variables with any
shifts removed and without boundary operator terms.
"""
function _partition_variables(ts::TearingState)
    vars = Set{SymbolicT}()
    for v in ts.fullvars
        _is_boundary_operator(v) && continue
        push!(vars, _strip_shifts(v))
    end
    return vars
end

"""
    $(TYPEDSIGNATURES)

The form in which the model variable `var` appears in the `fullvars` of the partition `ts`.
For discrete partitions this is the variable shifted forward by one step.
"""
function _fullvar_form(ts::TearingState, var::SymbolicT, iv::SymbolicT)
    shifted = MTKBase.simplify_shifts(Shift(iv, 1)(var))
    any(isequal(shifted), ts.fullvars) && return shifted
    any(isequal(var), ts.fullvars) && return var
    error("Variable $var was not found in the clock partition it was assigned to.")
end

"""
    $(TYPEDSIGNATURES)

Create the time-independent parameter that a user-specified input variable `var` is bound to.
The parameter has the same element type and shape as `var`. A time-independent parameter is
used because shifted discrete equations are shifted back to the current tick after
compilation, and a time-dependent parameter would be shifted along with them.
"""
function _input_parameter(var::SymbolicT, existing_names::Set{Symbol}, suffix::Symbol = :_input)
    name = Symbol(getname(var), suffix)
    while name in existing_names
        name = Symbol(name, :_)
    end
    push!(existing_names, name)
    return if symbolic_type(var) == SymbolicIndexingInterface.ArraySymbolic()
        T = eltype(symtype(var))
        unwrap(only(MTKBase.@parameters $name[SU.shape(var)...]::T))
    else
        T = symtype(var)
        unwrap(only(MTKBase.@parameters $name::T))
    end
end

"""
    $(TYPEDSIGNATURES)

Bind each variable in `inputs` to a new time-independent parameter by adding the equation
`input ~ parameter`. Returns the modified system and the parameters, in the order of `inputs`.
"""
function _bind_inputs_to_parameters(sys::System, inputs::Vector{SymbolicT})
    existing_names = Set{Symbol}(getname(p) for p in MTKBase.get_ps(sys))
    params = SymbolicT[]
    eqs = copy(MTKBase.get_eqs(sys))
    ics = copy(MTKBase.get_initial_conditions(sys))
    for var in inputs
        p = _input_parameter(var, existing_names)
        push!(params, p)
        push!(eqs, var ~ p)
        # A constant initial value of the input serves as the default operating-point value
        # of the parameter it is bound to.
        if haskey(ics, var) && _is_constant_value(ics[var])
            ics[p] = ics[var]
        end
    end
    @set! sys.eqs = eqs
    @set! sys.ps = [MTKBase.get_ps(sys); params]
    @set! sys.initial_conditions = ics
    return sys, params
end

_scalarized_vars(vars) = MTKBase.scalarized_vars(MTKBase.unwrap_vars(vars))

_is_constant_value(v) = v isa Union{Number, AbstractArray{<:Number}} || (v isa SymbolicT && SU.isconst(v))

"""
    $(TYPEDSIGNATURES)

Whether `clk` is a clock whose partition can be given a sample interval.
"""
_supported_clock(clk) = clk isa SciMLBase.ContinuousClock || (clk isa SciMLBase.PeriodicClock && clk.dt !== nothing)

_sample_interval(clk::SciMLBase.PeriodicClock) = clk.dt === nothing ? nothing : Float64(clk.dt)
_sample_interval(clk) = nothing

"""
    $(TYPEDEF)

Provider of operating-point values for [`linearize_hybrid`](@ref). Values are looked up in
the user-provided dictionary first, and then, when the operating point was given as a
[`LinearizationOpPoint`](@ref), evaluated from the solution at the requested time.
"""
struct HybridOperatingPoint{S}
    dict::Dict{SymbolicT, Any}
    """The symbols whose values were given explicitly by the user."""
    explicit::Vector{SymbolicT}
    sol::S
    t::Float64
end

function HybridOperatingPoint(op::AbstractDict, t)
    dict = Dict{SymbolicT, Any}()
    for (k, v) in op
        dict[unwrap(k)] = v
    end
    return HybridOperatingPoint(dict, collect(SymbolicT, keys(dict)), nothing, Float64(t))
end

function HybridOperatingPoint(op::LinearizationOpPoint, t)
    op.t isa AbstractVector && throw(
        ArgumentError(
            "`linearize_hybrid` does not support several time points in a `LinearizationOpPoint`."
        )
    )
    sol = op.sol
    t = Float64(op.t)
    sol_sys = MTKBase.indp_to_system(sol)
    dict = Dict{SymbolicT, Any}()
    u = sol(t)
    for (i, v) in enumerate(unknowns(sol_sys))
        dict[v] = u[i]
    end
    for p in parameters(sol_sys)
        if SymbolicIndexingInterface.is_timeseries_parameter(sol, p)
            val = _solution_value(sol, p, t)
            val === nothing && continue
            dict[p] = val
        else
            dict[p] = getp(sol, p)(sol)
        end
    end
    explicit = SymbolicT[]
    for (k, v) in op.op
        k = unwrap(k)
        dict[k] = v
        push!(explicit, k)
    end
    return HybridOperatingPoint(dict, explicit, sol, t)
end

# The value of `var` in `sol` at time `t`, or `nothing` when it is not available (a discrete
# variable before its first tick, or a variable the solution does not know).
function _solution_value(sol, var, t)
    try
        val = sol(t; idxs = var)
        return val
    catch
        return nothing
    end
end

"""
    $(TYPEDSIGNATURES)

The operating-point value of `var`, or `nothing` if none is available.
"""
function _op_value(op::HybridOperatingPoint, var::SymbolicT)
    haskey(op.dict, var) && return op.dict[var]
    op.sol === nothing && return nothing
    val = _solution_value(op.sol, var, op.t)
    val === nothing && return nothing
    op.dict[var] = val
    return val
end

"""
    $(TYPEDSIGNATURES)

Compile the tearing state `ts` of one clock partition. `boundary_inputs` are the operator
terms entering the partition, which become parameters, `outputs` are the model variables
that must remain accessible after compilation, and `input_params` are the parameters the
user-specified inputs of the partition are bound to. Returns the completed system and the
parameters that replace the boundary terms, in the order of `boundary_inputs`. Remaining
keyword arguments are forwarded to the structural simplification.
"""
function _compile_partition!(
        ts::TearingState, boundary_inputs::Vector{SymbolicT}, outputs::Vector{SymbolicT},
        input_params::Vector{SymbolicT}, discrete::Bool, iv::SymbolicT; kwargs...
    )
    if discrete
        outputs = [_fullvar_form(ts, o, iv) for o in outputs]
    else
        make_eqs_zero_equals!(ts)
    end
    # Boundary terms become parameters, as in `mtkcompile!` for the continuous partition of
    # a hybrid system. They are then replaced by plain parameters, since the compilation and
    # code generation of a discrete partition would otherwise have to treat every operator
    # that may cross a clock boundary, including operators defined outside ModelingToolkit,
    # as an atomic symbol.
    params = _boundary_parameters(ts.sys, boundary_inputs)
    if !isempty(boundary_inputs)
        inputs_to_parameters!(ts, OrderedSet{SymbolicT}(boundary_inputs), OrderedSet{SymbolicT}())
        rules = Dict{SymbolicT, SymbolicT}(zip(boundary_inputs, params))
        sys = ts.sys
        @set! sys.eqs = Equation[substitute(eq, rules) for eq in MTKBase.get_eqs(sys)]
        @set! sys.ps = SymbolicT[get(rules, p, p) for p in MTKBase.get_ps(sys)]
        binds = copy(parent(bindings(sys)))
        for (term, param) in rules
            haskey(binds, term) || continue
            binds[param] = binds[term]
            delete!(binds, term)
        end
        @set! sys.bindings = ROSymmapT(binds)
        @set! sys.initialization_eqs = Equation[
            substitute(eq, rules) for eq in MTKBase.get_initialization_eqs(sys)
        ]
        ts.sys = sys
        ts.original_eqs = Equation[substitute(eq, rules) for eq in ts.original_eqs]
    end
    # The boundary and input parameters are registered as the inputs of the partition. The
    # inline solution of linear subsystems derives the element type of its buffers from the
    # unknowns and the first input of the system, so the buffers take the type of the
    # perturbation when the Jacobians with respect to these parameters are evaluated.
    ssys = _mtkcompile!(
        ts; outputs = OrderedSet{SymbolicT}(outputs),
        input_parameters = OrderedSet{SymbolicT}([params; input_params]), kwargs...
    )
    # `complete` requires the registered inputs to follow the order of the parameters.
    input_set = Set{SymbolicT}([params; input_params])
    @set! ssys.inputs = OrderedSet{SymbolicT}(p for p in MTKBase.get_ps(ssys) if p in input_set)
    # The initialization equations of the model that refer to symbols of the partition only
    # are kept aside; those compatible with the operating point take part in the
    # initialization of the partition. Events are not accounted for, and the assertions of
    # the model may refer to variables of other partitions.
    init_eqs = _known_initialization_equations(ssys)
    @set! ssys.initialization_eqs = Equation[]
    @set! ssys.continuous_events = MTKBase.SymbolicContinuousCallback[]
    @set! ssys.discrete_events = MTKBase.SymbolicDiscreteCallback[]
    @set! ssys.assertions = Dict{SymbolicT, String}()
    # The parameters are inherited from the whole model. Parameters the partition does not use
    # would require operating-point values, and those bound to `missing` would be treated as
    # unknowns of the initialization problem.
    keep = SymbolicT[params; input_params]
    for eq in init_eqs
        union!(keep, _initialization_equation_symbols(eq))
    end
    ssys = _prune_parameters(ssys, keep)
    # Initial values are kept for the unknowns and parameters of the partition. Those of
    # observed variables would only constrain the initialization of the partition, and those
    # of variables of other partitions refer to symbols the partition does not know.
    known = Set{SymbolicT}()
    for v in Iterators.flatten((unknowns(ssys), parameters(ssys)))
        push!(known, _array_root(v))
    end
    @set! ssys.initial_conditions = filter(kv -> _array_root(first(kv)) in known, MTKBase.get_initial_conditions(ssys))
    for v in MTKBase.observables(ssys)
        push!(known, _array_root(v))
    end
    @set! ssys.guesses = filter(kv -> _array_root(first(kv)) in known, MTKBase.get_guesses(ssys))
    if discrete
        @set! ssys.is_discrete = true
    end
    return complete(ssys; split = true, allow_parameter_eqs = true), params, init_eqs
end

"""
    $(TYPEDSIGNATURES)

The symbols the initialization equation `eq` constrains: the variables and parameters
occurring in it, with the argument of a `Differential` or `Shift` term standing for the
term. Returns `nothing` if the equation contains an operator term that has no meaning in a
partition on its own, such as a boundary term of another partition.
"""
function _initialization_equation_symbols(eq)
    syms = Set{SymbolicT}()
    for v in get_variables(eq)
        v = unwrap(v)::SymbolicT
        if iscall(v) && operation(v) isa Operator
            op = operation(v)
            if op isa Union{Differential, Shift}
                push!(syms, _array_root(arguments(v)[1]))
            elseif op isa Union{Initial, MTKBase.Pre}
                push!(syms, v)
            else
                return nothing
            end
        else
            push!(syms, _array_root(v))
        end
    end
    return syms
end

"""
    $(TYPEDSIGNATURES)

The initialization equations of the compiled partition `ssys` whose symbols are all known to
the partition. Equations that refer to variables of other partitions are dropped.
"""
function _known_initialization_equations(ssys::System)
    known = Set{SymbolicT}()
    for v in Iterators.flatten((unknowns(ssys), MTKBase.observables(ssys), MTKBase.get_ps(ssys)))
        push!(known, _array_root(v))
    end
    iv = get_iv(ssys)
    iv === nothing || push!(known, iv)
    kept = Equation[]
    for eq in MTKBase.get_initialization_eqs(ssys)
        syms = _initialization_equation_symbols(eq)
        syms === nothing && continue
        all(in(known), syms) || continue
        push!(kept, eq)
    end
    return kept
end

"""
    $(TYPEDSIGNATURES)

Select, among the initialization equations `init_eqs` of the partition `ssys`, those the
initialization problem of the partition can use given its operating point `pop`: the
equations that determine at least one unknown or parameter bound to `missing` whose value
`pop` does not fix. Observed variables are expanded to the unknowns and parameters they
depend on for this purpose. Returns the selected equations and the set of symbols they
involve.
"""
function _usable_initialization_equations(
        ssys::System, init_eqs::Vector{Equation}, pop::Dict{SymbolicT, Any}
    )
    kept = Equation[]
    determined = Set{SymbolicT}()
    isempty(init_eqs) && return kept, determined
    fixed = Set{SymbolicT}(_array_root(k) for k in keys(pop))
    free = Set{SymbolicT}()
    for v in unknowns(ssys)
        r = _array_root(v)
        r in fixed || push!(free, r)
    end
    binds = bindings(ssys)
    for p in parameters(ssys)
        get(binds, p, nothing) === COMMON_MISSING || continue
        r = _array_root(p)
        r in fixed || push!(free, r)
    end
    rules = Dict{SymbolicT, SymbolicT}(eq.lhs => eq.rhs for eq in MTKBase.observed(ssys))
    maxiters = length(rules) + 1
    for eq in init_eqs
        expanded = fixpoint_sub(eq.lhs, rules; maxiters) ~ fixpoint_sub(eq.rhs, rules; maxiters)
        syms = _initialization_equation_symbols(expanded)
        syms === nothing && continue
        any(in(free), syms) || continue
        push!(kept, eq)
        union!(determined, syms)
    end
    return kept, determined
end

"""
    $(TYPEDSIGNATURES)

Values given by the initialization equations `init_eqs` of a discrete partition, keyed by the
symbol they determine. An equation `Shift(t, -1)(x) ~ expr` provides the value of the history
variable of `x`, an equation `x ~ expr` provides the value of `x` itself, which the history
variable takes as well under the assumption of an operating point that is stationary across
ticks.
"""
function _initialization_values(init_eqs::Vector{Equation})
    values = Dict{SymbolicT, Any}()
    for eq in reverse(init_eqs)
        lhs = eq.lhs
        if iscall(lhs) && operation(lhs) isa Shift
            values[MTKBase.default_toterm(lhs)] = eq.rhs
        else
            values[lhs] = eq.rhs
        end
    end
    return values
end

"""
    $(TYPEDSIGNATURES)

Give a value to the parameters of `ssys` bound to `missing` that neither `pop` fixes nor one
of the equations that involve the symbols in `determined` determines. The value is taken from
`op`, then from `fallback`, then from the guess of the parameter, and is zero otherwise, in
which case the parameter is recorded in `missing_vars`.
"""
function _undetermined_parameters!(
        pop::Dict{SymbolicT, Any}, ssys::System, op::HybridOperatingPoint,
        determined::Set{SymbolicT}, fallback::Dict{SymbolicT, Any}, missing_vars::Vector{SymbolicT}
    )
    binds = bindings(ssys)
    gs = MTKBase.guesses(ssys)
    for p in parameters(ssys)
        get(binds, p, nothing) === COMMON_MISSING || continue
        haskey(pop, p) && continue
        _array_root(p) in determined && continue
        val = _op_value(op, p)
        val === nothing && (val = get(fallback, p, nothing))
        val === nothing && (val = get(gs, p, nothing))
        if val === nothing
            push!(missing_vars, p)
            val = _zero_value(p)
        end
        pop[p] = val
    end
    return pop
end

"""
    $(TYPEDSIGNATURES)

Remove the parameters of `sys` that occur neither in its equations, observed equations,
initial values and guesses of its unknowns, nor in `keep`, nor in the bindings and initial
values of the parameters that are kept.
"""
function _prune_parameters(sys::System, keep::Vector{SymbolicT})
    used = Set{SymbolicT}()
    function add_variables!(expr)
        expr isa Union{SymbolicT, Equation} || return nothing
        for v in get_variables(expr)
            push!(used, _array_root(v))
        end
        return nothing
    end
    foreach(add_variables!, equations(sys))
    foreach(add_variables!, MTKBase.observed(sys))
    for p in keep
        push!(used, _array_root(p))
    end
    unknown_roots = Set{SymbolicT}(_array_root(v) for v in unknowns(sys))
    binds = bindings(sys)
    ics = MTKBase.get_initial_conditions(sys)
    gs = MTKBase.get_guesses(sys)
    for dict in (ics, gs), (k, v) in dict
        _array_root(k) in unknown_roots && add_variables!(v)
    end
    # Bindings and initial values of kept parameters may refer to further parameters.
    ps = MTKBase.get_ps(sys)
    changed = true
    while changed
        changed = false
        n = length(used)
        for p in ps
            _array_root(p) in used || continue
            for dict in (binds, ics, gs)
                add_variables!(get(dict, p, nothing))
            end
        end
        changed = length(used) != n
    end
    kept = SymbolicT[p for p in ps if _array_root(p) in used]
    length(kept) == length(ps) && return sys
    @set! sys.ps = kept
    newbinds = copy(parent(binds))
    filter!(kv -> _array_root(first(kv)) in used, newbinds)
    @set! sys.bindings = ROSymmapT(newbinds)
    return sys
end

"""
    $(TYPEDSIGNATURES)

The plain parameters that replace the boundary terms `terms` in a partition of `sys`. The
parameter of `Hold(u)` is named `u_hold`, that of `Sample(dt)(y)` `y_sample`, and so on; a
`SampleTime()` term becomes `sample_time`.
"""
function _boundary_parameters(sys::System, terms::Vector{SymbolicT})
    existing_names = Set{Symbol}(getname(p) for p in MTKBase.get_ps(sys))
    params = SymbolicT[]
    for term in terms
        op = operation(term)
        suffix = Symbol(:_, lowercase(String(nameof(typeof(op)))))
        if isempty(arguments(term))
            name = Symbol(lowercase(String(nameof(typeof(op)))))
            while name in existing_names
                name = Symbol(name, :_)
            end
            push!(existing_names, name)
            push!(params, unwrap(only(MTKBase.@parameters $name::Real)))
        else
            push!(params, _input_parameter(_boundary_source(term), existing_names, suffix))
        end
    end
    return params
end

"""
    $(TYPEDSIGNATURES)

Whether `k` is a symbol whose operating-point value can be set in the compiled system
`ssys`.
"""
function _settable_in(ssys::System, k::SymbolicT)
    return SymbolicIndexingInterface.is_variable(ssys, k) || is_parameter(ssys, k) || is_parameter(ssys, Initial(k))
end

"""
    $(TYPEDSIGNATURES)

The observed variables of `ssys` whose value is an expression of parameters only, such as an
alias of a boundary parameter. An operating-point value for such a variable does not determine
an unknown and would only add a constraint between parameters to the initialization problem.
"""
function _parameter_only_observables(ssys::System)
    obs = MTKBase.observed(ssys)
    result = Set{SymbolicT}()
    isempty(obs) && return result
    rules = Dict{SymbolicT, SymbolicT}(eq.lhs => eq.rhs for eq in obs)
    unknown_roots = Set{SymbolicT}(_array_root(v) for v in unknowns(ssys))
    for eq in obs
        rhs = fixpoint_sub(eq.rhs, rules; maxiters = length(rules) + 1)
        any(v -> _array_root(v) in unknown_roots, get_variables(rhs)) && continue
        push!(result, eq.lhs)
    end
    return result
end

_zero_value(var::SymbolicT) = SU.is_array_shape(SU.shape(var)) ? zeros(size(var)) : 0.0

"""
    $(TYPEDSIGNATURES)

Build the operating point of one partition. `boundary_values` maps the boundary parameters of
the partition to the model variable whose value they take, `input_params` maps the parameters
that user inputs are bound to, to the input variable. History variables of a discrete
partition take the value of the variable they are the history of, that is, the operating
point is assumed to be stationary across ticks. Returns the operating point and the boundary
parameters for which `op` provides no value, paired with the variable they take their value
from; these are set to zero for now and resolved later from the partition that variable
belongs to. Other symbols without a value are set to zero and recorded in `missing_vars`.
"""
function _partition_operating_point(
        ssys::System, op::HybridOperatingPoint, boundary_values::Vector{Pair{SymbolicT, SymbolicT}},
        input_params::Vector{Pair{SymbolicT, SymbolicT}}, init_values::Dict{SymbolicT, Any},
        discrete::Bool, missing_vars::Vector{SymbolicT}
    )
    result = Dict{SymbolicT, Any}()
    unresolved = Pair{SymbolicT, SymbolicT}[]
    input_vars = Set{SymbolicT}(var for (_, var) in input_params)
    parameter_only = _parameter_only_observables(ssys)
    # Values given explicitly by the user may address observed variables, which then
    # constrain the initialization of the partition.
    for k in op.explicit
        # User inputs are bound to parameters; their values are set through those.
        k in input_vars && continue
        k in parameter_only && continue
        _settable_in(ssys, k) || continue
        result[k] = op.dict[k]
    end
    # A solution provides a value for every variable it stores. Only the unknowns and the
    # parameters of the partition are taken from it; values of observed variables would
    # repeat the constraints the unknowns already satisfy.
    if op.sol !== nothing
        if !discrete
            for v in unknowns(ssys)
                haskey(result, v) && continue
                val = _op_value(op, v)
                val === nothing || (result[v] = val)
            end
        end
        for p in parameters(ssys)
            haskey(result, p) && continue
            iscall(p) && operation(p) isa Operator && continue
            val = _op_value(op, p)
            val === nothing || (result[p] = val)
        end
    end
    ics = initial_conditions(ssys)
    for (param, var) in boundary_values
        haskey(result, param) && continue
        val = _op_value(op, var)
        if val === nothing
            push!(unresolved, param => var)
            val = _zero_value(var)
        end
        result[param] = val
    end
    for (param, var) in input_params
        # Parameters of inputs of other partitions are pruned when the partition does not use them.
        is_parameter(ssys, param) || continue
        haskey(result, param) && continue
        val = _op_value(op, var)
        val === nothing && (val = get(ics, param, nothing))
        if val === nothing
            # The perturbation introduced by an analysis point is zero at the operating point.
            SU.hasmetadata(var, MTKBase.AnalysisVariable) || push!(missing_vars, var)
            val = _zero_value(var)
        end
        result[param] = val
    end
    if discrete
        for v in unknowns(ssys)
            haskey(result, v) && continue
            val = _history_value(op, ics, init_values, v)
            if val === nothing
                push!(missing_vars, v)
                val = 0.0
            end
            result[v] = val
        end
    end
    return result, unresolved
end

"""
    $(TYPEDSIGNATURES)

The operating-point value of the history variable `v` of a discrete partition, that is, the
value of the variable `v` is the history of, taken from `op`, from the initial values `ics`
or from the values `init_values` given by initialization equations. For an element of an
array variable, the value of the whole array is consulted as well. Returns `nothing` if no
value is available.
"""
function _history_value(op::HybridOperatingPoint, ics, init_values, v::SymbolicT)
    root = _array_root(v)
    base_root = MTKBase.getunshifted(root)
    base_root === nothing && (base_root = root)
    function lookup(sym)
        val = _op_value(op, sym)
        val === nothing && (val = get(ics, sym, nothing))
        val === nothing && (val = get(init_values, sym, nothing))
        return val
    end
    if root === v
        val = lookup(base_root)
        val === nothing && (val = get(ics, v, nothing))
        val === nothing && (val = get(init_values, v, nothing))
        return val
    end
    idxs = Int[Int(SU.unwrap_const(i)) for i in arguments(v)[2:end]]
    base = unwrap(wrap(base_root)[idxs...])
    val = lookup(base)
    val === nothing && (val = get(ics, v, nothing))
    val === nothing && (val = get(init_values, v, nothing))
    val === nothing || return val
    arr = lookup(base_root)
    arr === nothing && (arr = get(init_values, root, nothing))
    arr isa AbstractArray && return arr[idxs...]
    arr isa SymbolicT && SU.is_array_shape(SU.shape(arr)) && return unwrap(wrap(arr)[idxs...])
    return nothing
end

"""
    $(TYPEDSIGNATURES)

Set the boundary parameters listed in `unresolved` to the values of their source variables
recorded in `source_values`, keyed by the partition index of the source variable (given by
`partition_of`) and the variable. Entries whose value is not available yet are kept.
"""
function _resolve_boundary_values!(
        pop::Dict{SymbolicT, Any}, unresolved::Vector{Pair{SymbolicT, SymbolicT}},
        source_values::Dict{Tuple{Int, SymbolicT}, Any}, partition_of
    )
    filter!(unresolved) do (param, source)
        key = (partition_of(source), source)
        haskey(source_values, key) || return true
        pop[param] = source_values[key]
        return false
    end
    return pop
end

"""
    $(TYPEDSIGNATURES)

Warn if the state `x0` at which the continuous partition `ssys` was linearized differs from
the values its operating point `pop` requested for the differential unknowns, whose indices
are `diff_idxs`. This happens when the operating point also fixes signals the model
determines from the unknowns, so that the initialization problem is overdetermined and solved
in the least-squares sense. Requested values of algebraic unknowns serve as guesses and are
not compared.
"""
function _warn_operating_point_mismatch(
        ssys::System, pop::Dict{SymbolicT, Any}, x0::Vector{Float64}, diff_idxs::Vector{Int}, description
    )
    mismatched = String[]
    vars = unknowns(ssys)
    for i in diff_idxs
        v = vars[i]
        val = get(pop, v, nothing)
        val isa Number || continue
        isapprox(x0[i], val; rtol = 1.0e-6, atol = 1.0e-8) && continue
        push!(mismatched, "$v: requested $val, obtained $(x0[i])")
    end
    isempty(mismatched) && return nothing
    @warn "The initialization of the $description did not retain the requested values of $(join(mismatched, "; ")). This happens when the operating point also fixes signals the model determines from these unknowns. The linearization is taken at the obtained values; see the `x0` field of the partition."
    return nothing
end

"""
    $(TYPEDSIGNATURES)

Raise an error if the compiled partition `ssys`, its input parameters `inputs` or its
outputs `outputs` contain Boolean-valued variables. Linearization is not defined with respect
to a Boolean variable, and the generated functions of the partition cannot be evaluated with
real-valued perturbations of it. `description` names the partition in the message.
"""
function _check_boolean_variables(ssys::System, inputs::Vector{SymbolicT}, outputs::Vector{SymbolicT}, description)
    boolean = SymbolicT[]
    for v in Iterators.flatten((unknowns(ssys), inputs, outputs))
        T = symtype(v)
        T <: AbstractArray && (T = eltype(T))
        T <: Bool && push!(boolean, v)
    end
    isempty(boolean) && return nothing
    return error(
        """
        The $description has the Boolean variables $(join(string.(boolean), ", ")). \
        Linearization is not defined with respect to Boolean variables.
        """
    )
end

"""
    $(TYPEDSIGNATURES)

Whether `err` indicates a failure of automatic differentiation, such as a method that does not
accept dual numbers, as opposed to a structural or operating-point error.
"""
function _is_autodiff_error(err)
    err isa MethodError && return true
    msg = sprint(showerror, err)
    return occursin("Dual", msg) || occursin("ForwardDiff", msg)
end

"""
    $(TYPEDSIGNATURES)

Evaluate `f` with `autodiff`. If it throws an error that indicates a failure of automatic
differentiation and `fallback` differs from `autodiff`, emit a warning and evaluate `f` again
with `fallback`. `description` names the partition in the warning.
"""
function _with_autodiff_fallback(f, autodiff, fallback, description)
    fallback === nothing && return f(autodiff)
    typeof(fallback) === typeof(autodiff) && return f(autodiff)
    return try
        f(autodiff)
    catch err
        _is_autodiff_error(err) || rethrow()
        @warn "Linearization of the $description with $autodiff failed, retrying with $fallback." exception = (err, catch_backtrace())
        f(fallback)
    end
end

"""
    $(TYPEDSIGNATURES)

The values of `outputs` of the compiled partition `ssys` at the state `u` and parameters `p`.
"""
function _evaluate_outputs(ssys::System, outputs::Vector{SymbolicT}, u, p, t; eval_expression, eval_module)
    isempty(outputs) && return Float64[]
    h = build_explicit_observed_function(
        ssys, outputs, GeneratedFunctionOptions(; expression = Val{false}, eval_expression, eval_module)
    )
    return collect(h(u, p, t))
end

"""
    $(TYPEDSIGNATURES)

Linearize the compiled continuous partition `ssys` from the parameters `inputs` at the
operating point `op` with the prepared `lin_fun`.
"""
function _linearize_continuous_partition(
        ssys::System, lin_fun, inputs::Vector{SymbolicT}, op::Dict, t; allow_input_derivatives
    )
    # The Jacobians with respect to the inputs are evaluated at the input values stored in
    # the problem of `lin_fun`, which was built before the boundary values were resolved.
    prob = lin_fun.prob
    if !isempty(inputs)
        vals = map(inputs) do k
            v = op[k]
            v isa Union{Number, AbstractArray{<:Number}} ? v : getu(prob, unwrap(v))(prob)
        end
        setp(prob, inputs)(prob, vals)
    end
    mats, extras = linearize(ssys, lin_fun; op, t, allow_input_derivatives)
    x0 = extras.x === nothing ? Float64[] : collect(Float64, extras.x)
    return mats, x0
end

"""
    $(TYPEDSIGNATURES)

The update function, history state and parameters of the compiled discrete partition `ssys`
at the operating point `op`.
"""
function _discrete_partition_data(ssys::System, op::Dict, t; eval_expression, eval_module)
    if MTKBase.has_alg_equations(ssys)
        error(
            """
            The clock partition with unknowns $(unknowns(ssys)) contains algebraic equations \
            after structural simplification, which indicates an algebraic loop within a single \
            clock partition. Linearization of implicit discrete partitions is not supported.
            """
        )
    end
    return MTKBase.process_SciMLProblem(
        SciMLBase.DiscreteFunction{true}, ssys, op; build_initializeprob = false, t,
        eval_expression, eval_module
    )
end

"""
    $(TYPEDSIGNATURES)

Linearize the compiled discrete partition `ssys` from the parameters `inputs` to the
variables `outputs` at the operating point `op`. The state consists of the history
variables of the partition, the generated update function maps them to the values at the
current tick, and the observed function evaluates `outputs` at the current tick.
"""
function _linearize_discrete_partition(
        ssys::System, inputs::Vector{SymbolicT}, outputs::Vector{SymbolicT}, op::Dict, t;
        autodiff, eval_expression = false, eval_module = @__MODULE__
    )
    f, u0, p = _discrete_partition_data(ssys, op, t; eval_expression, eval_module)
    h = build_explicit_observed_function(
        ssys, outputs, GeneratedFunctionOptions(; expression = Val{false}, eval_expression, eval_module)
    )
    setter = setp_oop(ssys, inputs)
    input_vals = collect(Float64, getp(ssys, inputs)(p))
    ny = length(outputs)
    nu = length(inputs)
    # The Jacobians with respect to the inputs are only evaluated when there are inputs;
    # finite differences reject an empty vector of differentiation variables.
    if u0 === nothing || isempty(u0)
        A = zeros(0, 0)
        B = zeros(0, nu)
        C = zeros(ny, 0)
        D = if nu == 0
            zeros(ny, 0)
        else
            DI.jacobian(inp -> vec(collect(h(nothing, setter(p, inp), t))), autodiff, input_vals)
        end
        x0 = Float64[]
    else
        x0 = collect(Float64, u0)
        A = DI.jacobian((du, u) -> f(du, u, p, t), similar(x0), autodiff, x0)
        C = DI.jacobian(u -> vec(collect(h(u, p, t))), autodiff, x0)
        if nu == 0
            B = zeros(length(x0), 0)
            D = zeros(ny, 0)
        else
            B = DI.jacobian((du, inp) -> f(du, x0, setter(p, inp), t), similar(x0), autodiff, input_vals)
            D = DI.jacobian(inp -> vec(collect(h(x0, setter(p, inp), t))), autodiff, input_vals)
        end
    end
    return (; A = Matrix{Float64}(A), B = Matrix{Float64}(B), C = Matrix{Float64}(C), D = Matrix{Float64}(D)), x0
end

"""
    hl = linearize_hybrid(sys, inputs, outputs; op = Dict(), t = 0.0, kwargs...)

Linearize the hybrid (clocked) system `sys` partition by partition. The model may contain
periodic clock partitions (synchronous subsystems written with `ShiftIndex`, `Sample`,
`Hold` and components built on them) alongside the continuous equations. Each clock
partition and the continuous partition are linearized separately, and every signal that
crosses a clock boundary becomes an input of the partition it enters and an output of the
partition it leaves, in addition to the user-specified `inputs` and `outputs`. The result is
a [`HybridLinearization`](@ref) holding one [`ClockPartitionLinearization`](@ref) per
partition and the connections between them, from which the sampled-data loop can be
assembled with the advanced interface of `ControlSystemsBase.feedback` or with
`RobustAndOptimalControl.connect`.

The Jacobians are computed from the generated functions of each compiled partition with
`autodiff`, so models that call arbitrary Julia functions are supported. If differentiation
of a partition fails, it is retried with `fallback_autodiff`.

Only periodic clocks are supported. Partitions on other clocks, as well as continuous and
discrete events, assertions and state machines of the model, are not accounted for and
produce a warning. Linearization is not defined with respect to Boolean-valued variables,
and a partition containing them produces an error. Non-finite entries in the matrices of a
partition, which arise when an equation is not differentiable at the operating point, are
reported by a warning.

# Arguments

- `sys`: The unsimplified system.
- `inputs`: The input variables of the linearization, or analysis points.
- `outputs`: The output variables of the linearization, or analysis points.

# Keyword Arguments

- `op`: The operating point, a dictionary of variable and parameter values or a
  [`LinearizationOpPoint`](@ref) wrapping a solution and a time. The values of the state
  variables of every partition are taken from `op`; the history variables of a discrete
  partition take the value of the variable they are the history of, that is, the operating
  point is assumed to be stationary across ticks, or the value an initialization equation of
  the model gives. The value of a signal crossing a clock boundary is taken from `op` if
  present and otherwise evaluated from the operating point of the partition the signal
  originates from, with the boundary signals entering that partition that are not in `op` set
  to zero. The initialization equations of the model that refer to symbols of the continuous
  partition only and determine a value `op` does not fix take part in the initialization of
  the continuous partition, so parameters bound to `missing` and constant states that such
  equations determine take their values from them. Remaining values that are not available
  default to the guess of the symbol or to zero. From a `LinearizationOpPoint`, the values
  of the unknowns and parameters of each partition are taken; the values the solution holds
  for observed variables are not, since they would repeat the constraints the unknowns
  already satisfy.
- `t`: The time at which to linearize. Ignored if `op` is a `LinearizationOpPoint`.
- `autodiff`: The `ADTypes.AbstractADType` used to compute the Jacobians. Defaults to
  `AutoForwardDiff()`.
- `fallback_autodiff`: The AD type used when `autodiff` fails for a partition. Defaults to
  `AutoFiniteDiff()`; pass `nothing` to disable the fallback.
- `allow_input_derivatives`: See [`linearize`](@ref). Applies to the continuous partition.
- `warn_missing_op`: Whether to warn about operating-point values that default to zero.
- `warn_unsupported`: Whether to warn about constructs that are not accounted for.
- `initialize`, `initializealg`, `initialization_abstol`, `initialization_reltol`,
  `initialization_solver_alg`, `guesses`, `missing_guess_value`, `eval_expression`,
  `eval_module`: See [`linearization_function`](@ref). They apply to the continuous partition.
- `loop_openings`, `system_modifier`: See [`get_sensitivity`](@ref). Only for analysis
  points.
- Remaining keyword arguments are forwarded to the structural simplification of each
  partition.

See also [`linearize`](@ref), [`HybridLinearization`](@ref) and
[`ClockPartitionLinearization`](@ref).
"""
function linearize_hybrid(
        sys::AbstractSystem, inputs, outputs;
        op = Dict{SymbolicT, SymbolicT}(), t = 0.0,
        autodiff = AutoForwardDiff(), fallback_autodiff = AutoFiniteDiff(),
        allow_input_derivatives = false, warn_missing_op = true, warn_unsupported = true,
        loop_opening_params = SymbolicT[],
        initialize = true, initializealg = nothing,
        initialization_abstol = 1.0e-5, initialization_reltol = 1.0e-3,
        initialization_solver_alg = nothing, guesses = Dict{SymbolicT, SymbolicT}(),
        missing_guess_value = MTKBase.default_missing_guess_value(),
        eval_expression = false, eval_module = @__MODULE__,
        kwargs...
    )
    inputs isa AbstractVector || (inputs = [inputs])
    outputs isa AbstractVector || (outputs = [outputs])
    inputs = _scalarized_vars(inputs)
    outputs = _scalarized_vars(outputs)
    hop = HybridOperatingPoint(op, t)
    t = hop.t
    loop_opening_params = collect(SymbolicT, loop_opening_params)
    _check_loop_opening_op(loop_opening_params, hop.dict)

    if warn_unsupported
        if !isempty(continuous_events(sys)) || !isempty(discrete_events(sys))
            @warn "The system has continuous or discrete events. Events are not accounted for by `linearize_hybrid`."
        end
        if !isempty(assertions(sys))
            @warn "The system has assertions. Assertions are not accounted for by `linearize_hybrid`."
        end
    end

    # The model is preprocessed as in `mtkcompile`; the inputs are then bound to parameters
    # before clock inference, so that they stay variables of their partition.
    sys, statemachines = extract_top_level_statemachines(sys)
    if warn_unsupported && !isempty(statemachines)
        @warn "The system has state machines. State machines are not accounted for by `linearize_hybrid`."
    end
    sys = expand_connections(sys)
    sys = MTKBase.discover_maybe_zeros(sys)
    sys = MTKBase.apply_limited_lowering(sys)
    sys, input_params = _bind_inputs_to_parameters(sys, inputs)
    iv = get_iv(sys)::SymbolicT
    state = TearingState(sys; defer_scalarization = true)
    ci = MTKTearing.ClockInference(state)
    ci = MTKTearing.infer_clocks!(ci)
    tss, clocked_inputs, continuous_id, id_to_clock = MTKTearing.split_system(ci)
    npart = length(tss)
    for i in 1:npart
        # `split_system` shares this buffer between partitions; tearing appends to it.
        @set! tss[i].additional_observed = copy(tss[i].additional_observed)
        # The incidence structure of a partition is only complete once its array equations
        # have been scalarized.
        MTKTearing.scalarize_tearing_state_eqs!(tss[i])
    end
    # `split_system` also lists shifted history variables as time-domain conversions; only
    # operator terms other than shifts cross a clock boundary.
    boundary_inputs = [filter(_is_boundary_operator, ins) for ins in clocked_inputs]

    if warn_unsupported
        for clk in id_to_clock
            _supported_clock(clk) && continue
            @warn "The clock $clk is not a periodic clock. Its partition is linearized as a discrete map without a sample interval."
        end
    end

    # Assignment of model variables to partitions.
    partition_vars = map(_partition_variables, tss)
    function partition_of(var)
        for (i, vars) in enumerate(partition_vars)
            var in vars && return i
        end
        throw(IONotFoundError("inputs or outputs", nameof(sys), [var]))
    end

    # User I/O per partition, in the order given by the user.
    user_inputs = [SymbolicT[] for _ in 1:npart]
    user_input_params = [SymbolicT[] for _ in 1:npart]
    user_outputs = [SymbolicT[] for _ in 1:npart]
    for (var, param) in zip(inputs, input_params)
        i = partition_of(var)
        push!(user_inputs[i], var)
        push!(user_input_params[i], param)
    end
    for var in outputs
        push!(user_outputs[partition_of(var)], var)
    end

    # Boundary analysis: each boundary term entering partition `i` is driven by a variable of
    # another partition `j`. Terms that no equation of their partition reads are compiled as
    # parameters but not exposed. `SampleTime()` terms are parameters of their partition whose
    # value is the sample interval of its clock.
    struct_conn = @NamedTuple{from::Int, source::SymbolicT, to::Int, term::SymbolicT}
    connections = struct_conn[]
    sample_time_terms = [SymbolicT[] for _ in 1:npart]
    for i in 1:npart
        for term in boundary_inputs[i]
            if operation(term) isa SampleTime
                _sample_interval(id_to_clock[i]) === nothing && error(
                    "`SampleTime()` is used on the clock $(id_to_clock[i]), which has no fixed sample interval."
                )
                push!(sample_time_terms[i], term)
                continue
            end
            _is_dangling(tss[i], term) && continue
            source = _boundary_source(term)
            j = findfirst(vars -> source in vars, partition_vars)
            if j === nothing
                error("The signal $source crossing a clock boundary through $term could not be assigned to a clock partition.")
            end
            push!(connections, (; from = j, source, to = i, term))
        end
    end
    # Group boundary inputs by the partition they come from, and boundary outputs by the
    # partition they go to. An output driving several partitions appears once, and an output
    # that is also a user output is not repeated.
    boundary_in = [SymbolicT[] for _ in 1:npart]
    boundary_out = [SymbolicT[] for _ in 1:npart]
    for i in 1:npart
        for j in 1:npart
            for c in connections
                (c.to == i && c.from == j) || continue
                push!(boundary_in[i], c.term)
            end
            for c in connections
                (c.from == i && c.to == j) || continue
                any(isequal(c.source), boundary_out[i]) && continue
                any(isequal(c.source), user_outputs[i]) && continue
                push!(boundary_out[i], c.source)
            end
        end
    end

    # Compile every partition.
    partitions_data = map(1:npart) do i
        discrete = i != continuous_id
        outs = [user_outputs[i]; boundary_out[i]]
        ssys, boundary_params, init_eqs = _compile_partition!(
            tss[i], boundary_inputs[i], outs, user_input_params[i], discrete, iv; kwargs...
        )
        term_to_param = Dict{SymbolicT, SymbolicT}(zip(boundary_inputs[i], boundary_params))
        ins = [user_input_params[i]; SymbolicT[term_to_param[term] for term in boundary_in[i]]]
        boundary_values = Pair{SymbolicT, SymbolicT}[
            term_to_param[term] => _boundary_source(term)
                for term in boundary_inputs[i] if !(operation(term) isa SampleTime)
        ]
        sample_time_values = Pair{SymbolicT, Float64}[
            term_to_param[term] => _sample_interval(id_to_clock[i]) for term in sample_time_terms[i]
        ]
        (; ssys, discrete, outs, ins, boundary_values, sample_time_values, init_eqs, clock = id_to_clock[i])
    end
    for pd in partitions_data
        description = pd.discrete ? "clock partition on $(pd.clock)" : "continuous partition"
        _check_boolean_variables(pd.ssys, pd.ins, pd.outs, description)
    end

    # Operating points of the partitions. Boundary signals without a value in `op` are resolved
    # from the partition they originate from, in the order continuous partition first.
    missing_vars = SymbolicT[]
    all_input_params = Pair{SymbolicT, SymbolicT}[param => var for (var, param) in zip(inputs, input_params)]
    pops = Vector{Dict{SymbolicT, Any}}(undef, npart)
    unresolved = Vector{Vector{Pair{SymbolicT, SymbolicT}}}(undef, npart)
    # The initialization equations of the model that the initialization problem of the
    # continuous partition uses, such as those determining parameters bound to `missing`.
    init_eqs = Vector{Vector{Equation}}(undef, npart)
    for i in 1:npart
        pd = partitions_data[i]
        init_values = pd.discrete ? _initialization_values(pd.init_eqs) : Dict{SymbolicT, Any}()
        pops[i], unresolved[i] = _partition_operating_point(
            pd.ssys, hop, pd.boundary_values, all_input_params, init_values, pd.discrete, missing_vars
        )
        for (param, val) in pd.sample_time_values
            pops[i][param] = val
        end
        if pd.discrete
            # A discrete partition has no initialization problem.
            init_eqs[i] = Equation[]
            determined = Set{SymbolicT}()
        else
            init_eqs[i], determined = _usable_initialization_equations(pd.ssys, pd.init_eqs, pops[i])
        end
        _undetermined_parameters!(pops[i], pd.ssys, hop, determined, init_values, missing_vars)
    end
    order = continuous_id == 0 ? collect(1:npart) : [continuous_id; filter(!=(continuous_id), 1:npart)]
    codegen = (; eval_expression, eval_module)
    linearization_kwargs = (;
        initialize, initializealg, initialization_abstol, initialization_reltol,
        initialization_solver_alg, guesses, missing_guess_value, eval_expression, eval_module,
        loop_opening_params,
    )
    source_values = Dict{Tuple{Int, SymbolicT}, Any}()
    lin_funs = Vector{Any}(nothing, npart)
    for i in order
        pd = partitions_data[i]
        _resolve_boundary_values!(pops[i], unresolved[i], source_values, partition_of)
        if pd.discrete
            _, u0, p = _discrete_partition_data(pd.ssys, pops[i], t; codegen...)
            values = _evaluate_outputs(pd.ssys, pd.outs, u0, p, t; codegen...)
        else
            lin_fun, _ = _linearization_function_compiled(
                pd.ssys, pd.ins, pd.outs; op = pops[i], t, autodiff,
                initialization_eqs = init_eqs[i], linearization_kwargs...
            )
            lin_funs[i] = lin_fun
            prob = lin_fun.prob
            values = _evaluate_outputs(pd.ssys, pd.outs, prob.u0, prob.p, t; codegen...)
        end
        for (k, source) in enumerate(pd.outs)
            source_values[(i, source)] = values[k]
        end
    end
    for i in 1:npart
        _resolve_boundary_values!(pops[i], unresolved[i], source_values, partition_of)
    end

    # Linearize every partition.
    results = Vector{ClockPartitionLinearization}(undef, npart)
    for i in 1:npart
        pd = partitions_data[i]
        description = pd.discrete ? "clock partition on $(pd.clock)" : "continuous partition"
        mats, x0 = _with_autodiff_fallback(autodiff, fallback_autodiff, description) do ad
            if pd.discrete
                _linearize_discrete_partition(pd.ssys, pd.ins, pd.outs, pops[i], t; autodiff = ad, codegen...)
            else
                lin_fun = if ad === autodiff
                    lin_funs[i]
                else
                    first(
                        _linearization_function_compiled(
                            pd.ssys, pd.ins, pd.outs; op = pops[i], t, autodiff = ad,
                            initialization_eqs = init_eqs[i], linearization_kwargs...
                        )
                    )
                end
                _linearize_continuous_partition(pd.ssys, lin_fun, pd.ins, pops[i], t; allow_input_derivatives)
            end
        end
        if !all(m -> all(isfinite, m), (mats.A, mats.B, mats.C, mats.D))
            @warn "The linearization of the $description has non-finite entries. The equations of the partition are not differentiable at the operating point."
        end
        pd.discrete || _warn_operating_point_mismatch(pd.ssys, pops[i], x0, lin_funs[i].diff_idxs, description)
        results[i] = ClockPartitionLinearization(
            mats.A, mats.B, mats.C, mats.D, unknowns(pd.ssys),
            [user_inputs[i]; boundary_in[i]], pd.outs, length(user_inputs[i]), length(user_outputs[i]),
            pd.clock, _sample_interval(pd.clock), x0, pd.ssys
        )
    end
    if warn_missing_op && !isempty(missing_vars)
        unique!(missing_vars)
        @warn "No operating-point value was found for the following variables; zero is used: $(join(string.(missing_vars), ", ")). Pass values in `op` or disable this warning with `warn_missing_op = false`."
    end

    # Order the partitions with the continuous one first and translate the connections to
    # indices into the input and output lists.
    position = Dict(orig => new for (new, orig) in enumerate(order))
    partitions = results[order]
    index_connections = @NamedTuple{from::Int, output::Int, to::Int, input::Int}[]
    for c in connections
        from = position[c.from]
        to = position[c.to]
        output = findfirst(isequal(c.source), partitions[from].outputs)::Int
        input = findfirst(isequal(c.term), partitions[to].inputs)::Int
        push!(index_connections, (; from, output, to, input))
    end
    sort!(index_connections; by = c -> (c.to, c.input))
    continuous_index = continuous_id == 0 ? nothing : position[continuous_id]
    return HybridLinearization(partitions, index_connections, continuous_index)
end
