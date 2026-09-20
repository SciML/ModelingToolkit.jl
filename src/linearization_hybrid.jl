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
            io, ": ", size(p.A, 1), " state", size(p.A, 1) == 1 ? "" : "s", ", ",
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
        # An initial value of the input serves as the default operating-point value of the
        # parameter it is bound to.
        if haskey(ics, var)
            ics[p] = ics[var]
        end
    end
    @set! sys.eqs = eqs
    @set! sys.ps = [MTKBase.get_ps(sys); params]
    @set! sys.initial_conditions = ics
    return sys, params
end

_scalarized_vars(vars) = MTKBase.scalarized_vars(MTKBase.unwrap_vars(vars))

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
    sol::S
    t::Float64
end

function HybridOperatingPoint(op::AbstractDict, t)
    dict = Dict{SymbolicT, Any}()
    for (k, v) in op
        dict[unwrap(k)] = v
    end
    return HybridOperatingPoint(dict, nothing, Float64(t))
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
    for (k, v) in op.op
        dict[unwrap(k)] = v
    end
    return HybridOperatingPoint(dict, sol, t)
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
terms entering the partition, which become parameters, and `outputs` are the model variables
that must remain accessible after compilation. Remaining keyword arguments are forwarded to
the structural simplification.
"""
function _compile_partition!(
        ts::TearingState, boundary_inputs::Vector{SymbolicT}, outputs::Vector{SymbolicT},
        discrete::Bool, iv::SymbolicT; kwargs...
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
        ts.sys = sys
        ts.original_eqs = Equation[substitute(eq, rules) for eq in ts.original_eqs]
    end
    ssys = _mtkcompile!(ts; outputs = OrderedSet{SymbolicT}(outputs), kwargs...)
    # Initialization equations of the model, such as initial values of held signals, do not
    # concern the operating point of the linearization, and events are not accounted for.
    @set! ssys.initialization_eqs = Equation[]
    @set! ssys.continuous_events = MTKBase.SymbolicContinuousCallback[]
    @set! ssys.discrete_events = MTKBase.SymbolicDiscreteCallback[]
    # Initial values and guesses of variables of other partitions are inherited from the
    # model and would be treated as initial conditions of unknown symbols.
    known = Set{SymbolicT}()
    for v in Iterators.flatten((unknowns(ssys), MTKBase.observables(ssys), parameters(ssys)))
        push!(known, first(MTKBase.split_indexed_var(v)))
    end
    isknown(kv) = first(MTKBase.split_indexed_var(first(kv))) in known
    @set! ssys.initial_conditions = filter(isknown, MTKBase.get_initial_conditions(ssys))
    @set! ssys.guesses = filter(isknown, MTKBase.get_guesses(ssys))
    if discrete
        @set! ssys.is_discrete = true
    end
    return complete(ssys; split = true, allow_parameter_eqs = true), params
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

Build the operating point of one partition. `boundary_values` maps the boundary parameters of
the partition to the model variable whose value they take, `input_params` maps the parameters
that user inputs are bound to, to the input variable. History variables of a discrete
partition take the value of the variable they are the history of, that is, the operating point
is assumed to be stationary across ticks. Variables without a value are set to zero and
recorded in `missing_vars`.
"""
function _partition_operating_point(
        ssys::System, op::HybridOperatingPoint, boundary_values::Vector{Pair{SymbolicT, SymbolicT}},
        input_params::Vector{Pair{SymbolicT, SymbolicT}}, discrete::Bool, missing_vars::Vector{SymbolicT}
    )
    result = Dict{SymbolicT, Any}()
    input_vars = Set{SymbolicT}(var for (_, var) in input_params)
    for (k, v) in op.dict
        # User inputs are bound to parameters; their values are set through those.
        k in input_vars && continue
        _settable_in(ssys, k) || continue
        result[k] = v
    end
    ics = initial_conditions(ssys)
    function value_or_zero(var)
        val = _op_value(op, var)
        val === nothing && (val = get(ics, var, nothing))
        if val === nothing
            push!(missing_vars, var)
            val = SU.is_array_shape(SU.shape(var)) ? zeros(size(var)) : 0.0
        end
        return val
    end
    for (param, var) in boundary_values
        haskey(result, param) && continue
        result[param] = value_or_zero(var)
    end
    for (param, var) in input_params
        haskey(result, param) && continue
        result[param] = value_or_zero(var)
    end
    if discrete
        for v in unknowns(ssys)
            haskey(result, v) && continue
            base = MTKBase.getunshifted(v)
            base === nothing && (base = v)
            val = _op_value(op, base)
            if val === nothing
                val = get(ics, v, nothing)
                val === nothing && (val = get(ics, base, nothing))
            end
            if val === nothing
                push!(missing_vars, v)
                val = 0.0
            end
            result[v] = val
        end
    end
    return result
end

"""
    $(TYPEDSIGNATURES)

Evaluate `f`, and if it throws, evaluate it again with `fallback` in place of `autodiff`
after emitting a warning. `description` names the partition in the warning.
"""
function _with_autodiff_fallback(f, autodiff, fallback, description)
    fallback === nothing && return f(autodiff)
    typeof(fallback) === typeof(autodiff) && return f(autodiff)
    return try
        f(autodiff)
    catch err
        @warn "Linearization of the $description with $autodiff failed, retrying with $fallback." exception = (err, catch_backtrace())
        f(fallback)
    end
end

"""
    $(TYPEDSIGNATURES)

Linearize the compiled continuous partition `ssys` from the parameters `inputs` to the
variables `outputs` at the operating point `op`.
"""
function _linearize_continuous_partition(
        ssys::System, inputs::Vector{SymbolicT}, outputs::Vector{SymbolicT}, op::Dict, t;
        autodiff, allow_input_derivatives, kwargs...
    )
    lin_fun, _ = _linearization_function_compiled(
        ssys, inputs, outputs; op, t, autodiff, kwargs...
    )
    mats, extras = linearize(ssys, lin_fun; op, t, allow_input_derivatives)
    x0 = extras.x === nothing ? Float64[] : collect(Float64, extras.x)
    return mats, x0
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
    if MTKBase.has_alg_equations(ssys)
        error(
            """
            The clock partition with unknowns $(unknowns(ssys)) contains algebraic equations \
            after structural simplification, which indicates an algebraic loop within a single \
            clock partition. Linearization of implicit discrete partitions is not supported.
            """
        )
    end
    f, u0, p = MTKBase.process_SciMLProblem(
        SciMLBase.DiscreteFunction{true}, ssys, op; build_initializeprob = false, t,
        eval_expression, eval_module
    )
    h = build_explicit_observed_function(
        ssys, outputs, GeneratedFunctionOptions(; expression = Val{false}, eval_expression, eval_module)
    )
    setter = setp_oop(ssys, inputs)
    input_vals = collect(Float64, getp(ssys, inputs)(p))
    ny = length(outputs)
    nu = length(inputs)
    if u0 === nothing || isempty(u0)
        A = zeros(0, 0)
        B = zeros(0, nu)
        C = zeros(ny, 0)
        D = DI.jacobian(inp -> vec(collect(h(nothing, setter(p, inp), t))), autodiff, input_vals)
        x0 = Float64[]
    else
        x0 = collect(Float64, u0)
        A = DI.jacobian((du, u) -> f(du, u, p, t), similar(x0), autodiff, x0)
        B = DI.jacobian((du, inp) -> f(du, x0, setter(p, inp), t), similar(x0), autodiff, input_vals)
        C = DI.jacobian(u -> vec(collect(h(u, p, t))), autodiff, x0)
        D = DI.jacobian(inp -> vec(collect(h(x0, setter(p, inp), t))), autodiff, input_vals)
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
discrete events of the model, are not accounted for and produce a warning.

# Arguments

- `sys`: The unsimplified system.
- `inputs`: The input variables of the linearization, or analysis points.
- `outputs`: The output variables of the linearization, or analysis points.

# Keyword Arguments

- `op`: The operating point, a dictionary of variable and parameter values or a
  [`LinearizationOpPoint`](@ref) wrapping a solution and a time. Signals crossing a clock
  boundary take the value of the variable they are derived from, and the history variables
  of a discrete partition take the value of the variable they are the history of, that is,
  the operating point is assumed to be stationary across ticks. Values that are not available
  default to zero.
- `t`: The time at which to linearize. Ignored if `op` is a `LinearizationOpPoint`.
- `autodiff`: The `ADTypes.AbstractADType` used to compute the Jacobians. Defaults to
  `AutoForwardDiff()`.
- `fallback_autodiff`: The AD type used when `autodiff` fails for a partition. Defaults to
  `AutoFiniteDiff()`; pass `nothing` to disable the fallback.
- `allow_input_derivatives`: See [`linearize`](@ref). Applies to the continuous partition.
- `warn_missing_op`: Whether to warn about operating-point values that default to zero.
- `warn_unsupported`: Whether to warn about events and clocks that are not accounted for.
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
    end

    # The model is flattened and its connections expanded as in `mtkcompile`; the inputs are
    # then bound to parameters before clock inference, so that they stay variables of their
    # partition.
    sys = expand_connections(sys)
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
    exposed_boundary_inputs = [SymbolicT[] for _ in 1:npart]
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
            push!(exposed_boundary_inputs[i], term)
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

    # Compile and linearize every partition.
    missing_vars = SymbolicT[]
    results = Vector{ClockPartitionLinearization}(undef, npart)
    for i in 1:npart
        discrete = i != continuous_id
        clk = id_to_clock[i]
        outs = [user_outputs[i]; boundary_out[i]]
        ssys, boundary_params = _compile_partition!(tss[i], boundary_inputs[i], outs, discrete, iv; kwargs...)
        term_to_param = Dict{SymbolicT, SymbolicT}(zip(boundary_inputs[i], boundary_params))
        ins = [user_input_params[i]; SymbolicT[term_to_param[term] for term in boundary_in[i]]]
        boundary_values = Pair{SymbolicT, SymbolicT}[
            term_to_param[term] => _boundary_source(term)
                for term in boundary_inputs[i] if !(operation(term) isa SampleTime)
        ]
        # Parameters of inputs belonging to other partitions are part of every partition's
        # parameter set and need a value as well.
        all_input_params = Pair{SymbolicT, SymbolicT}[param => var for (var, param) in zip(inputs, input_params)]
        pop = _partition_operating_point(ssys, hop, boundary_values, all_input_params, discrete, missing_vars)
        for term in sample_time_terms[i]
            pop[term_to_param[term]] = _sample_interval(clk)
        end
        description = discrete ? "clock partition on $clk" : "continuous partition"
        if discrete
            mats, x0 = _with_autodiff_fallback(autodiff, fallback_autodiff, description) do ad
                _linearize_discrete_partition(ssys, ins, outs, pop, t; autodiff = ad, eval_expression, eval_module)
            end
        else
            mats, x0 = _with_autodiff_fallback(autodiff, fallback_autodiff, description) do ad
                _linearize_continuous_partition(
                    ssys, ins, outs, pop, t; autodiff = ad, allow_input_derivatives,
                    initialize, initializealg, initialization_abstol, initialization_reltol,
                    initialization_solver_alg, guesses, missing_guess_value, eval_expression,
                    eval_module, loop_opening_params
                )
            end
        end
        results[i] = ClockPartitionLinearization(
            mats.A, mats.B, mats.C, mats.D, unknowns(ssys),
            [user_inputs[i]; boundary_in[i]], outs, length(user_inputs[i]), length(user_outputs[i]),
            clk, _sample_interval(clk), x0, ssys
        )
    end
    if warn_missing_op && !isempty(missing_vars)
        unique!(missing_vars)
        @warn "No operating-point value was found for the following variables; zero is used: $(join(string.(missing_vars), ", ")). Pass values in `op` or disable this warning with `warn_missing_op = false`."
    end

    # Order the partitions with the continuous one first and translate the connections to
    # indices into the input and output lists.
    order = continuous_id == 0 ? collect(1:npart) : [continuous_id; filter(!=(continuous_id), 1:npart)]
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
