module MTKFMIExt

using ModelingToolkit
using SymbolicIndexingInterface: NotSymbolic, symbolic_type
using ModelingToolkitBase: t_nounits, D_nounits
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ModelingToolkit as MTK
import ModelingToolkitBase as MTKBase
import SciMLBase
import SymbolicUtils
import FMIImport as FMI
import EzXML

# Spelled as `const` rather than `t_nounits as t`, because ExplicitImports reports a
# renaming import as stale (it tracks the pre-rename name, which never appears in use).
const t = t_nounits
const D = D_nounits

"""
    fmu_status_is_error(::Val{fmi_version}, status, fnname)

Whether `status`, returned by the FMI function `fnname`, is a failure. `Val{2}` and `Val{3}`
select the status constants of the FMI version in use. A warning status is not a failure: the
standard allows the computation to continue, so it is only logged.
"""
function fmu_status_is_error end

function fmu_status_is_error(::Val{2}, status, fnname)
    status === nothing && return false
    status == FMI.fmi2StatusOK && return false
    if status == FMI.fmi2StatusWarning
        @warn "FMU function $fnname returned `fmi2StatusWarning`." maxlog = 10
        return false
    end
    return true
end

function fmu_status_is_error(::Val{3}, status, fnname)
    status === nothing && return false
    status == FMI.fmi3StatusOK && return false
    if status == FMI.fmi3StatusWarning
        @warn "FMU function $fnname returned `fmi3StatusWarning`." maxlog = 10
        return false
    end
    return true
end

"""
    @statuscheck FMI.fmiXFunction(...)

Check the status an FMI.jl function returned, where `X` is `2` or `3`. On a failing status this
frees the instance and errors, terminating it first unless the status is fatal.

Evaluates to whatever the wrapped call returned, so callers keep the extra outputs of functions
like `FMI.fmi2CompletedIntegratorStep`. Those return their status as the first element of a
tuple, which is where this looks for it.
"""
macro statuscheck(expr)
    @assert Meta.isexpr(expr, :call)
    fn = expr.args[1]
    @assert Meta.isexpr(fn, :.)
    @assert fn.args[1] == :FMI
    # the qualified name is parsed as `Expr(:., :FMI, QuoteNode(:fmiXFunction))`
    fnname = fn.args[2] isa QuoteNode ? fn.args[2].value : fn.args[2]

    is_v2 = startswith(string(fnname), "fmi2")

    version = Val(is_v2 ? 2 : 3)
    return quote
        result = $expr
        check_status!(wrapper, $version, result isa Tuple ? result[1] : result, $(QuoteNode(fnname)))
        result
    end |> esc
end

"""
    $(TYPEDSIGNATURES)

The function behind `@statuscheck`, for calls whose FMI function is chosen at run time. On a
failing `status` this frees the instance and errors, terminating it first unless the status
is fatal.
"""
function check_status!(wrapper, version::Val{V}, status, fnname) where {V}
    fmu_status_is_error(version, status, fnname) || return
    # a fatal status leaves the instance in a state where no further call is allowed
    if status != (V == 2 ? FMI.fmi2StatusFatal : FMI.fmi3StatusFatal)
        V == 2 ? FMI.fmi2Terminate(wrapper.instance) : FMI.fmi3Terminate(wrapper.instance)
    end
    V == 2 ? FMI.fmi2FreeInstance!(wrapper.instance) : FMI.fmi3FreeInstance!(wrapper.instance)
    wrapper.instance = nothing
    error("FMU Error in $fnname: status $status")
end

"""
For each FMI 2 type FMIComponent can carry as a number: the element type FMI uses for it,
and its C-level setter and getter (`(instance, vrs, nvr, values)`).
"""
const FMI2_ACCESSORS = Dict(
    :Real => (FMI.fmi2Real, FMI.fmi2SetReal, FMI.fmi2GetReal!),
    :Integer => (FMI.fmi2Integer, FMI.fmi2SetInteger, FMI.fmi2GetInteger!),
    :Boolean => (FMI.fmi2Boolean, FMI.fmi2SetBoolean, FMI.fmi2GetBoolean!),
)

"""
The FMI 3 counterpart of `FMI2_ACCESSORS`; the C-level functions take
`(instance, vrs, nvr, values, nvalues)`. FMI 3 reads and writes enumerations as `Int64`.
"""
const FMI3_ACCESSORS = Dict(
    :Float64 => (FMI.fmi3Float64, FMI.fmi3SetFloat64, FMI.fmi3GetFloat64!),
    :Float32 => (FMI.fmi3Float32, FMI.fmi3SetFloat32, FMI.fmi3GetFloat32!),
    :Int8 => (FMI.fmi3Int8, FMI.fmi3SetInt8, FMI.fmi3GetInt8!),
    :Int16 => (FMI.fmi3Int16, FMI.fmi3SetInt16, FMI.fmi3GetInt16!),
    :Int32 => (FMI.fmi3Int32, FMI.fmi3SetInt32, FMI.fmi3GetInt32!),
    :UInt8 => (FMI.fmi3UInt8, FMI.fmi3SetUInt8, FMI.fmi3GetUInt8!),
    :UInt16 => (FMI.fmi3UInt16, FMI.fmi3SetUInt16, FMI.fmi3GetUInt16!),
    :UInt32 => (FMI.fmi3UInt32, FMI.fmi3SetUInt32, FMI.fmi3GetUInt32!),
    :Boolean => (FMI.fmi3Boolean, FMI.fmi3SetBoolean, FMI.fmi3GetBoolean!),
    :Enumeration => (FMI.fmi3Int64, FMI.fmi3SetInt64, FMI.fmi3GetInt64!),
)

function fmi_kind(::FMI.FMU2, mv)
    mv.Real === nothing || return :Real
    (mv.Integer === nothing && mv.Enumeration === nothing) || return :Integer
    mv.Boolean === nothing || return :Boolean
    return :String
end

fmi_kind(::FMI.FMU3, mv) = Symbol(chopprefix(string(nameof(typeof(mv))), "fmi3Variable"))
fmi_key(fmu, mv) = (fmi_kind(fmu, mv), mv.valueReference)

real_kind(::FMI.FMU2) = :Real
real_kind(::FMI.FMU3) = :Float64
accessor_table(::FMI.FMU2) = FMI2_ACCESSORS
accessor_table(::FMI.FMU3) = FMI3_ACCESSORS

causality_codes(::FMI.FMU2) =
    FMI.fmi2CausalityInput, FMI.fmi2CausalityOutput, FMI.fmi2CausalityParameter
causality_codes(::FMI.FMU3) =
    FMI.fmi3CausalityInput, FMI.fmi3CausalityOutput, FMI.fmi3CausalityParameter

const FMIReference = Tuple{Symbol, UInt32}

"""
    $(TYPEDEF)

Value references grouped by FMI type, with each group's positions in the MTK buffer.
"""
struct TypedValueReferences
    groups::Vector{Tuple{Symbol, Vector{Int}, Vector{UInt32}}}
end

function TypedValueReferences(keys::Vector{FMIReference})
    kinds = first.(keys)
    vrs = last.(keys)
    groups = Tuple{Symbol, Vector{Int}, Vector{UInt32}}[]
    for kind in unique(kinds)
        positions = findall(==(kind), kinds)
        push!(groups, (kind, positions, vrs[positions]))
    end
    return TypedValueReferences(groups)
end

to_fmi(::Type{T}, x) where {T <: AbstractFloat} = T(x)
to_fmi(::Type{T}, x) where {T <: Integer} = round(T, x)

"""
    $(TYPEDSIGNATURES)

Whether a model variable changes only at events: every non-float variable, and a float one the FMU
declares `discrete`. A Model Exchange FMU accepts a new value for such an input only in
Event Mode.
"""
is_discrete(fmu::FMI.FMU2, mv) =
    fmi_kind(fmu, mv) !== :Real || mv.variability == FMI.fmi2VariabilityDiscrete
is_discrete(fmu::FMI.FMU3, mv) =
    !(fmi_kind(fmu, mv) in (:Float64, :Float32)) || mv.variability == FMI.fmi3VariabilityDiscrete

"""
    $(TYPEDEF)

The inputs of a Model Exchange FMU, split into those that may change continuously and
those that change only at events.

# Fields

$(TYPEDFIELDS)
"""
struct InputSplit
    """
    Positions of the continuous inputs in the full input vector.
    """
    continuous_positions::Vector{Int}
    continuous::TypedValueReferences
    """
    Positions of the discrete inputs in the full input vector.
    """
    discrete_positions::Vector{Int}
    discrete::TypedValueReferences
    """
    The full input vector as last written to the FMU, to tell when a discrete input changed.
    """
    last::Vector{Float64}
end

function InputSplit(refs, discrete)
    cpos, dpos = findall(!, discrete), findall(discrete)
    return InputSplit(
        cpos, TypedValueReferences(refs[cpos]), dpos, TypedValueReferences(refs[dpos]),
        fill(NaN, length(refs))
    )
end

"""
    $(TYPEDSIGNATURES)

How many event functions a Model Exchange wrapper returns after the derivatives and
outputs: one per event indicator of the FMU, then one whose zero is the FMU's next time
event.
"""
number_of_events(fmu) = Int(something(FMI.getNumberOfEventIndicators(fmu), 0)) + 1

"""
    $(TYPEDSIGNATURES)

Write `values` to `wrapper.instance` with one typed FMI call per group.
"""
function set_values!(wrapper, refs::TypedValueReferences, values)
    version = fmi_version(wrapper)
    accessors = accessor_table(wrapper.fmu)
    for (kind, positions, vrs) in refs.groups
        T, set, _ = accessors[kind]
        converted = values isa Vector{T} && length(positions) == length(values) ?
            values : T[to_fmi(T, values[i]) for i in positions]
        n = Csize_t(length(vrs))
        status = version isa Val{2} ? set(wrapper.instance, vrs, n, converted) :
            set(wrapper.instance, vrs, n, converted, n)
        check_status!(wrapper, version, status, nameof(set))
    end
    return
end

"""
    $(TYPEDSIGNATURES)

Read the variables in `refs` from `wrapper.instance` into `out` with one typed FMI call per group.
"""
function get_values!(out, wrapper, refs::TypedValueReferences)
    version = fmi_version(wrapper)
    accessors = accessor_table(wrapper.fmu)
    for (kind, positions, vrs) in refs.groups
        T, _, get! = accessors[kind]
        values = out isa Vector{T} && length(positions) == length(out) ? out : zeros(T, length(vrs))
        n = Csize_t(length(vrs))
        status = version isa Val{2} ? get!(wrapper.instance, vrs, n, values) :
            get!(wrapper.instance, vrs, n, values, n)
        check_status!(wrapper, version, status, nameof(get!))
        values === out && continue
        for (j, i) in enumerate(positions)
            value = Float64(values[j])
            T !== Int64 || value == values[j] || throw(ArgumentError("FMI value $(values[j]) for $(kind) reference $(vrs[j]) cannot be represented exactly as Float64."))
            out[i] = value
        end
    end
    return out
end

"""
    $(TYPEDSIGNATURES)

A component that wraps an FMU loaded via FMI.jl. The FMI version (2 or 3) should be
provided as a `Val` to the function. Supports Model Exchange and CoSimulation FMUs.
Continuous states must be `FMI.fmi2Real` or `FMI.fmi3Float64`. FMI 2 Real, Integer,
Enumeration and Boolean ports and parameters and FMI 3 Float32, Float64, Int8, Int16, Int32,
UInt8, UInt16, UInt32, Boolean and exactly representable Enumeration values are carried
as `Float64`. FMI 3 Int64 and UInt64 ports and parameters are omitted with a warning because
their full ranges cannot be represented as `Float64`. String, binary and clock variables
are also omitted with a warning.

Model Exchange supports state, time, accepted-step and discrete input events. Discrete
outputs are updated at events. An FMU without continuous states gets a placeholder state.
Use a finite-difference solver option such as `Rodas5P(autodiff = AutoFiniteDiff())` for
solvers that need derivatives through a Model Exchange FMU.

Parameters default to their FMU start values. By default, supply continuous state values
in the problem map or determine them with initialization equations on the system or problem.
A Co-Simulation FMU accepts initial values only for states with a start value. It computes
states declared with `initial = "calculated"` itself. Hierarchical names such as
`namespace.variable` become `namespace__variable`.

# Keyword Arguments

- `fmu`: The FMU loaded via `FMI.loadFMU`.
- `tolerance`: The tolerance to provide to the FMU. Not used for v3 FMUs since it is not
  supported by FMI.jl.
- `communication_step_size`: The periodic interval at which communication with CoSimulation
  FMUs will occur. Required for Co-Simulation. A shorter final interval is only supported
  if the FMU declares variable communication step sizes.
- `reinitializealg`: The DAE initialization algorithm to use for the callback managing the
  FMU. For CoSimulation FMUs whose states/outputs are used in algebraic equations of the
  system, this needs to be an algorithm that will solve for the new algebraic variables.
  For example, `OrdinaryDiffEqCore.BrownFullBasicInit()`.
- `stop_time`: The FMU experiment limit. It defaults to `tspan[2]` for each solve and must
  be at or after that endpoint. The last communication occurs at the solve endpoint.
- `fmu_initial_states`: If `true`, the FMU's computed initial states are used as default
  equations for states without a problem-map value. Leave this off when supplying
  initialization equations that determine those states, since both sets of equations
  would constrain the initialization.
- `type`: Either `:ME` or `:CS` depending on whether `fmu` is a Model Exchange or
  CoSimulation FMU respectively.
- `name`: The name of the system.
"""
function MTK.FMIComponent(
        ::Val{Ver}; fmu = nothing, tolerance = 1.0e-6,
        communication_step_size = nothing, reinitializealg = nothing,
        stop_time = nothing, fmu_initial_states = false, type, name
    ) where {Ver}
    if Ver != 2 && Ver != 3
        throw(ArgumentError("FMI Version must be `2` or `3`"))
    end
    if type == :CS && !(communication_step_size isa Real && isfinite(communication_step_size) && communication_step_size > 0)
        throw(ArgumentError("`communication_step_size` must be positive and finite for Co-Simulation FMUs."))
    end
    md = fmu.modelDescription
    model_variables = Dict{FMIReference, Vector{eltype(md.modelVariables)}}()
    for mv in md.modelVariables
        push!(get!(() -> eltype(md.modelVariables)[], model_variables, fmi_key(fmu, mv)), mv)
    end
    function port_keys(causality)
        return unique(FMIReference[fmi_key(fmu, mv) for mv in md.modelVariables if mv.causality == causality])
    end
    input_causality, output_causality, parameter_causality = causality_codes(fmu)
    output_keys = port_keys(output_causality)
    state_keys = [(real_kind(fmu), vr) for vr in md.stateValueReferences]
    value_references = Dict{Any, FMIReference}()
    # defaults
    defs = Dict()
    # unknowns of the system
    states = []
    # differential variables of the system
    # this is a subset of `states` in the case where the FMU has multiple names for
    # the same value reference.
    diffvars = []
    # variables that are derivatives of diffvars
    dervars = []
    # observed equations
    observed = Equation[]
    # need to separate observed equations for duplicate derivative variables
    # since they aren't included in CS FMUs
    der_observed = Equation[]

    # parse states, in the order of the model description (`ModelStructure/Derivatives`
    # in FMI 2, `ContinuousStateDerivative` in FMI 3): the FMU's continuous-state functions
    # are positional, and `__mtk_internal_u` below is passed to them in this order.
    # A state the FMU also lists as an output is tagged as an output here rather than
    # declared again when the outputs are parsed.
    output_key_set = Set(output_keys)
    fmi_variables_to_mtk_variables!(
        fmu, model_variables, state_keys,
        value_references, diffvars, states, observed;
        postprocess_variable = (v, vr) -> vr in output_key_set ?
            MTKBase.setoutput(v, true) : v
    )
    # create a symbolic variable __mtk_internal_u to pass to the relevant registered
    # functions as the state vector
    if isempty(diffvars)
        # no differential variables
        __mtk_internal_u = Float64[]
    elseif type == :ME
        # to avoid running into `structural_simplify` warnings about array variables
        # and some unfortunate circular dependency issues, ME FMUs use an array of
        # symbolics instead. This is also not worse off in performance
        # because the former approach would allocate anyway.
        # TODO: Can we avoid an allocation here using static arrays?
        __mtk_internal_u = copy(diffvars)
    elseif type == :CS
        # CS FMUs do their own independent integration in a periodic callback, so their
        # unknowns are discrete variables in the `ODESystem`. A default of `missing` allows
        # them to be solved for during initialization.
        @discretes __mtk_internal_u(t)[1:length(diffvars)] = missing [guess = diffvars]
        push!(observed, __mtk_internal_u ~ copy(diffvars))
    end

    # parse derivatives of states, in the order of the model description: the FMU returns
    # all derivatives at once, in that order, whatever order their references are given in.
    # The variables passed to `postprocess_variable` haven't been differentiated yet, so they
    # should match one variable in states. That's the one this is the derivative of, and we
    # keep track of this ordering
    derivative_order = []
    function derivative_order_postprocess(var, _)
        idx = findfirst(isequal(var), states)
        idx === nothing || push!(derivative_order, states[idx])
        return var
    end
    derivative_keys = [(real_kind(fmu), vr) for vr in md.derivativeValueReferences]
    fmi_variables_to_mtk_variables!(
        fmu, model_variables, derivative_keys,
        value_references, dervars, states, der_observed;
        postprocess_variable = derivative_order_postprocess
    )
    @assert length(derivative_order) == length(dervars)

    # parse the inputs to the FMU
    inputs = []
    fmi_variables_to_mtk_variables!(
        fmu, model_variables, port_keys(input_causality),
        value_references, inputs, states, observed; postprocess_variable = (v, _) -> MTKBase.setinput(
            v, true
        )
    )
    # create a symbolic variable for the input buffer
    __mtk_internal_x = copy(inputs)
    if isempty(__mtk_internal_x)
        __mtk_internal_x = Float64[]
    end

    # parse the outputs of the FMU
    outputs = []
    fmi_variables_to_mtk_variables!(
        fmu, model_variables, setdiff(output_keys, state_keys),
        value_references, outputs, states, observed; postprocess_variable = (v, _) -> MTKBase.setoutput(
            v, true
        )
    )
    # A Model Exchange FMU's discrete outputs change only at events, which can depend on
    # discrete state inside the FMU that the system cannot see. They are kept as discrete
    # variables, set when the system starts and after every FMU event, so a solution
    # reports what the FMU said during the solve rather than asking a fresh instance.
    discrete_outputs = []
    is_discrete_port(v) = is_discrete(fmu, first(model_variables[value_references[v]]))
    if type == :ME
        discrete_outputs = filter(is_discrete_port, outputs)
        outputs = filter(!is_discrete_port, outputs)
        if !isempty(discrete_outputs)
            @discretes __mtk_internal_do(t)[1:length(discrete_outputs)] = missing [guess = zeros(length(discrete_outputs))]
            push!(observed, __mtk_internal_do ~ discrete_outputs)
        end
    end

    # create the output buffer. This is only required for CoSimulation to pass it to
    # the callback affect
    if type == :CS
        if isempty(outputs)
            __mtk_internal_o = Float64[]
        else
            @discretes __mtk_internal_o(t)[1:length(outputs)] = missing [guess = zeros(length(outputs))]
            push!(observed, __mtk_internal_o ~ outputs)
        end
    end

    # parse the parameters
    params = []
    # multiple names for the same parameter are treated as parameter dependencies.
    parameter_dependencies = Equation[]
    fmi_variables_to_mtk_variables!(
        fmu, model_variables, port_keys(parameter_causality), value_references,
        params, [], parameter_dependencies, defs; parameters = true
    )
    # create a symbolic variable for the parameter buffer
    __mtk_internal_p = copy(params)
    if isempty(__mtk_internal_p)
        __mtk_internal_p = Float64[]
    end

    state_value_references = UInt32[last(value_references[var]) for var in diffvars]
    output_value_references = FMIReference[value_references[var] for var in outputs]
    discrete_output_value_references = FMIReference[value_references[var] for var in discrete_outputs]
    input_value_references = FMIReference[value_references[var] for var in inputs]
    param_value_references = FMIReference[value_references[var] for var in params]
    discrete_inputs = Bool[is_discrete_port(var) for var in inputs]

    # create a parameter for the instance wrapper
    # this manages the creation and deallocation of FMU instances
    buffer_length = length(diffvars) + length(outputs) + length(discrete_outputs) +
        number_of_events(fmu)
    if Ver == 2
        @parameters (wrapper::FMI2InstanceWrapper)(..)[1:buffer_length] = FMI2InstanceWrapper(
            fmu, state_value_references, output_value_references,
            param_value_references, input_value_references, discrete_inputs, tolerance,
            discrete_output_value_references
        )
    else
        @parameters (wrapper::FMI3InstanceWrapper)(..)[1:buffer_length] = FMI3InstanceWrapper(
            fmu, state_value_references,
            output_value_references, param_value_references, input_value_references, discrete_inputs,
            discrete_output_value_references
        )
    end

    if fmu_initial_states && !isempty(diffvars)
        @parameters (initial_states::FMUInitialStates)(..)[1:length(diffvars)] = FMUInitialStates(type)
        initial_states_expr = initial_states(wrapper, __mtk_internal_x, __mtk_internal_p, t)
        for (i, var) in enumerate(diffvars)
            defs[var] = initial_states_expr[i]
        end
    end

    # any additional initialization equations for the system
    initialization_eqs = Equation[]
    continuous_events = []
    discrete_events = []

    if type == :ME
        # the wrapper is a callable struct which returns the state derivative and
        # output values
        # symbolic expression for calling the wrapper
        call_expr = wrapper(__mtk_internal_u, __mtk_internal_x, __mtk_internal_p, t)

        # differential and observed equations
        diffeqs = Equation[]
        for (i, var) in enumerate([dervars; outputs])
            push!(diffeqs, var ~ call_expr[i])
        end
        for (var, dervar) in zip(derivative_order, dervars)
            push!(diffeqs, D(var) ~ dervar)
        end
        # The event functions below need a continuous callback, which ModelingToolkit
        # cannot compile for a system with no unknowns. An FMU with no continuous states
        # gets one placeholder state that never changes, so it can run on its own. It is
        # marked irreducible because `mtkcompile` otherwise removes a constant state.
        if isempty(diffvars)
            placeholder = only(@variables __mtk_internal_placeholder(t) = 0.0 [irreducible = true])
            push!(diffeqs, D(placeholder) ~ 0)
            push!(states, placeholder)
        end

        # the discrete outputs take their first values from the FMU during initialization
        offset = length(dervars) + length(outputs)
        for i in eachindex(discrete_outputs)
            push!(initialization_eqs, __mtk_internal_do[i] ~ call_expr[offset + i])
        end
        discrete_modified = isempty(discrete_outputs) ? (;) :
            (; discrete_outputs = __mtk_internal_do)

        # event functions: the FMU's event indicators, then one that crosses zero at its
        # next time event. A zero crossing of any runs the FMU's event iteration.
        offset += length(discrete_outputs)
        event_vars = [
            only(@variables $(Symbol(:__mtk_internal_event_, k))(t))
                for k in 1:number_of_events(fmu)
        ]
        for (k, var) in enumerate(event_vars)
            push!(diffeqs, var ~ call_expr[offset + k])
            push!(states, var)
        end
        me_observed = (; wrapper, inputs = __mtk_internal_x, params = copy(params))
        me_modified = (; (isempty(diffvars) ? (;) : (; states = diffvars))..., discrete_modified...)
        for (k, var) in enumerate(event_vars)
            # A signed index encodes the crossing direction; zero denotes the time event.
            root = k == length(event_vars) ? 0 : k
            affect = MTKBase.ImperativeAffect(
                fmiMEEvent!; observed = me_observed, modified = me_modified, ctx = root
            )
            affect_neg = root == 0 ? nothing : MTKBase.ImperativeAffect(
                    fmiMEEvent!; observed = me_observed, modified = me_modified, ctx = -root
                )
            push!(
                continuous_events, MTKBase.SymbolicContinuousCallback(
                    [var ~ 0], affect; affect_neg, rootfind = SciMLBase.RightRootFind, reinitializealg
                )
            )
        end
        isempty(discrete_outputs) || push!(params, __mtk_internal_do)

        @parameters (step_due::MEStepDue)(..) = MEStepDue()
        push!(params, step_due)
        initialize_affect = MTKBase.ImperativeAffect(
            fmiMEInitialize!; observed = me_observed, modified = discrete_modified
        )
        finalize_affect = MTKBase.ImperativeAffect(fmiMEFinalize!; observed = (; wrapper))
        step_affect = MTKBase.ImperativeAffect(
            fmiMEEvent!; observed = me_observed, modified = me_modified
        )
        instance_management_callback = MTKBase.SymbolicDiscreteCallback(
            step_due(wrapper, __mtk_internal_x, t) > 0.5, step_affect; initialize = initialize_affect,
            finalize = finalize_affect, reinitializealg
        )

        push!(params, wrapper)
        append!(observed, der_observed)
    elseif type == :CS
        _functor = CSFunctor(length(output_value_references), something(stop_time, NaN))
        @parameters (functor::CSFunctor)(..)[1:length(__mtk_internal_o)] = _functor
        # for co-simulation, we need to ensure the output buffer is solved for
        # during initialization
        for (i, x) in enumerate(collect(__mtk_internal_o))
            push!(
                initialization_eqs,
                x ~ functor(
                    wrapper, __mtk_internal_u, __mtk_internal_x, __mtk_internal_p, t
                )[i]
            )
        end

        diffeqs = Equation[]

        # use `ImperativeAffect` for instance management here
        cb_observed = (;
            inputs = __mtk_internal_x, params = copy(params),
            t, wrapper, dt = communication_step_size,
        )
        cb_modified = (;)
        # modify the outputs if present
        if symbolic_type(__mtk_internal_o) != NotSymbolic()
            cb_modified = (cb_modified..., outputs = __mtk_internal_o)
        end
        # modify the continuous state if present
        if symbolic_type(__mtk_internal_u) != NotSymbolic()
            cb_modified = (cb_modified..., states = __mtk_internal_u)
        end
        initialize_affect = MTKBase.ImperativeAffect(
            fmiCSInitialize!; observed = cb_observed,
            modified = cb_modified, ctx = _functor
        )
        finalize_affect = MTKBase.ImperativeAffect(fmiFinalize!; observed = (; wrapper))
        # the callback affect performs the stepping
        step_affect = MTKBase.ImperativeAffect(
            fmiCSStep!; observed = cb_observed, modified = cb_modified, ctx = _functor
        )
        instance_management_callback = MTKBase.SymbolicDiscreteCallback(
            communication_step_size, step_affect; initialize = initialize_affect,
            finalize = finalize_affect, reinitializealg
        )
        # the periodic callback does not fire at the final time; this takes the last step
        @parameters (final_step_due::CSFinalStepDue)(..) = CSFinalStepDue()
        push!(params, final_step_due)
        push!(
            discrete_events, MTKBase.SymbolicDiscreteCallback(
                final_step_due(functor, t) > 0.5,
                MTKBase.ImperativeAffect(
                    fmiCSStep!; observed = cb_observed, modified = cb_modified, ctx = _functor
                ); reinitializealg
            )
        )

        # guarded in case there are no outputs/states and the variable is `[]`.
        symbolic_type(__mtk_internal_o) == NotSymbolic() || push!(params, __mtk_internal_o)
        symbolic_type(__mtk_internal_u) == NotSymbolic() || push!(params, __mtk_internal_u)

        push!(params, wrapper, functor)
    end

    # added last: the branches above copy `params` as the FMU's own parameter values
    fmu_initial_states && !isempty(diffvars) && push!(params, initial_states)
    eqs = [observed; diffeqs]
    bindings = [eq.lhs => eq.rhs for eq in parameter_dependencies]
    return System(
        eqs, t, states, params; bindings, initial_conditions = defs,
        discrete_events = [instance_management_callback; discrete_events], continuous_events, name,
        initialization_eqs
    )
end

"""
    $(TYPEDSIGNATURES)

A symbolic variable is created for every FMU variable name associated with an entry in
`keys`. The typed value reference is stored in `value_references`. For aliases, the
symbolic variable for the first name is pushed
to `truevars`. All of the created symbolic variables are pushed to `allvars`. Observed
equations equating identical variables are pushed to `obseqs`. `defs` is a dictionary of
defaults.

# Keyword Arguments
- `parameters`: A boolean indicating whether to use `@parameters` for the symbolic
  variables instead of `@variables`.
- `postprocess_variable`: A function of the created variable and its value reference, applied to each created variable that should
  return the updated variable. This is useful to add metadata to variables.
"""
function fmi_variables_to_mtk_variables!(
        fmu, model_variables, keys,
        value_references::AbstractDict, truevars, allvars,
        obseqs, defs = Dict(); parameters = false, postprocess_variable = (v, _) -> v
    )
    for valRef in keys
        records = model_variables[valRef]
        varnames = [mv.name for mv in records]
        kind = first(valRef)
        if kind in (:Int64, :UInt64)
            @warn "FMI 3 $kind variables $(join(varnames, ", ")) are unsupported because the Float64 transport cannot represent all $kind values."
            continue
        elseif !haskey(accessor_table(fmu), kind)
            @warn "FMU variable $(join(varnames, ", ")) is not numeric and is left out of the component."
            continue
        end
        stateT = Float64
        snames = Symbol[]
        ders = Int[]
        for name in varnames
            sname, der = parseFMIVariableName(name)
            push!(snames, sname)
            push!(ders, der)
        end
        if parameters
            vars = [
                postprocess_variable(SymbolicUtils.unwrap(only(@parameters $sname::stateT)), valRef)
                    for sname in snames
            ]
        else
            vars = [
                postprocess_variable(SymbolicUtils.unwrap(only(@variables $sname(t)::stateT)), valRef)
                    for sname in snames
            ]
        end
        for i in eachindex(vars)
            der = ders[i]
            vars[i] = SymbolicUtils.unwrap(vars[i])
            for j in 1:der
                vars[i] = D(vars[i])
            end
            vars[i] = MTKBase.default_toterm(vars[i])
        end
        for i in eachindex(vars)
            if i == 1
                push!(truevars, vars[i])
            else
                push!(obseqs, vars[i] ~ vars[1])
            end
            value_references[vars[i]] = valRef
        end
        append!(allvars, vars)
        if parameters
            defval = FMI.getStartValue(first(records))
            if defval isa Integer && Float64(defval) != defval
                throw(ArgumentError("The start value of FMU parameter $(first(varnames)) cannot be represented exactly as Float64."))
            end
            defs[vars[1]] = defval isa Number ? Float64(defval) : defval
        end
    end
    return
end

"""
    $(TYPEDSIGNATURES)

Parse the string name of an FMI variable into a `Symbol` name for the corresponding
MTK variable. Return the `Symbol` name and the number of times it is differentiated.
"""
function parseFMIVariableName(name::AbstractString)
    name = replace(name, "." => "__")
    der = 0
    if startswith(name, "der(")

        # account for multi-dimensional array variable derivatives, e.g. der(x[1,2], 2)
        array_variable_pattern = r"\[\d+?,\d+?\]"
        patternmatches = match(array_variable_pattern, name)
        if (patternmatches !== nothing)
            safe_array_index_str = replace(String(patternmatches.match), "," => "_")
            safe_name = replace(name, array_variable_pattern => safe_array_index_str)
        else
            safe_name = name
        end


        idx = findfirst(',', safe_name)
        if idx === nothing
            name = @view name[5:(end - 1)]
            der = 1
        else

            der = parse(Int, @view name[(idx + 1):(end - 1)])
            name = @view name[5:(idx - 1)]
        end
    end
    return Symbol(name), der
end

"""
    $(TYPEDEF)

A struct which manages instance creation and deallocation for v2 FMUs.

# Fields

$(TYPEDFIELDS)
"""
mutable struct FMI2InstanceWrapper
    """
    The FMU from `FMI.loadFMU`.
    """
    const fmu::FMI.FMU2
    """
    The value references for the states of the FMU.
    """
    const state_value_references::Vector{FMI.fmi2ValueReference}
    """
    The tolerance with which to setup the FMU instance.
    """
    const tolerance::FMI.fmi2Real
    """
    The FMU instance, if present, and `nothing` otherwise.
    """
    instance::Union{FMI.FMU2Component{FMI.FMU2}, Nothing}
    """
    Continuous state buffer pre-allocated to avoid simulation allocations.
    """
    const states_buffer::Vector{FMI.fmi2Real}
    """
    Output buffer pre-allocated to avoid simulation allocations.
    """
    const outputs_buffer::Vector{FMI.fmi2Real}
    """
    Return buffer caching state and output arrays to eliminate return vcat allocations.
    """
    const res_buffer::Vector{FMI.fmi2Real}
    """
    Ports and parameters grouped by FMI type.
    """
    const outputs::TypedValueReferences
    const params::TypedValueReferences
    const inputs::TypedValueReferences
    const input_split::InputSplit
    """
    Buffer for the FMU's event indicators (Model Exchange).
    """
    const event_indicators::Vector{FMI.fmi2Real}
    """
    The discrete outputs of a Model Exchange FMU, grouped by FMI type, and a buffer for
    their values. They change only at events and are kept by the system as discrete
    variables, so that a solution reports the values the FMU had during the solve.
    """
    const discrete_outputs::TypedValueReferences
    const discrete_outputs_buffer::Vector{Float64}
    """
    The time of the FMU's next time event, or `NaN` if it has none scheduled (Model
    Exchange).
    """
    next_event_time::Float64
    """
    The parameter values the current instance was created with.
    """
    const instance_params::Vector{Float64}
    """
    Whether a solve is running (Model Exchange): between solves, an instance whose
    parameters differ from the ones asked for is replaced rather than reused.
    """
    in_solve::Bool
    const needs_completed_step::Bool
    last_completed_step::Float64
end

"""
    $(TYPEDSIGNATURES)

Create an `FMI2InstanceWrapper` with no instance.
"""
function FMI2InstanceWrapper(
        fmu, states, outputs, params, inputs, discrete_inputs, tolerance, discrete_outputs
    )
    n = length(states) + length(outputs) + length(discrete_outputs) + number_of_events(fmu)
    return FMI2InstanceWrapper(
        fmu, states, tolerance, nothing,
        zeros(FMI.fmi2Real, length(states)), zeros(FMI.fmi2Real, length(outputs)),
        zeros(FMI.fmi2Real, n),
        TypedValueReferences(outputs), TypedValueReferences(params),
        TypedValueReferences(inputs), InputSplit(inputs, discrete_inputs),
        zeros(FMI.fmi2Real, number_of_events(fmu) - 1),
        TypedValueReferences(discrete_outputs), zeros(length(discrete_outputs)), NaN,
        zeros(length(params)), false, needs_completed_step(fmu), NaN
    )
end

Base.nameof(::FMI2InstanceWrapper) = :FMI2InstanceWrapper
fmi_version(::FMI2InstanceWrapper) = Val(2)

"""
    $(TYPEDSIGNATURES)

Common functionality for creating an instance of a v2 FMU. Does not check if
`wrapper.instance` is already present, and overwrites the existing value with
a new instance. `inputs` should be in the order used when the component was constructed, as should `params`. `t` is the
current time. Returns the created instance, which is also stored in `wrapper.instance`.
"""
function get_instance_common!(wrapper::FMI2InstanceWrapper, inputs, params, t, t_end = nothing)
    wrapper.instance = FMI.fmi2Instantiate!(wrapper.fmu)::FMI.FMU2Component
    wrapper.instance_params .= params
    set_values!(wrapper, wrapper.inputs, inputs)
    wrapper.input_split.last .= inputs
    set_values!(wrapper, wrapper.params, params)
    stop_time_defined = t_end === nothing ? FMI.fmi2False : FMI.fmi2True
    stop_time = something(t_end, t)
    @statuscheck FMI.fmi2SetupExperiment(
        wrapper.instance, FMI.fmi2True, wrapper.tolerance, t, stop_time_defined, stop_time
    )
    @statuscheck FMI.fmi2EnterInitializationMode(wrapper.instance)
    return wrapper.instance
end

"""
    $(TYPEDSIGNATURES)

Create an instance of a Model Exchange FMU. Use the existing instance in `wrapper` if
present and create a new one otherwise. Return the instance.

See `get_instance_common!` for a description of the arguments.
"""
function get_instance_ME!(wrapper::FMI2InstanceWrapper, inputs, params, t)
    replace_stale_instance!(wrapper, params)
    if wrapper.instance === nothing
        get_instance_common!(wrapper, inputs, params, t)
        @statuscheck FMI.fmi2ExitInitializationMode(wrapper.instance)
        # the instance is now in Event Mode
        leave_event_mode!(wrapper)
    end

    return wrapper.instance
end

"""
    $(TYPEDSIGNATURES)

Create an instance of a CoSimulation FMU. Use the existing instance in `wrapper` if
present and create a new one otherwise. Return the instance.

See `get_instance_common!` for a description of the arguments.
"""
function get_instance_CS!(wrapper::FMI2InstanceWrapper, states, inputs, params, t, t_end)
    if wrapper.instance === nothing
        get_instance_common!(wrapper, inputs, params, t, t_end)
        settable = settable_states(wrapper)
        if !isempty(states) && !isempty(settable)
            @statuscheck FMI.fmi2SetReal(
                wrapper.instance, wrapper.state_value_references[settable],
                Csize_t(length(settable)), collect(Float64, states[settable])
            )
        end
        @statuscheck FMI.fmi2ExitInitializationMode(wrapper.instance)
    else
        # A Co-Simulation FMU integrates its own states: they are written only to start it.
        # Between steps the FMI standard lets the caller write inputs (and tunable
        # parameters), not states, and FMUs such as BouncingBall reject the write.
        set_values!(wrapper, wrapper.inputs, inputs)
    end
    return wrapper.instance
end

function completed_step!(wrapper::FMI2InstanceWrapper)
    _, enter_event, terminate = @statuscheck FMI.fmi2CompletedIntegratorStep(wrapper.instance, FMI.fmi2False)
    return enter_event == FMI.fmi2True, terminate == FMI.fmi2True
end

"""
    $(TYPEDSIGNATURES)

If `wrapper.instance !== nothing`, terminate and free the instance. Also set
`wrapper.instance` to `nothing`.
"""
function reset_instance!(wrapper::FMI2InstanceWrapper)
    wrapper.instance === nothing && return
    FMI.fmi2Terminate(wrapper.instance)
    FMI.fmi2FreeInstance!(wrapper.instance)
    return wrapper.instance = nothing
end

"""
    $(TYPEDEF)

A struct which manages instance creation and deallocation for v3 FMUs.

# Fields

$(TYPEDFIELDS)
"""
mutable struct FMI3InstanceWrapper
    """
    The FMU from `FMI.loadFMU`.
    """
    const fmu::FMI.FMU3
    """
    The value references for the states of the FMU.
    """
    const state_value_references::Vector{FMI.fmi3ValueReference}
    """
    The FMU instance, if present, and `nothing` otherwise.
    """
    instance::Union{FMI.FMU3Instance{FMI.FMU3}, Nothing}
    """
    Continuous state buffer pre-allocated to avoid simulation allocations.
    """
    const states_buffer::Vector{FMI.fmi3Float64}
    """
    Output buffer pre-allocated to avoid simulation allocations.
    """
    const outputs_buffer::Vector{FMI.fmi3Float64}
    """
    Return buffer caching state and output arrays to eliminate return vcat allocations.
    """
    const res_buffer::Vector{FMI.fmi3Float64}
    """
    Ports and parameters grouped by FMI type.
    """
    const outputs::TypedValueReferences
    const params::TypedValueReferences
    const inputs::TypedValueReferences
    const input_split::InputSplit
    """
    Buffer for the FMU's event indicators (Model Exchange).
    """
    const event_indicators::Vector{FMI.fmi3Float64}
    """
    The discrete outputs of a Model Exchange FMU, grouped by FMI type, and a buffer for
    their values. They change only at events and are kept by the system as discrete
    variables, so that a solution reports the values the FMU had during the solve.
    """
    const discrete_outputs::TypedValueReferences
    const discrete_outputs_buffer::Vector{Float64}
    """
    The time of the FMU's next time event, or `NaN` if it has none scheduled (Model
    Exchange).
    """
    next_event_time::Float64
    """
    The parameter values the current instance was created with.
    """
    const instance_params::Vector{Float64}
    """
    Whether a solve is running (Model Exchange): between solves, an instance whose
    parameters differ from the ones asked for is replaced rather than reused.
    """
    in_solve::Bool
    const needs_completed_step::Bool
    last_completed_step::Float64
end

"""
    $(TYPEDSIGNATURES)

Create an `FMI3InstanceWrapper` with no instance.
"""
function FMI3InstanceWrapper(
        fmu, states, outputs, params, inputs, discrete_inputs, discrete_outputs
    )
    n = length(states) + length(outputs) + length(discrete_outputs) + number_of_events(fmu)
    return FMI3InstanceWrapper(
        fmu, states, nothing,
        zeros(FMI.fmi3Float64, length(states)), zeros(FMI.fmi3Float64, length(outputs)),
        zeros(FMI.fmi3Float64, n),
        TypedValueReferences(outputs), TypedValueReferences(params),
        TypedValueReferences(inputs), InputSplit(inputs, discrete_inputs),
        zeros(FMI.fmi3Float64, number_of_events(fmu) - 1),
        TypedValueReferences(discrete_outputs), zeros(length(discrete_outputs)), NaN,
        zeros(length(params)), false, needs_completed_step(fmu), NaN
    )
end

Base.nameof(::FMI3InstanceWrapper) = :FMI3InstanceWrapper
fmi_version(::FMI3InstanceWrapper) = Val(3)

"""
    $(TYPEDSIGNATURES)

Common functionality for creating an instance of a v3 FMU. Since v3 FMUs need to be
instantiated differently depending on the type, this assumes `wrapper.instance` is a
freshly instantiated FMU. `inputs` and `params` follow the order used to construct the
component. `t` is the current time. Returns `wrapper.instance`.
"""
function get_instance_common!(wrapper::FMI3InstanceWrapper, inputs, params, t, t_end = nothing)
    wrapper.instance_params .= params
    set_values!(wrapper, wrapper.params, params)
    stop_time_defined = t_end === nothing ? FMI.fmi3False : FMI.fmi3True
    stop_time = something(t_end, t)
    @statuscheck FMI.fmi3EnterInitializationMode(
        wrapper.instance, FMI.fmi3False, zero(FMI.fmi3Float64), t, stop_time_defined, stop_time
    )
    set_values!(wrapper, wrapper.inputs, inputs)
    wrapper.input_split.last .= inputs

    return wrapper.instance
end

"""
    $(TYPEDSIGNATURES)

Create an instance of a Model Exchange FMU. Use the existing instance in `wrapper` if
present and create a new one otherwise. Return the instance.

See `get_instance_common!` for a description of the arguments.
"""
function get_instance_ME!(wrapper::FMI3InstanceWrapper, inputs, params, t)
    replace_stale_instance!(wrapper, params)
    if wrapper.instance === nothing
        wrapper.instance = FMI.fmi3InstantiateModelExchange!(wrapper.fmu)::FMI.FMU3Instance
        get_instance_common!(wrapper, inputs, params, t)
        @statuscheck FMI.fmi3ExitInitializationMode(wrapper.instance)
        # the instance is now in Event Mode
        leave_event_mode!(wrapper)
    end

    return wrapper.instance
end

"""
    $(TYPEDSIGNATURES)

Create an instance of a CoSimulation FMU. Use the existing instance in `wrapper` if
present and create a new one otherwise. Return the instance.

See `get_instance_common!` for a description of the arguments.
"""
function get_instance_CS!(wrapper::FMI3InstanceWrapper, states, inputs, params, t, t_end)
    if wrapper.instance === nothing
        wrapper.instance = FMI.fmi3InstantiateCoSimulation!(
            wrapper.fmu; eventModeUsed = false
        )::FMI.FMU3Instance
        get_instance_common!(wrapper, inputs, params, t, t_end)
        settable = settable_states(wrapper)
        if !isempty(states) && !isempty(settable)
            @statuscheck FMI.fmi3SetFloat64(
                wrapper.instance, wrapper.state_value_references[settable],
                collect(Float64, states[settable])
            )
        end
        @statuscheck FMI.fmi3ExitInitializationMode(wrapper.instance)
    else
        # see the FMI 2 method: states are written only to start the FMU
        set_values!(wrapper, wrapper.inputs, inputs)
    end
    return wrapper.instance
end

"""
    $(TYPEDSIGNATURES)
"""
function completed_step!(wrapper::FMI3InstanceWrapper)
    enterEventMode = Ref(FMI.fmi3False)
    terminateSimulation = Ref(FMI.fmi3False)
    @statuscheck FMI.fmi3CompletedIntegratorStep!(
        wrapper.instance, FMI.fmi3False, enterEventMode, terminateSimulation
    )
    return enterEventMode[] == FMI.fmi3True, terminateSimulation[] == FMI.fmi3True
end

"""
    $(TYPEDSIGNATURES)
"""
function reset_instance!(wrapper::FMI3InstanceWrapper)
    wrapper.instance === nothing && return
    FMI.fmi3Terminate(wrapper.instance)
    FMI.fmi3FreeInstance!(wrapper.instance)
    return wrapper.instance = nothing
end

@register_array_symbolic (fn::FMI2InstanceWrapper)(
    states::Vector{<:Real}, inputs::Vector{<:Real}, params::Vector{<:Real}, t::Real
) begin
    size = (length(fn.res_buffer),)
    eltype = eltype(states)
    ndims = 1
end

@register_array_symbolic (fn::FMI3InstanceWrapper)(
    states::Vector{<:Real}, inputs::Vector{<:Real}, params::Vector{<:Real}, t::Real
) begin
    size = (length(fn.res_buffer),)
    eltype = eltype(states)
    ndims = 1
end

function model_exchange_flag(fmu, name, default)
    md = fmu.modelDescription.modelExchange
    md === nothing && return default
    key = Symbol(name)
    if isdefined(md, key) && getproperty(md, key) !== nothing
        return getproperty(md, key)
    end
    # FMIImport may leave this capability unset in the parsed description.
    doc = EzXML.readxml(joinpath(fmu.path, "modelDescription.xml"))
    mode = findfirst("/fmiModelDescription/ModelExchange", doc)
    return mode === nothing || !haskey(mode, name) ? default : mode[name] in ("true", "1")
end

needs_completed_step(fmu::FMI.FMU2) = !model_exchange_flag(fmu, "completedIntegratorStepNotNeeded", false)
needs_completed_step(fmu::FMI.FMU3) = model_exchange_flag(fmu, "needsCompletedIntegratorStep", false)

function enter_event_mode!(wrapper::FMI2InstanceWrapper, step_event, root)
    return @statuscheck FMI.fmi2EnterEventMode(wrapper.instance)
end
function enter_event_mode!(wrapper::FMI3InstanceWrapper, step_event, root)
    time_event = root === 0
    state_event = root !== nothing && !time_event
    roots = zeros(FMI.fmi3Int32, length(wrapper.event_indicators))
    state_event && (roots[abs(root)] = sign(root))
    return @statuscheck FMI.fmi3EnterEventMode(
        wrapper.instance, FMI.fmi3Boolean(step_event), FMI.fmi3Boolean(state_event), roots,
        Csize_t(length(roots)), FMI.fmi3Boolean(time_event)
    )
end

"""
    $(TYPEDSIGNATURES)

Iterate the FMU's discrete states until they settle (the event iteration of the FMI
standard), then return to Continuous-Time Mode. Returns whether the FMU changed the values
of its continuous states.
"""
function leave_event_mode!(wrapper::FMI2InstanceWrapper, integrator = nothing)
    info = FMI.fmi2EventInfo()
    info.newDiscreteStatesNeeded = FMI.fmi2True
    changed = false
    while info.newDiscreteStatesNeeded == FMI.fmi2True
        @statuscheck FMI.fmi2NewDiscreteStates!(wrapper.instance, info)
        changed |= info.valuesOfContinuousStatesChanged == FMI.fmi2True
        if info.terminateSimulation == FMI.fmi2True
            integrator === nothing && error("The FMU asked to terminate during initialization.")
            SciMLBase.terminate!(integrator)
            return changed
        end
    end
    wrapper.next_event_time = info.nextEventTimeDefined == FMI.fmi2True ? info.nextEventTime : NaN
    @statuscheck FMI.fmi2EnterContinuousTimeMode(wrapper.instance)
    return changed
end
function leave_event_mode!(wrapper::FMI3InstanceWrapper, integrator = nothing)
    needed = Ref(FMI.fmi3True)
    terminate = Ref(FMI.fmi3False)
    nominals = Ref(FMI.fmi3False)
    values = Ref(FMI.fmi3False)
    next_defined = Ref(FMI.fmi3False)
    next_time = Ref(zero(FMI.fmi3Float64))
    changed = false
    while needed[] == FMI.fmi3True
        @statuscheck FMI.fmi3UpdateDiscreteStates(
            wrapper.instance, needed, terminate, nominals, values, next_defined, next_time
        )
        changed |= values[] == FMI.fmi3True
        if terminate[] == FMI.fmi3True
            integrator === nothing && error("The FMU asked to terminate during initialization.")
            SciMLBase.terminate!(integrator)
            return changed
        end
    end
    wrapper.next_event_time = next_defined[] == FMI.fmi3True ? next_time[] : NaN
    @statuscheck FMI.fmi3EnterContinuousTimeMode(wrapper.instance)
    return changed
end

function discrete_inputs_changed(wrapper, inputs)
    split = wrapper.input_split
    return any(i -> !isequal(inputs[i], split.last[i]), split.discrete_positions)
end

struct MEStepDue end

function (::MEStepDue)(wrapper, inputs, time)
    needs_completion = wrapper.needs_completed_step && time != wrapper.last_completed_step
    return (needs_completion || discrete_inputs_changed(wrapper, inputs)) ? 1.0 : 0.0
end

"""
    $(TYPEDSIGNATURES)

Read the FMU's event indicators into `wrapper.event_indicators`.
"""
function get_event_indicators!(wrapper::FMI2InstanceWrapper)
    n = length(wrapper.event_indicators)
    n == 0 && return
    return @statuscheck FMI.fmi2GetEventIndicators!(wrapper.instance, wrapper.event_indicators, Csize_t(n))
end
function get_event_indicators!(wrapper::FMI3InstanceWrapper)
    n = length(wrapper.event_indicators)
    n == 0 && return
    return @statuscheck FMI.fmi3GetEventIndicators!(wrapper.instance, wrapper.event_indicators, Csize_t(n))
end

"""
    $(TYPEDSIGNATURES)

Synchronize a Model Exchange FMU with the solver at an accepted step or a localized
event. Notify the FMU of the completed step if required, then handle any event and
write back the FMU's state resets and discrete outputs.
"""
function fmiMEEvent!(m, o, root, integrator)
    wrapper = o.wrapper
    time = Float64(integrator.t)
    inputs = as_fmi_vector(o.inputs)
    states = haskey(m, :states) ? as_fmi_vector(m.states) : Float64[]
    get_instance_ME!(wrapper, inputs, as_fmi_vector(o.params), time)
    sync_me_point!(wrapper, states, inputs, time)

    step_event = false
    if time != wrapper.last_completed_step
        wrapper.last_completed_step = time
        if wrapper.needs_completed_step
            step_event, terminate = completed_step!(wrapper)
            if terminate
                SciMLBase.terminate!(integrator)
                return (;)
            end
        end
    end
    input_event = discrete_inputs_changed(wrapper, inputs)
    if !(step_event || input_event || root !== nothing)
        SciMLBase.derivative_discontinuity!(integrator, false)
        return (;)
    end

    enter_event_mode!(wrapper, step_event, root)
    if input_event
        split = wrapper.input_split
        set_values!(wrapper, split.discrete, view(inputs, split.discrete_positions))
        split.last .= inputs
    end
    changed = leave_event_mode!(wrapper, integrator)
    updates = with_discrete_outputs(m, wrapper)
    if changed && haskey(m, :states)
        reset_states = copy(states)
        get_me_states!(wrapper, reset_states)
        return (; updates..., states = reset_states)
    end
    return updates
end

function with_discrete_outputs(m, wrapper)
    haskey(m, :discrete_outputs) || return (;)
    values = get_values!(zeros(length(m.discrete_outputs)), wrapper, wrapper.discrete_outputs)
    return (; discrete_outputs = values)
end

"""
    $(TYPEDSIGNATURES)

The positions of the states the FMI standard lets the caller set before a CoSimulation FMU
is initialized: those with a `start` value (`initial` exact or approx). The others are
calculated by the FMU (`initial = "calculated"`) and must not be written.
"""
settable_states(wrapper) = findall(
    vr -> any(
        mv -> fmi_kind(wrapper.fmu, mv) === real_kind(wrapper.fmu) && FMI.getStartValue(mv) !== nothing,
        FMI.modelVariablesForValueReference(wrapper.fmu, vr)
    ), wrapper.state_value_references
)

"""
    $(TYPEDEF)

A callable parameter giving the FMU's own initial continuous states: initializes a fresh
instance with the given inputs and parameters and reads its states back. `FMIComponent`
uses it as the default initial value of each state, so a state with no `start` value gets
the value the FMU calculates, and one whose start depends on parameters follows them.
"""
struct FMUInitialStates
    """
    `:ME` or `:CS`, which decides how the instance is created.
    """
    type::Symbol
end

function (fn::FMUInitialStates)(wrapper, inputs, params, t)
    inputs, params = as_fmi_vector(inputs), as_fmi_vector(params)
    out = zeros(length(wrapper.state_value_references))
    if fn.type == :ME
        # a running solve owns the instance; between solves start from a fresh one
        wrapper.in_solve || reset_instance!(wrapper)
        get_instance_ME!(wrapper, inputs, params, Float64(t))
        get_me_states!(wrapper, out)
    else
        # the CoSimulation instance is created again, with these values, when the solve
        # starts; this one only reports the initial states
        reset_instance!(wrapper)
        get_instance_CS!(wrapper, Float64[], inputs, params, Float64(t), nothing)
        get_cs_states!(wrapper, out)
        reset_instance!(wrapper)
    end
    return out
end

get_me_states!(wrapper::FMI2InstanceWrapper, out) = @statuscheck FMI.fmi2GetContinuousStates!(
    wrapper.instance, out, Csize_t(length(out))
)
get_me_states!(wrapper::FMI3InstanceWrapper, out) = @statuscheck FMI.fmi3GetContinuousStates!(
    wrapper.instance, out, Csize_t(length(out))
)

get_cs_states!(wrapper::FMI2InstanceWrapper, out) = @statuscheck FMI.fmi2GetReal!(
    wrapper.instance, wrapper.state_value_references, Csize_t(length(out)), out
)
get_cs_states!(wrapper::FMI3InstanceWrapper, out) = @statuscheck FMI.fmi3GetFloat64!(
    wrapper.instance, wrapper.state_value_references, Csize_t(length(out)), out, Csize_t(length(out))
)

"""
    $(TYPEDSIGNATURES)

Between solves, free a Model Exchange instance created with other parameter values than
`params`. Reading a solution after a solve evaluates the FMU, which creates an instance
with that solution's parameters; a later solve of a remade problem would otherwise reuse
it and silently run with the old values. During a solve the instance is kept: a parameter
the solve changes must not discard the FMU's internal state.
"""
function replace_stale_instance!(wrapper, params)
    if wrapper.instance !== nothing && !wrapper.in_solve &&
            !isequal(wrapper.instance_params, params)
        reset_instance!(wrapper)
    end
    return nothing
end

"""
    $(TYPEDSIGNATURES)

An affect function for use inside an `ImperativeAffect`, run when a solve of a Model
Exchange FMU starts, resetting the instance and reading initial discrete outputs.
`o` must contain the wrapper, inputs, and parameters.
"""
function fmiMEInitialize!(m, o, ctx, integrator)
    wrapper = o.wrapper
    reset_instance!(wrapper)
    wrapper.in_solve = true
    wrapper.last_completed_step = Float64(integrator.t)
    get_instance_ME!(wrapper, as_fmi_vector(o.inputs), as_fmi_vector(o.params), wrapper.last_completed_step)
    return with_discrete_outputs(m, wrapper)
end

"""
    $(TYPEDSIGNATURES)

The Model Exchange finalizer: marks the solve as over and frees the instance.
"""
function fmiMEFinalize!(m, o, ctx, integrator)
    wrapper = o.wrapper
    wrapper.in_solve = false
    reset_instance!(wrapper)
    return (;)
end

"""
    $(TYPEDSIGNATURES)

`v` as the `Vector{Float64}` FMIBase requires of the values it writes to an FMU. Solvers
can hand other vectors to the registered functions, such as the fixed-size `SVector`s used
for small systems, which FMIBase refuses.
"""
as_fmi_vector(v::Vector{Float64}) = v
as_fmi_vector(v) = Vector{Float64}(v)

function set_me_point!(wrapper::FMI2InstanceWrapper, states, t)
    @statuscheck FMI.fmi2SetTime(wrapper.instance, Float64(t))
    isempty(states) || @statuscheck FMI.fmi2SetContinuousStates(wrapper.instance, states)
    return
end

function set_me_point!(wrapper::FMI3InstanceWrapper, states, t)
    @statuscheck FMI.fmi3SetTime(wrapper.instance, Float64(t))
    isempty(states) || @statuscheck FMI.fmi3SetContinuousStates(wrapper.instance, states)
    return
end

function sync_me_point!(wrapper, states, inputs, t)
    set_me_point!(wrapper, states, t)
    split = wrapper.input_split
    continuous = isempty(split.discrete_positions) ? inputs : view(inputs, split.continuous_positions)
    set_values!(wrapper, split.continuous, continuous)
    return
end

get_derivatives!(wrapper::FMI2InstanceWrapper, out) =
    @statuscheck FMI.fmi2GetDerivatives!(wrapper.instance, out)
get_derivatives!(wrapper::FMI3InstanceWrapper, out) =
    @statuscheck FMI.fmi3GetContinuousStateDerivatives!(wrapper.instance, out)

"""
    $(TYPEDSIGNATURES)

Update the internal state of the ME FMU and return a vector of updated values
for continuous state derivatives and output variables respectively. Needs to be a
callable struct to enable symbolic registration with an inferred return size.
"""
function (wrapper::Union{FMI2InstanceWrapper, FMI3InstanceWrapper})(
        states, inputs, params, t
    )
    states, inputs, params = as_fmi_vector(states), as_fmi_vector(inputs), as_fmi_vector(params)
    get_instance_ME!(wrapper, inputs, params, t)

    states_buffer = wrapper.states_buffer
    outputs_buffer = wrapper.outputs_buffer
    @assert length(states_buffer) == length(states)
    sync_me_point!(wrapper, states, inputs, t)
    # OpenModelica FMUs can clear their update flag when computing derivatives alone.
    get_values!(outputs_buffer, wrapper, wrapper.outputs)
    get_values!(wrapper.discrete_outputs_buffer, wrapper, wrapper.discrete_outputs)
    isempty(states_buffer) || get_derivatives!(wrapper, states_buffer)
    get_event_indicators!(wrapper)
    ns, no = length(states_buffer), length(outputs_buffer)
    nd, ne = length(wrapper.discrete_outputs_buffer), length(wrapper.event_indicators)
    wrapper.res_buffer[1:ns] .= states_buffer
    wrapper.res_buffer[(ns + 1):(ns + no)] .= outputs_buffer
    wrapper.res_buffer[(ns + no + 1):(ns + no + nd)] .= wrapper.discrete_outputs_buffer
    wrapper.res_buffer[(ns + no + nd + 1):(ns + no + nd + ne)] .= wrapper.event_indicators
    # crosses zero at the next time event; stays negative while none is scheduled
    wrapper.res_buffer[end] = isnan(wrapper.next_event_time) ? -1.0 : t - wrapper.next_event_time
    return wrapper.res_buffer
end

"""
    $(TYPEDSIGNATURES)

An affect function for use inside an `ImperativeAffect`. This should be triggered at the
end of the solve, regardless of whether it succeeded or failed. Expects `p` to be a
1-length array containing the index of the instance wrapper (`FMI2InstanceWrapper` or
`FMI3InstanceWrapper`) in the parameter object.
"""
function fmiFinalize!(m, o, ctx, integrator)
    wrapper = o.wrapper
    reset_instance!(wrapper)
    return (;)
end

"""
    $(TYPEDEF)

The condition of a CoSimulation FMU's last step, as a callable parameter: 1.0 once the solve
has reached its final time and the FMU has not been stepped there. The periodic callback
that steps the FMU does not fire at the final time, so without this the solution's last
values would be one communication step behind (or short of a final time that is not a
multiple of the step).
"""
struct CSFinalStepDue end

function (::CSFinalStepDue)(functor, t)
    return (t == functor.solve_end[] && t > functor.t_last[]) ? 1.0 : 0.0
end

variable_communication_step(wrapper) =
    wrapper.fmu.modelDescription.coSimulation.canHandleVariableCommunicationStepSize === true

function initialize_cs_times!(ctx, o, integrator)
    start, finish = integrator.sol.prob.tspan
    experiment_end = isnan(ctx.stop_time) ? finish : ctx.stop_time
    experiment_end >= finish || throw(ArgumentError("The FMU experiment stop_time must be at or after the solve endpoint."))
    if !variable_communication_step(o.wrapper)
        duration = finish - start
        steps = round(duration / o.dt)
        grid_duration = steps * o.dt
        if duration < 0 || (
                iszero(steps) ? !iszero(duration) :
                    !(start + grid_duration == finish || isapprox(duration, grid_duration; rtol = 8eps(Float64), atol = 0))
            )
            throw(ArgumentError("The solve endpoint must align with communication_step_size for an FMU that requires fixed communication steps."))
        end
    end
    ctx.solve_end[] = finish
    ctx.t_last[] = o.t
    return experiment_end
end

"""
    $(TYPEDEF)

Provides Co-Simulation output values to the initialization problem and tracks the
communication point and experiment limit during a solve.
"""
struct CSFunctor
    output_count::Int
    solve_end::Base.RefValue{Float64}
    t_last::Base.RefValue{Float64}
    stop_time::Float64
end

CSFunctor(output_count, stop_time) = CSFunctor(output_count, Ref(NaN), Ref(NaN), stop_time)

function (fn::CSFunctor)(wrapper, states, inputs, params, t)
    stop_time = isnan(fn.stop_time) ? nothing : fn.stop_time
    get_instance_CS!(wrapper, as_fmi_vector(states), as_fmi_vector(inputs), as_fmi_vector(params), t, stop_time)
    return get_values!(zeros(eltype(states), fn.output_count), wrapper, wrapper.outputs)
end

@register_array_symbolic (fn::CSFunctor)(
    wrapper::Union{FMI2InstanceWrapper, FMI3InstanceWrapper}, states::Vector{<:Real},
    inputs::Vector{<:Real}, params::Vector{<:Real}, t::Real
) begin
    size = (fn.output_count,)
    eltype = eltype(states)
    ndims = 1
end

function fmiCSInitialize!(m, o, ctx, integrator)
    wrapper = o.wrapper
    reset_instance!(wrapper)
    experiment_end = initialize_cs_times!(ctx, o, integrator)
    states = haskey(m, :states) ? as_fmi_vector(m.states) : Float64[]
    get_instance_CS!(wrapper, states, as_fmi_vector(o.inputs), as_fmi_vector(o.params), o.t, experiment_end)
    return read_cs_values!(m, wrapper)
end

function read_cs_values!(m, wrapper)
    haskey(m, :states) && get_cs_states!(wrapper, m.states)
    haskey(m, :outputs) && get_values!(m.outputs, wrapper, wrapper.outputs)
    return m
end

function fmiCSStep!(m, o, ctx, integrator)
    from = ctx.t_last[]
    to = Float64(integrator.t)
    to > from || return (;)
    wrapper = o.wrapper
    h = variable_communication_step(wrapper) ? to - from : o.dt
    set_values!(wrapper, wrapper.inputs, o.inputs)
    do_cs_step!(wrapper, from, h)
    ctx.t_last[] = to
    return read_cs_values!(m, wrapper)
end

function do_cs_step!(wrapper::FMI2InstanceWrapper, from, h)
    @statuscheck FMI.fmi2DoStep(wrapper.instance, from, h, FMI.fmi2True)
    return
end

function do_cs_step!(wrapper::FMI3InstanceWrapper, from, h)
    event_encountered = Ref(FMI.fmi3False)
    terminate_simulation = Ref(FMI.fmi3False)
    early_return = Ref(FMI.fmi3False)
    last_successful_time = Ref(zero(FMI.fmi3Float64))
    @statuscheck FMI.fmi3DoStep!(
        wrapper.instance, from, h, FMI.fmi3True, event_encountered,
        terminate_simulation, early_return, last_successful_time
    )
    @assert event_encountered[] == FMI.fmi3False
    @assert terminate_simulation[] == FMI.fmi3False
    @assert early_return[] == FMI.fmi3False
    return
end

end
