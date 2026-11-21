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
    :Int64 => (FMI.fmi3Int64, FMI.fmi3SetInt64, FMI.fmi3GetInt64!),
    :UInt8 => (FMI.fmi3UInt8, FMI.fmi3SetUInt8, FMI.fmi3GetUInt8!),
    :UInt16 => (FMI.fmi3UInt16, FMI.fmi3SetUInt16, FMI.fmi3GetUInt16!),
    :UInt32 => (FMI.fmi3UInt32, FMI.fmi3SetUInt32, FMI.fmi3GetUInt32!),
    :UInt64 => (FMI.fmi3UInt64, FMI.fmi3SetUInt64, FMI.fmi3GetUInt64!),
    :Boolean => (FMI.fmi3Boolean, FMI.fmi3SetBoolean, FMI.fmi3GetBoolean!),
    :Enumeration => (FMI.fmi3Int64, FMI.fmi3SetInt64, FMI.fmi3GetInt64!),
)

"""
    $(TYPEDSIGNATURES)

The FMI type the FMU declares for `vr`, as a key of `FMI2_ACCESSORS`, or `nothing` when
FMIComponent cannot carry it as a number (a string).
"""
function fmi_kind(fmu::FMI.FMU2, vr)
    mv = FMI.modelVariablesForValueReference(fmu.modelDescription, vr)[1]
    mv.Real === nothing || return :Real
    (mv.Integer === nothing && mv.Enumeration === nothing) || return :Integer
    mv.Boolean === nothing || return :Boolean
    return nothing
end

"""
    $(TYPEDSIGNATURES)

The FMI type the FMU declares for `vr`, as a key of `FMI3_ACCESSORS`, or `nothing` when
FMIComponent cannot carry it as a number (a string, binary or clock).
"""
function fmi_kind(fmu::FMI.FMU3, vr)
    mv = FMI.modelVariablesForValueReference(fmu.modelDescription, vr)[1]
    kind = Symbol(chopprefix(string(nameof(typeof(mv))), "fmi3Variable"))
    return haskey(FMI3_ACCESSORS, kind) ? kind : nothing
end

real_kind(::FMI.FMU2) = :Real
real_kind(::FMI.FMU3) = :Float64

"""
    $(TYPEDEF)

Value references in the order the caller passes or expects their values, grouped by the FMI
type the FMU declares for each. ModelingToolkit carries every FMU variable as a `Float64`; the
groups say which typed FMI function reads or writes each one, so the conversion happens only
at the FMU boundary.

# Fields

$(TYPEDFIELDS)
"""
struct TypedValueReferences{V}
    """
    Every value reference, in the caller's order.
    """
    all::Vector{V}
    """
    Whether every variable has the FMI version's 64-bit float type, so that `all` can be
    read or written in one call without conversion.
    """
    all_real::Bool
    """
    For each FMI type present: its key in the accessor table, the positions of its variables
    in `all`, and their value references.
    """
    groups::Vector{Tuple{Symbol, Vector{Int}, Vector{V}}}
end

function TypedValueReferences(fmu, vrs::Vector{V}) where {V}
    kinds = [fmi_kind(fmu, vr) for vr in vrs]
    groups = Tuple{Symbol, Vector{Int}, Vector{V}}[]
    for kind in unique(kinds)
        positions = findall(==(kind), kinds)
        push!(groups, (kind, positions, vrs[positions]))
    end
    return TypedValueReferences{V}(vrs, all(==(real_kind(fmu)), kinds), groups)
end

to_fmi(::Type{T}, x) where {T <: AbstractFloat} = T(x)
to_fmi(::Type{T}, x) where {T <: Integer} = round(T, x)

"""
    $(TYPEDSIGNATURES)

Whether `vr` changes only at events: every non-float variable, and a float one the FMU
declares `discrete`. A Model Exchange FMU accepts a new value for such an input only in
Event Mode.
"""
function is_discrete(fmu::FMI.FMU2, vr)
    mv = FMI.modelVariablesForValueReference(fmu.modelDescription, vr)[1]
    return fmi_kind(fmu, vr) !== :Real || mv.variability == FMI.fmi2VariabilityDiscrete
end
function is_discrete(fmu::FMI.FMU3, vr)
    mv = FMI.modelVariablesForValueReference(fmu.modelDescription, vr)[1]
    return !(fmi_kind(fmu, vr) in (:Float64, :Float32)) ||
        mv.variability == FMI.fmi3VariabilityDiscrete
end

"""
    $(TYPEDEF)

The inputs of a Model Exchange FMU, split into those that may change continuously and
those that change only at events.

# Fields

$(TYPEDFIELDS)
"""
struct InputSplit{V}
    """
    Positions of the continuous inputs in the full input vector.
    """
    continuous_positions::Vector{Int}
    continuous::TypedValueReferences{V}
    """
    Positions of the discrete inputs in the full input vector.
    """
    discrete_positions::Vector{Int}
    discrete::TypedValueReferences{V}
    """
    The full input vector as last written to the FMU, to tell when a discrete input changed.
    """
    last::Vector{Float64}
    """
    The full input vector as it was when the discrete outputs were last read, to tell when
    they must be read again.
    """
    reported::Vector{Float64}
end

function InputSplit(fmu, vrs::Vector{V}) where {V}
    discrete = [is_discrete(fmu, vr) for vr in vrs]
    cpos, dpos = findall(!, discrete), findall(discrete)
    return InputSplit{V}(
        cpos, TypedValueReferences(fmu, vrs[cpos]), dpos, TypedValueReferences(fmu, vrs[dpos]),
        fill(NaN, length(vrs)), fill(NaN, length(vrs))
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

Write `values` (in the order of `refs.all`) to `wrapper.instance`, one typed FMI call per
group.
"""
function set_values!(wrapper, refs::TypedValueReferences, values)
    version = fmi_version(wrapper)
    accessors = version isa Val{2} ? FMI2_ACCESSORS : FMI3_ACCESSORS
    for (kind, positions, vrs) in refs.groups
        T, set, _ = accessors[kind]
        converted = T[to_fmi(T, values[i]) for i in positions]
        n = Csize_t(length(vrs))
        status = version isa Val{2} ? set(wrapper.instance, vrs, n, converted) :
            set(wrapper.instance, vrs, n, converted, n)
        check_status!(wrapper, version, status, nameof(set))
    end
    return
end

"""
    $(TYPEDSIGNATURES)

Read the variables in `refs` from `wrapper.instance` into `out` (in the order of
`refs.all`), one typed FMI call per group. Returns `out`.
"""
function get_values!(out, wrapper, refs::TypedValueReferences)
    version = fmi_version(wrapper)
    accessors = version isa Val{2} ? FMI2_ACCESSORS : FMI3_ACCESSORS
    for (kind, positions, vrs) in refs.groups
        T, _, get! = accessors[kind]
        values = zeros(T, length(vrs))
        n = Csize_t(length(vrs))
        status = version isa Val{2} ? get!(wrapper.instance, vrs, n, values) :
            get!(wrapper.instance, vrs, n, values, n)
        check_status!(wrapper, version, status, nameof(get!))
        for (j, i) in enumerate(positions)
            out[i] = values[j]
        end
    end
    return out
end

@static if !hasmethod(FMI.getValueReferencesAndNames, Tuple{FMI.fmi3ModelDescription})
    """
        $(TYPEDSIGNATURES)

        This is type piracy, but FMI.jl is missing this implementation. It allows
        `FMI.getStateValueReferencesAndNames` to work.
    """
    function FMI.getValueReferencesAndNames(
            md::FMI.fmi3ModelDescription; vrs = md.valueReferences
        )
        dict = Dict{FMI.fmi3ValueReference, Array{String}}()
        for vr in vrs
            dict[vr] = FMI.valueReferenceToString(md, vr)
        end
        return dict
    end
end

"""
    $(TYPEDSIGNATURES)

A component that wraps an FMU loaded via FMI.jl. The FMI version (2 or 3) should be
provided as a `Val` to the function. Supports Model Exchange and CoSimulation FMUs.
Continuous states must be `FMI.fmi2Real` or `FMI.fmi3Float64`. Inputs, outputs and
parameters of any numeric FMI type (integers, enumerations and booleans included) are carried
as `Float64` variables and converted when they are passed to or read from the FMU; string,
binary and clock variables are left out with a warning. Does not support events in the FMU. Does not support automatic
differentiation. Parameters of the FMU will have defaults corresponding to their initial
values in the FMU specification. All other variables will not have a default. Hierarchical
names in the FMU of the form `namespace.variable` are transformed into symbolic variables
with the name `namespace__variable`.

# Keyword Arguments

- `fmu`: The FMU loaded via `FMI.loadFMU`.
- `tolerance`: The tolerance to provide to the FMU. Not used for v3 FMUs since it is not
  supported by FMI.jl.
- `communication_step_size`: The periodic interval at which communication with CoSimulation
  FMUs will occur. Must be provided for CoSimulation FMU components.
- `reinitializealg`: The DAE initialization algorithm to use for the callback managing the
  FMU. For CoSimulation FMUs whose states/outputs are used in algebraic equations of the
  system, this needs to be an algorithm that will solve for the new algebraic variables.
  For example, `OrdinaryDiffEqCore.BrownFullBasicInit()`.
- `stop_time`: The `stopTime` for the FMU, as defined in the FMI spec. By default, this will
  be set to `tspan[2]` when the `ODEProblem` is solved. An explicit value will overwrite this.
- `type`: Either `:ME` or `:CS` depending on whether `fmu` is a Model Exchange or
  CoSimulation FMU respectively.
- `name`: The name of the system.
"""
function MTK.FMIComponent(
        ::Val{Ver}; fmu = nothing, tolerance = 1.0e-6,
        communication_step_size = nothing, reinitializealg = nothing,
        stop_time = nothing, type, name
    ) where {Ver}
    if Ver != 2 && Ver != 3
        throw(ArgumentError("FMI Version must be `2` or `3`"))
    end
    if type == :CS && communication_step_size === nothing
        throw(ArgumentError("`communication_step_size` must be specified for Co-Simulation FMUs."))
    end
    # mapping from MTK variable to value reference
    value_references = Dict()
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
    state_value_references = fmu.modelDescription.stateValueReferences
    output_value_reference_set = Set(keys(FMI.getOutputValueReferencesAndNames(fmu)))
    state_names = FMI.getStateValueReferencesAndNames(fmu)
    fmi_variables_to_mtk_variables!(
        fmu, [vr => state_names[vr] for vr in state_value_references],
        value_references, diffvars, states, observed;
        postprocess_variable = (v, vr) -> vr in output_value_reference_set ?
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

    # parse derivatives of states
    # the variables passed to `postprocess_variable` haven't been differentiated yet, so they
    # should match one variable in states. That's the one this is the derivative of, and we
    # keep track of this ordering
    derivative_order = []
    function derivative_order_postprocess(var, _)
        idx = findfirst(isequal(var), states)
        idx === nothing || push!(derivative_order, states[idx])
        return var
    end
    fmi_variables_to_mtk_variables!(
        fmu, FMI.getDerivateValueReferencesAndNames(fmu), value_references, dervars,
        states, der_observed; postprocess_variable = derivative_order_postprocess
    )
    @assert length(derivative_order) == length(dervars)

    # parse the inputs to the FMU
    inputs = []
    fmi_variables_to_mtk_variables!(
        fmu, FMI.getInputValueReferencesAndNames(fmu),
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
        fmu, filter(p -> !(first(p) in state_value_references), FMI.getOutputValueReferencesAndNames(fmu)),
        value_references, outputs, states, observed; postprocess_variable = (v, _) -> MTKBase.setoutput(
            v, true
        )
    )
    # A Model Exchange FMU's discrete outputs change only at events, which can depend on
    # discrete state inside the FMU that the system cannot see. They are kept as discrete
    # variables, set when the system starts and after every FMU event, so a solution
    # reports what the FMU said during the solve rather than asking a fresh instance.
    discrete_outputs = []
    if type == :ME
        discrete_outputs = filter(v -> is_discrete(fmu, value_references[v]), outputs)
        outputs = filter(v -> !is_discrete(fmu, value_references[v]), outputs)
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
        fmu, FMI.getParameterValueReferencesAndNames(fmu), value_references,
        params, [], parameter_dependencies, defs; parameters = true
    )
    # create a symbolic variable for the parameter buffer
    __mtk_internal_p = copy(params)
    if isempty(__mtk_internal_p)
        __mtk_internal_p = Float64[]
    end

    derivative_value_references = UInt32[value_references[var] for var in dervars]
    state_value_references = UInt32[value_references[var] for var in diffvars]
    output_value_references = UInt32[value_references[var] for var in outputs]
    discrete_output_value_references = UInt32[value_references[var] for var in discrete_outputs]
    input_value_references = UInt32[value_references[var] for var in inputs]
    param_value_references = UInt32[value_references[var] for var in params]

    # create a parameter for the instance wrapper
    # this manages the creation and deallocation of FMU instances
    buffer_length = length(diffvars) + length(outputs) + length(discrete_outputs) +
        number_of_events(fmu)
    if Ver == 2
        @parameters (wrapper::FMI2InstanceWrapper)(..)[1:buffer_length] = FMI2InstanceWrapper(
            fmu, derivative_value_references, state_value_references, output_value_references,
            param_value_references, input_value_references, tolerance,
            discrete_output_value_references
        )
    else
        @parameters (wrapper::FMI3InstanceWrapper)(..)[1:buffer_length] = FMI3InstanceWrapper(
            fmu, derivative_value_references, state_value_references,
            output_value_references, param_value_references, input_value_references,
            discrete_output_value_references
        )
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
        event_affect = MTKBase.ImperativeAffect(
            fmiMEEvent!; observed = (; wrapper),
            modified = (;
                (isempty(diffvars) ? (;) : (; states = diffvars))..., discrete_modified...,
            )
        )
        # the FMI standard has the environment step past the zero crossing, so the FMU sees
        # the indicator's new sign when it enters Event Mode
        push!(continuous_events, MTKBase.SymbolicContinuousCallback(
            [var ~ 0 for var in event_vars], event_affect;
            rootfind = SciMLBase.RightRootFind
        ))
        # a changed discrete input can change the discrete outputs with no FMU event; this
        # reads them again at the end of the step where the input changed
        if !isempty(discrete_outputs) && any(vr -> is_discrete(fmu, vr), input_value_references)
            @parameters (input_check::DiscreteInputCheck)(..) = DiscreteInputCheck()
            push!(params, input_check)
            push!(discrete_events, MTKBase.SymbolicDiscreteCallback(
                input_check(wrapper, __mtk_internal_x) > 0.5,
                MTKBase.ImperativeAffect(
                    fmiMEInputEvent!; observed = (; wrapper, inputs = __mtk_internal_x),
                    modified = discrete_modified
                )
            ))
        end
        isempty(discrete_outputs) || push!(params, __mtk_internal_do)

        # instance management callback which deallocates the instance when
        # necessary and notifies the FMU of completed integrator steps
        initialize_affect = MTKBase.ImperativeAffect(fmiMEInitialize!; observed = (; wrapper))
        finalize_affect = MTKBase.ImperativeAffect(fmiMEFinalize!; observed = (; wrapper))
        step_affect = MTKBase.ImperativeAffect(Returns((;)))
        instance_management_callback = MTKBase.SymbolicDiscreteCallback(
            (t == t - 1), step_affect; initialize = initialize_affect,
            finalize = finalize_affect, reinitializealg = SciMLBase.NoInit()
        )

        push!(params, wrapper)
        append!(observed, der_observed)
    elseif type == :CS
        _functor = if Ver == 2
            FMI2CSFunctor(state_value_references, output_value_references, Base.Ref{FMI.fmi2Real}(something(stop_time, NaN)))
        else
            FMI3CSFunctor(state_value_references, output_value_references, Base.Ref{FMI.fmi3Float64}(something(stop_time, NaN)))
        end
        @parameters (functor::(typeof(_functor)))(..)[1:(length(__mtk_internal_u) + length(__mtk_internal_o))] = _functor
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
        push!(discrete_events, MTKBase.SymbolicDiscreteCallback(
            final_step_due(functor, t) > 0.5,
            MTKBase.ImperativeAffect(
                fmiCSFinalStep!; observed = cb_observed, modified = cb_modified, ctx = _functor
            ); reinitializealg
        ))

        # guarded in case there are no outputs/states and the variable is `[]`.
        symbolic_type(__mtk_internal_o) == NotSymbolic() || push!(params, __mtk_internal_o)
        symbolic_type(__mtk_internal_u) == NotSymbolic() || push!(params, __mtk_internal_u)

        push!(params, wrapper, functor)
    end

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

A utility function which accepts an FMU `fmu` and a mapping from value reference to a
list of associated names `varmap`. A symbolic variable is created for each name. The
associated value reference is kept track of in `value_references`. In case there are
multiple names for a value reference, the symbolic variable for the first name is pushed
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
        fmu::Union{FMI.FMU2, FMI.FMU3}, varmap,
        value_references::AbstractDict, truevars, allvars,
        obseqs, defs = Dict(); parameters = false, postprocess_variable = (v, _) -> v
    )
    for (valRef, varnames) in varmap
        # every numeric FMU variable is carried as a `Float64`; `set_values!` and
        # `get_values!` convert to and from the FMU's declared type
        if fmi_kind(fmu, valRef) === nothing
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
        defval = FMI.getStartValue(fmu, valRef)
        defs[vars[1]] = defval isa Number ? Float64(defval) : defval
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
    The value references for derivatives of states of the FMU, in the order that the
    caller expects them to be returned when calling this struct.
    """
    const derivative_value_references::Vector{FMI.fmi2ValueReference}
    """
    The value references for the states of the FMU.
    """
    const state_value_references::Vector{FMI.fmi2ValueReference}
    """
    The value references for outputs of the FMU, in the order that the caller expects
    them to be returned when calling this struct.
    """
    const output_value_references::Vector{FMI.fmi2ValueReference}
    """
    The parameter value references. These should be in the same order as the parameter
    vector passed to functions involving this wrapper.
    """
    const param_value_references::Vector{FMI.fmi2ValueReference}
    """
    The input value references. These should be in the same order as the inputs passed
    to functions involving this wrapper.
    """
    const input_value_references::Vector{FMI.fmi2ValueReference}
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
    `output_value_references`, `param_value_references` and `input_value_references`
    grouped by the FMI type of each variable.
    """
    const outputs::TypedValueReferences{FMI.fmi2ValueReference}
    const params::TypedValueReferences{FMI.fmi2ValueReference}
    const inputs::TypedValueReferences{FMI.fmi2ValueReference}
    const input_split::InputSplit{FMI.fmi2ValueReference}
    """
    Buffer for the FMU's event indicators (Model Exchange).
    """
    const event_indicators::Vector{FMI.fmi2Real}
    """
    The discrete outputs of a Model Exchange FMU, grouped by FMI type, and a buffer for
    their values. They change only at events and are kept by the system as discrete
    variables, so that a solution reports the values the FMU had during the solve.
    """
    const discrete_outputs::TypedValueReferences{FMI.fmi2ValueReference}
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
end

"""
    $(TYPEDSIGNATURES)

Create an `FMI2InstanceWrapper` with no instance.
"""
function FMI2InstanceWrapper(
        fmu, ders, states, outputs, params, inputs, tolerance,
        discrete_outputs = FMI.fmi2ValueReference[]
    )
    n = length(states) + length(outputs) + length(discrete_outputs) + number_of_events(fmu)
    return FMI2InstanceWrapper(
        fmu, ders, states, outputs, params, inputs, tolerance, nothing,
        zeros(FMI.fmi2Real, length(states)), zeros(FMI.fmi2Real, length(outputs)),
        zeros(FMI.fmi2Real, n),
        TypedValueReferences(fmu, outputs), TypedValueReferences(fmu, params),
        TypedValueReferences(fmu, inputs), InputSplit(fmu, inputs),
        zeros(FMI.fmi2Real, number_of_events(fmu) - 1),
        TypedValueReferences(fmu, discrete_outputs), zeros(length(discrete_outputs)), NaN,
        zeros(length(params)), false
    )
end

Base.nameof(::FMI2InstanceWrapper) = :FMI2InstanceWrapper
fmi_version(::FMI2InstanceWrapper) = Val(2)

"""
    $(TYPEDSIGNATURES)

Common functionality for creating an instance of a v2 FMU. Does not check if
`wrapper.instance` is already present, and overwrites the existing value with
a new instance. `inputs` should be in the order of `wrapper.input_value_references`.
`params` should be in the order of `wrapper.param_value_references`. `t` is the current
time. Returns the created instance, which is also stored in `wrapper.instance`.
"""
function get_instance_common!(wrapper::FMI2InstanceWrapper, inputs, params, t, t_end = nothing)
    wrapper.instance = FMI.fmi2Instantiate!(wrapper.fmu)::FMI.FMU2Component
    wrapper.instance_params .= params
    set_values!(wrapper, wrapper.inputs, inputs)
    wrapper.input_split.last .= inputs
    wrapper.input_split.reported .= inputs
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
        if !isempty(states)
            @statuscheck FMI.fmi2SetReal(
                wrapper.instance, wrapper.state_value_references,
                Csize_t(length(wrapper.state_value_references)), states
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

"""
    $(TYPEDSIGNATURES)

Call `fmiXCompletedIntegratorStep` with `noSetFMUStatePriorToCurrentPoint` as false.
"""
function partiallyCompleteIntegratorStep(wrapper::FMI2InstanceWrapper)
    return @statuscheck FMI.fmi2CompletedIntegratorStep(wrapper.instance, FMI.fmi2False)
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
    The value references for derivatives of states of the FMU, in the order that the
    caller expects them to be returned when calling this struct.
    """
    const derivative_value_references::Vector{FMI.fmi3ValueReference}
    const state_value_references::Vector{FMI.fmi3ValueReference}
    """
    The value references for outputs of the FMU, in the order that the caller expects
    them to be returned when calling this struct.
    """
    const output_value_references::Vector{FMI.fmi3ValueReference}
    """
    The parameter value references. These should be in the same order as the parameter
    vector passed to functions involving this wrapper.
    """
    const param_value_references::Vector{FMI.fmi3ValueReference}
    """
    The input value references. These should be in the same order as the inputs passed
    to functions involving this wrapper.
    """
    const input_value_references::Vector{FMI.fmi3ValueReference}
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
    `output_value_references`, `param_value_references` and `input_value_references`
    grouped by the FMI type of each variable.
    """
    const outputs::TypedValueReferences{FMI.fmi3ValueReference}
    const params::TypedValueReferences{FMI.fmi3ValueReference}
    const inputs::TypedValueReferences{FMI.fmi3ValueReference}
    const input_split::InputSplit{FMI.fmi3ValueReference}
    """
    Buffer for the FMU's event indicators (Model Exchange).
    """
    const event_indicators::Vector{FMI.fmi3Float64}
    """
    The discrete outputs of a Model Exchange FMU, grouped by FMI type, and a buffer for
    their values. They change only at events and are kept by the system as discrete
    variables, so that a solution reports the values the FMU had during the solve.
    """
    const discrete_outputs::TypedValueReferences{FMI.fmi3ValueReference}
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
end

"""
    $(TYPEDSIGNATURES)

Create an `FMI3InstanceWrapper` with no instance.
"""
function FMI3InstanceWrapper(
        fmu, ders, states, outputs, params, inputs, discrete_outputs = FMI.fmi3ValueReference[]
    )
    n = length(states) + length(outputs) + length(discrete_outputs) + number_of_events(fmu)
    return FMI3InstanceWrapper(
        fmu, ders, states, outputs, params, inputs, nothing,
        zeros(FMI.fmi3Float64, length(states)), zeros(FMI.fmi3Float64, length(outputs)),
        zeros(FMI.fmi3Float64, n),
        TypedValueReferences(fmu, outputs), TypedValueReferences(fmu, params),
        TypedValueReferences(fmu, inputs), InputSplit(fmu, inputs),
        zeros(FMI.fmi3Float64, number_of_events(fmu) - 1),
        TypedValueReferences(fmu, discrete_outputs), zeros(length(discrete_outputs)), NaN,
        zeros(length(params)), false
    )
end

Base.nameof(::FMI3InstanceWrapper) = :FMI3InstanceWrapper
fmi_version(::FMI3InstanceWrapper) = Val(3)

"""
    $(TYPEDSIGNATURES)

Common functionality for creating an instance of a v3 FMU. Since v3 FMUs need to be
instantiated differently depending on the type, this assumes `wrapper.instance` is a
freshly instantiated FMU which needs to be initialized. `inputs` should be in the order
of `wrapper.input_value_references`. `params` should be in the order of
`wrapper.param_value_references`. `t` is the current time. Returns `wrapper.instance`.
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
    wrapper.input_split.reported .= inputs

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
        if !isempty(states)
            @statuscheck FMI.fmi3SetFloat64(
                wrapper.instance, wrapper.state_value_references, states
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
function partiallyCompleteIntegratorStep(wrapper::FMI3InstanceWrapper)
    enterEventMode = Ref(FMI.fmi3False)
    terminateSimulation = Ref(FMI.fmi3False)
    @statuscheck FMI.fmi3CompletedIntegratorStep!(
        wrapper.instance, FMI.fmi3False, enterEventMode, terminateSimulation
    )
    @assert enterEventMode[] == FMI.fmi3False
    return @assert terminateSimulation[] == FMI.fmi3False
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
    wrapper::FMI3InstanceWrapper, states::Vector{<:Real},
    inputs::Vector{<:Real}, params::Vector{<:Real}, t::Real
) begin
    size = (length(fn.res_buffer),)
    eltype = eltype(states)
    ndims = 1
end

"""
    $(TYPEDSIGNATURES)

Move a Model Exchange instance from Continuous-Time Mode into Event Mode.
"""
function enter_event_mode!(wrapper::FMI2InstanceWrapper)
    @statuscheck FMI.fmi2EnterEventMode(wrapper.instance)
end
function enter_event_mode!(wrapper::FMI3InstanceWrapper)
    @statuscheck FMI.fmi3EnterEventMode(
        wrapper.instance, FMI.fmi3False, FMI.fmi3False, FMI.fmi3Int32[], Csize_t(0), FMI.fmi3False
    )
end

"""
    $(TYPEDSIGNATURES)

Iterate the FMU's discrete states until they settle (the event iteration of the FMI
standard), then return to Continuous-Time Mode. Returns whether the FMU changed the values
of its continuous states.
"""
function leave_event_mode!(wrapper::FMI2InstanceWrapper)
    info = FMI.fmi2EventInfo()
    info.newDiscreteStatesNeeded = FMI.fmi2True
    changed = false
    while info.newDiscreteStatesNeeded == FMI.fmi2True
        @statuscheck FMI.fmi2NewDiscreteStates!(wrapper.instance, info)
        info.terminateSimulation == FMI.fmi2True && error("The FMU asked to terminate the simulation.")
        changed |= info.valuesOfContinuousStatesChanged == FMI.fmi2True
    end
    wrapper.next_event_time = info.nextEventTimeDefined == FMI.fmi2True ? info.nextEventTime : NaN
    @statuscheck FMI.fmi2EnterContinuousTimeMode(wrapper.instance)
    return changed
end
function leave_event_mode!(wrapper::FMI3InstanceWrapper)
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
        terminate[] == FMI.fmi3True && error("The FMU asked to terminate the simulation.")
        changed |= values[] == FMI.fmi3True
    end
    wrapper.next_event_time = next_defined[] == FMI.fmi3True ? next_time[] : NaN
    @statuscheck FMI.fmi3EnterContinuousTimeMode(wrapper.instance)
    return changed
end

"""
    $(TYPEDSIGNATURES)

Write the discrete inputs of a Model Exchange instance if any differs from the value last
written. The FMI standard lets them change only in Event Mode, so this enters it, writes
them, iterates and returns to Continuous-Time Mode.
"""
function update_discrete_inputs!(wrapper, inputs)
    split = wrapper.input_split
    any(i -> !isequal(inputs[i], split.last[i]), split.discrete_positions) || return
    enter_event_mode!(wrapper)
    set_values!(wrapper, split.discrete, view(inputs, split.discrete_positions))
    leave_event_mode!(wrapper)
    split.last .= inputs
    return
end

"""
    $(TYPEDSIGNATURES)

Read the FMU's event indicators into `wrapper.event_indicators`.
"""
function get_event_indicators!(wrapper::FMI2InstanceWrapper)
    n = length(wrapper.event_indicators)
    n == 0 && return
    @statuscheck FMI.fmi2GetEventIndicators!(wrapper.instance, wrapper.event_indicators, Csize_t(n))
end
function get_event_indicators!(wrapper::FMI3InstanceWrapper)
    n = length(wrapper.event_indicators)
    n == 0 && return
    @statuscheck FMI.fmi3GetEventIndicators!(wrapper.instance, wrapper.event_indicators, Csize_t(n))
end

"""
    $(TYPEDSIGNATURES)

An affect function for use inside an `ImperativeAffect`, triggered when one of a Model
Exchange FMU's event functions (an event indicator or the next time event) crosses zero.
Brings the FMU to the event's time and states, runs its event iteration, and writes back
any state the FMU reset (BouncingBall reverses its velocity here). `o` must contain the
`wrapper`; `m` contains the `states` if the FMU has any.
"""
function fmiMEEvent!(m, o, ctx, integrator)
    wrapper = o.wrapper
    wrapper.instance === nothing && return m
    states = isdefined(m, :states) ? as_fmi_vector(m.states) : Float64[]
    wrapper.instance(; x = states, t = Float64(integrator.t))
    enter_event_mode!(wrapper)
    changed = leave_event_mode!(wrapper)
    if changed && isdefined(m, :states)
        get_states!(wrapper, states)
        m = (; m..., states)
    end
    return with_discrete_outputs(m, wrapper)
end

"""
    $(TYPEDSIGNATURES)

`m` with its `discrete_outputs`, if it has them, read afresh from the FMU.
"""
function with_discrete_outputs(m, wrapper)
    isdefined(m, :discrete_outputs) || return m
    values = get_values!(zeros(length(m.discrete_outputs)), wrapper, wrapper.discrete_outputs)
    wrapper.input_split.reported .= wrapper.input_split.last
    return (; m..., discrete_outputs = values)
end

"""
    $(TYPEDEF)

A callable struct, used as a callable parameter so that it can take the instance wrapper:
returns `1.0` if a discrete input of a Model Exchange FMU differs from its value when the
FMU's discrete outputs were last read (so they must be read again), and `0.0` otherwise.
"""
struct DiscreteInputCheck end

function (::DiscreteInputCheck)(wrapper, inputs)
    split = wrapper.input_split
    changed = any(i -> !isequal(inputs[i], split.reported[i]), split.discrete_positions)
    return changed ? 1.0 : 0.0
end

"""
    $(TYPEDSIGNATURES)

An affect function for use inside an `ImperativeAffect`, triggered when a discrete input of
a Model Exchange FMU has changed. Writes it to the FMU (in Event Mode) if the derivative
evaluation has not already, then reads the discrete outputs. `o` must contain the
`wrapper` and the `inputs`; `m` contains the `discrete_outputs`.
"""
function fmiMEInputEvent!(m, o, ctx, integrator)
    wrapper = o.wrapper
    wrapper.instance === nothing && return m
    update_discrete_inputs!(wrapper, as_fmi_vector(o.inputs))
    return with_discrete_outputs(m, wrapper)
end

get_states!(wrapper::FMI2InstanceWrapper, out) = @statuscheck FMI.fmi2GetReal!(
    wrapper.instance, wrapper.state_value_references, Csize_t(length(out)), out
)
get_states!(wrapper::FMI3InstanceWrapper, out) = @statuscheck FMI.fmi3GetFloat64!(
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
Exchange FMU starts: frees any instance left from earlier evaluations (such as reading a
previous solution), so the solve starts from a freshly initialized FMU, and marks the
solve as running. `o` must contain the `wrapper`.
"""
function fmiMEInitialize!(m, o, ctx, integrator)
    wrapper = o.wrapper
    reset_instance!(wrapper)
    wrapper.in_solve = true
    return (;)
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
    instance = get_instance_ME!(wrapper, inputs, params, t)

    states_buffer = wrapper.states_buffer
    outputs_buffer = wrapper.outputs_buffer
    #buffer size matches
    @assert length(states_buffer) == length(states)
    # Defined in FMIBase.jl/src/eval.jl
    # Doesn't seem to be documented, but somehow this is the only way to
    # propagate inputs to the FMU consistently. I have no idea why.
    # FMIBase's evaluation reads and writes only 64-bit floats; any other input or output
    # goes through the typed calls instead. Discrete inputs are written only when they
    # change, in Event Mode.
    split = wrapper.input_split
    continuous_inputs = isempty(split.discrete_positions) ? inputs :
        inputs[split.continuous_positions]
    if split.continuous.all_real
        instance(;
            x = states, u = continuous_inputs, u_refs = split.continuous.all,
            p = params, p_refs = wrapper.param_value_references, t = t
        )
    else
        set_values!(wrapper, split.continuous, continuous_inputs)
        instance(; x = states, t = t)
    end
    update_discrete_inputs!(wrapper, inputs)
    # the spec requires completing the step before getting updated derivative/output values
    partiallyCompleteIntegratorStep(wrapper)
    if !wrapper.outputs.all_real
        instance(;
            dx = states_buffer, dx_refs = wrapper.derivative_value_references,
            y = similar(outputs_buffer, 0), y_refs = similar(wrapper.output_value_references, 0)
        )
        get_values!(outputs_buffer, wrapper, wrapper.outputs)
    else
        instance(;
            dx = states_buffer, dx_refs = wrapper.derivative_value_references,
            y = outputs_buffer, y_refs = wrapper.output_value_references
        )
    end
    get_values!(wrapper.discrete_outputs_buffer, wrapper, wrapper.discrete_outputs)
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
    from = functor.t_last[]
    tol = sqrt(eps(Float64)) * max(1.0, abs(t))
    return (!isnan(from) && t >= functor.t_end[] - tol && t - from > tol) ? 1.0 : 0.0
end

"""
    $(TYPEDSIGNATURES)

An affect function for use inside an `ImperativeAffect`, triggered by
[`CSFinalStepDue`](@ref): steps the CoSimulation FMU from where it was last stepped to the
final time. Has the same requirements for `m` and `o` as `fmiCSStep!`.
"""
function fmiCSFinalStep!(m, o, ctx, integrator)
    from = ctx.t_last[]
    return cs_step!(m, o, ctx, integrator, from, integrator.t - from)
end

"""
    $(TYPEDEF)

A callable struct useful for initializing v2 CoSimulation FMUs. When called, updates the
internal state of the FMU and gets updated values for output variables.

# Fields

$(TYPEDFIELDS)
"""
struct FMI2CSFunctor
    """
    The value references of state variables in the FMU.
    """
    state_value_references::Vector{FMI.fmi2ValueReference}
    """
    The value references of output variables in the FMU.
    """
    output_value_references::Vector{FMI.fmi2ValueReference}
    """
    `tspan[2]`, set during callback initialization
    """
    t_end::Base.RefValue{FMI.fmi2Real}
    """
    The time the FMU has been stepped to, so the step to the end of the solve can be
    taken when the periodic callback does not fire there.
    """
    t_last::Base.RefValue{FMI.fmi2Real}
    """
    The `stop_time` given to `FMIComponent`, or `NaN` to use the end of each solve's `tspan`.
    """
    stop_time::FMI.fmi2Real
end

FMI2CSFunctor(states, outputs, t_end) = FMI2CSFunctor(
    states, outputs, t_end, Base.Ref{FMI.fmi2Real}(NaN), t_end[]
)

function (fn::FMI2CSFunctor)(wrapper::FMI2InstanceWrapper, states, inputs, params, t)
    states = states isa SubArray ? copy(states) : states
    inputs = inputs isa SubArray ? copy(inputs) : inputs
    params = params isa SubArray ? copy(params) : params
    instance = get_instance_CS!(wrapper, states, inputs, params, t, fn.t_end[])
    if isempty(fn.output_value_references)
        return eltype(states)[]
    else
        return get_values!(zeros(eltype(states), length(fn.output_value_references)), wrapper, wrapper.outputs)
    end
end

@register_array_symbolic (fn::FMI2CSFunctor)(
    wrapper::FMI2InstanceWrapper, states::Vector{<:Real},
    inputs::Vector{<:Real}, params::Vector{<:Real}, t::Real
) begin
    size = (length(states) + length(fn.output_value_references),)
    eltype = eltype(states)
    ndims = 1
end

"""
    $(TYPEDSIGNATURES)

An affect function designed for use with `ImperativeAffect`. Should be triggered during
callback initialization. `m` should contain the key `:states` with the value being the
state vector if the FMU has continuous states. `m` should contain the key `:outputs` with
the value being the output vector if the FMU has output variables. `o` should contain the
`:inputs`, `:params`, `:t` and `:wrapper` where the latter contains the `FMI2InstanceWrapper`.

Initializes the FMU. Only for use with CoSimulation FMUs.
"""
function fmiCSInitialize!(m, o, ctx::FMI2CSFunctor, integrator)
    states = isdefined(m, :states) ? m.states : ()
    inputs = o.inputs
    params = o.params
    t = o.t
    wrapper = o.wrapper
    # every solve tells the FMU its own end time, unless one was fixed when the system
    # was built; a remade problem can have a different `tspan`
    ctx.t_end[] = isnan(ctx.stop_time) ? integrator.sol.prob.tspan[2] : ctx.stop_time
    if wrapper.instance !== nothing
        reset_instance!(wrapper)
    end

    instance = get_instance_CS!(wrapper, states, inputs, params, t, ctx.t_end[])
    ctx.t_last[] = t
    if isdefined(m, :states)
        @statuscheck FMI.fmi2GetReal!(instance, ctx.state_value_references, m.states)
    end
    if isdefined(m, :outputs)
        get_values!(m.outputs, wrapper, wrapper.outputs)
    end

    return m
end

"""
    $(TYPEDSIGNATURES)

An affect function designed for use with `ImperativeAffect`. Should be triggered
periodically to communicte with the CoSimulation FMU. Has the same requirements as
`fmiCSInitialize!` for `m` and `o`, with the addition that `o` should have a key
`:dt` with the value being the communication step size.
"""
function fmiCSStep!(m, o, ctx::FMI2CSFunctor, integrator)
    return cs_step!(m, o, ctx, integrator, integrator.t - o.dt, o.dt)
end

"""
    $(TYPEDSIGNATURES)

Step the CoSimulation FMU from `from` by `h` and read its states and outputs into `m`.
"""
function cs_step!(m, o, ctx::FMI2CSFunctor, integrator, from, h)
    wrapper = o.wrapper
    states = isdefined(m, :states) ? m.states : ()
    inputs = o.inputs
    params = o.params

    instance = get_instance_CS!(wrapper, states, inputs, params, integrator.t, ctx.t_end[])
    set_values!(wrapper, wrapper.inputs, inputs)
    @statuscheck FMI.fmi2DoStep(instance, from, h, FMI.fmi2True)
    ctx.t_last[] = from + h

    if isdefined(m, :states)
        @statuscheck FMI.fmi2GetReal!(instance, ctx.state_value_references, m.states)
    end
    if isdefined(m, :outputs)
        get_values!(m.outputs, wrapper, wrapper.outputs)
    end

    return m
end

"""
    $(TYPEDEF)

A callable struct useful for initializing v3 CoSimulation FMUs. When called, updates the
internal state of the FMU and gets updated values for output variables.

# Fields

$(TYPEDFIELDS)
"""
struct FMI3CSFunctor
    """
    The value references of state variables in the FMU.
    """
    state_value_references::Vector{FMI.fmi3ValueReference}
    """
    The value references of output variables in the FMU.
    """
    output_value_references::Vector{FMI.fmi3ValueReference}
    """
    `tspan[2]`, set during callback initialization
    """
    t_end::Base.RefValue{FMI.fmi3Float64}
    """
    The time the FMU has been stepped to, so the step to the end of the solve can be
    taken when the periodic callback does not fire there.
    """
    t_last::Base.RefValue{FMI.fmi3Float64}
    """
    The `stop_time` given to `FMIComponent`, or `NaN` to use the end of each solve's `tspan`.
    """
    stop_time::FMI.fmi3Float64
end

FMI3CSFunctor(states, outputs, t_end) = FMI3CSFunctor(
    states, outputs, t_end, Base.Ref{FMI.fmi3Float64}(NaN), t_end[]
)

function (fn::FMI3CSFunctor)(wrapper::FMI3InstanceWrapper, states, inputs, params, t)
    states = states isa SubArray ? copy(states) : states
    inputs = inputs isa SubArray ? copy(inputs) : inputs
    params = params isa SubArray ? copy(params) : params
    instance = get_instance_CS!(wrapper, states, inputs, params, t, fn.t_end[])

    if isempty(fn.output_value_references)
        return eltype(states)[]
    else
        return get_values!(zeros(eltype(states), length(fn.output_value_references)), wrapper, wrapper.outputs)
    end
end

@register_array_symbolic (fn::FMI3CSFunctor)(
    wrapper::FMI3InstanceWrapper, states::Vector{<:Real},
    inputs::Vector{<:Real}, params::Vector{<:Real}, t::Real
) begin
    size = (length(states) + length(fn.output_value_references),)
    eltype = eltype(states)
    ndims = 1
end

"""
    $(TYPEDSIGNATURES)
"""
function fmiCSInitialize!(m, o, ctx::FMI3CSFunctor, integrator)
    states = isdefined(m, :states) ? m.states : ()
    inputs = o.inputs
    params = o.params
    t = o.t
    wrapper = o.wrapper
    # every solve tells the FMU its own end time, unless one was fixed when the system
    # was built; a remade problem can have a different `tspan`
    ctx.t_end[] = isnan(ctx.stop_time) ? integrator.sol.prob.tspan[2] : ctx.stop_time
    if wrapper.instance !== nothing
        reset_instance!(wrapper)
    end
    instance = get_instance_CS!(wrapper, states, inputs, params, t, ctx.t_end[])
    ctx.t_last[] = t
    if isdefined(m, :states)
        @statuscheck FMI.fmi3GetFloat64!(instance, ctx.state_value_references, m.states)
    end
    if isdefined(m, :outputs)
        get_values!(m.outputs, wrapper, wrapper.outputs)
    end

    return m
end

"""
    $(TYPEDSIGNATURES)
"""
function fmiCSStep!(m, o, ctx::FMI3CSFunctor, integrator)
    return cs_step!(m, o, ctx, integrator, integrator.t - o.dt, o.dt)
end

"""
    $(TYPEDSIGNATURES)

Step the CoSimulation FMU from `from` by `h` and read its states and outputs into `m`.
"""
function cs_step!(m, o, ctx::FMI3CSFunctor, integrator, from, h)
    wrapper = o.wrapper
    states = isdefined(m, :states) ? m.states : ()
    inputs = o.inputs
    params = o.params

    instance = get_instance_CS!(wrapper, states, inputs, params, integrator.t, ctx.t_end[])
    set_values!(wrapper, wrapper.inputs, inputs)
    eventEncountered = Ref(FMI.fmi3False)
    terminateSimulation = Ref(FMI.fmi3False)
    earlyReturn = Ref(FMI.fmi3False)
    lastSuccessfulTime = Ref(zero(FMI.fmi3Float64))
    @statuscheck FMI.fmi3DoStep!(
        instance, from, h, FMI.fmi3True, eventEncountered,
        terminateSimulation, earlyReturn, lastSuccessfulTime
    )
    ctx.t_last[] = from + h
    @assert eventEncountered[] == FMI.fmi3False
    @assert terminateSimulation[] == FMI.fmi3False
    @assert earlyReturn[] == FMI.fmi3False

    if isdefined(m, :states)
        @statuscheck FMI.fmi3GetFloat64!(instance, ctx.state_value_references, m.states)
    end
    if isdefined(m, :outputs)
        get_values!(m.outputs, wrapper, wrapper.outputs)
    end

    return m
end

end # module
