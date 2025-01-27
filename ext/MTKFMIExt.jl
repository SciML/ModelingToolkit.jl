module MTKFMIExt

using ModelingToolkit
using SymbolicIndexingInterface
using ModelingToolkit: t_nounits as t, D_nounits as D
using DocStringExtensions
import ModelingToolkit as MTK
import SciMLBase
import FMI

"""
    $(TYPEDSIGNATURES)

A utility macro for FMI.jl functions that return a status. Will terminate on
fatal statuses. Must be used as `@statuscheck FMI.fmiXFunction(...)` where
`X` should be `2` or `3`. Has an edge case for handling tuples for
`FMI.fmi2CompletedIntegratorStep`.
"""
macro statuscheck(expr)
    @assert Meta.isexpr(expr, :call)
    fn = expr.args[1]
    @assert Meta.isexpr(fn, :.)
    @assert fn.args[1] == :FMI
    fnname = fn.args[2]

    instance = expr.args[2]
    is_v2 = startswith("fmi2", string(fnname))

    fmiTrue = is_v2 ? FMI.fmi2True : FMI.fmi3True
    fmiStatusOK = is_v2 ? FMI.fmi2StatusOK : FMI.fmi3StatusOK
    fmiStatusWarning = is_v2 ? FMI.fmi2StatusWarning : FMI.fmi3StatusWarning
    fmiStatusFatal = is_v2 ? FMI.fmi2StatusFatal : FMI.fmi3StatusFatal
    fmiTerminate = is_v2 ? FMI.fmi2Terminate : FMI.fmi3Terminate
    fmiFreeInstance! = is_v2 ? FMI.fmi2FreeInstance! : FMI.fmi3FreeInstance!
    return quote
        status = $expr
        fnname = $fnname
        if status !== nothing && ((status isa Tuple && status[1] == $fmiTrue) ||
            (!(status isa Tuple) && status != $fmiStatusOK &&
             status != $fmiStatusWarning))
            if status != $fmiStatusFatal
                $fmiTerminate(wrapper.instance)
            end
            $fmiFreeInstance!(wrapper.instance)
            wrapper.instance = nothing
            error("FMU Error in $fnname: status $status")
        end
    end |> esc
end

@static if !hasmethod(FMI.getValueReferencesAndNames, Tuple{FMI.fmi3ModelDescription})
    """
        $(TYPEDSIGNATURES)

        This is type piracy, but FMI.jl is missing this implementation. It allows
        `FMI.getStateValueReferencesAndNames` to work.
    """
    function FMI.getValueReferencesAndNames(
            md::FMI.fmi3ModelDescription; vrs = md.valueReferences)
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
All inputs, continuous variables and outputs must be `FMI.fmi2Real` or `FMI.fmi3Float64`.
Does not support events or discrete variables in the FMU. Does not support automatic
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
- `type`: Either `:ME` or `:CS` depending on whether `fmu` is a Model Exchange or
  CoSimulation FMU respectively.
- `name`: The name of the system.
"""
function MTK.FMIComponent(::Val{Ver}; fmu = nothing, tolerance = 1e-6,
        communication_step_size = nothing, reinitializealg = SciMLBase.NoInit(), type, name) where {Ver}
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

    # parse states
    fmi_variables_to_mtk_variables!(fmu, FMI.getStateValueReferencesAndNames(fmu),
        value_references, diffvars, states, observed)
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
        @parameters __mtk_internal_u(t)[1:length(diffvars)]=missing [guess = diffvars]
        push!(observed, __mtk_internal_u ~ copy(diffvars))
    end

    # parse derivatives of states
    # the variables passed to `postprocess_variable` haven't been differentiated yet, so they
    # should match one variable in states. That's the one this is the derivative of, and we
    # keep track of this ordering
    derivative_order = []
    function derivative_order_postprocess(var)
        idx = findfirst(isequal(var), states)
        idx === nothing || push!(derivative_order, states[idx])
        return var
    end
    fmi_variables_to_mtk_variables!(
        fmu, FMI.getDerivateValueReferencesAndNames(fmu), value_references, dervars,
        states, der_observed; postprocess_variable = derivative_order_postprocess)
    @assert length(derivative_order) == length(dervars)

    # parse the inputs to the FMU
    inputs = []
    fmi_variables_to_mtk_variables!(fmu, FMI.getInputValueReferencesAndNames(fmu),
        value_references, inputs, states, observed; postprocess_variable = v -> MTK.setinput(
            v, true))
    # create a symbolic variable for the input buffer
    __mtk_internal_x = copy(inputs)
    if isempty(__mtk_internal_x)
        __mtk_internal_x = Float64[]
    end

    # parse the outputs of the FMU
    outputs = []
    fmi_variables_to_mtk_variables!(fmu, FMI.getOutputValueReferencesAndNames(fmu),
        value_references, outputs, states, observed; postprocess_variable = v -> MTK.setoutput(
            v, true))
    # create the output buffer. This is only required for CoSimulation to pass it to
    # the callback affect
    if type == :CS
        if isempty(outputs)
            __mtk_internal_o = Float64[]
        else
            @parameters __mtk_internal_o(t)[1:length(outputs)]=missing [guess = zeros(length(outputs))]
            push!(observed, __mtk_internal_o ~ outputs)
        end
    end

    # parse the parameters
    params = []
    # multiple names for the same parameter are treated as parameter dependencies.
    parameter_dependencies = Equation[]
    fmi_variables_to_mtk_variables!(
        fmu, FMI.getParameterValueReferencesAndNames(fmu), value_references,
        params, [], parameter_dependencies, defs; parameters = true)
    # create a symbolic variable for the parameter buffer
    __mtk_internal_p = copy(params)
    if isempty(__mtk_internal_p)
        __mtk_internal_p = Float64[]
    end

    derivative_value_references = UInt32[value_references[var] for var in dervars]
    state_value_references = UInt32[value_references[var] for var in diffvars]
    output_value_references = UInt32[value_references[var] for var in outputs]
    input_value_references = UInt32[value_references[var] for var in inputs]
    param_value_references = UInt32[value_references[var] for var in params]

    # create a parameter for the instance wrapper
    # this manages the creation and deallocation of FMU instances
    buffer_length = length(diffvars) + length(outputs)
    if Ver == 2
        @parameters (wrapper::FMI2InstanceWrapper)(..)[1:buffer_length] = FMI2InstanceWrapper(
            fmu, derivative_value_references, state_value_references, output_value_references,
            param_value_references, input_value_references, tolerance)
    else
        @parameters (wrapper::FMI3InstanceWrapper)(..)[1:buffer_length] = FMI3InstanceWrapper(
            fmu, derivative_value_references, state_value_references,
            output_value_references, param_value_references, input_value_references)
    end

    # any additional initialization equations for the system
    initialization_eqs = Equation[]

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

        # instance management callback which deallocates the instance when
        # necessary and notifies the FMU of completed integrator steps
        finalize_affect = MTK.FunctionalAffect(fmiFinalize!, [], [wrapper], [])
        step_affect = MTK.FunctionalAffect(Returns(nothing), [], [], [])
        instance_management_callback = MTK.SymbolicDiscreteCallback(
            (t != t - 1), step_affect; finalize = finalize_affect, reinitializealg = reinitializealg)

        push!(params, wrapper)
        append!(observed, der_observed)
    elseif type == :CS
        _functor = if Ver == 2
            FMI2CSFunctor(state_value_references, output_value_references)
        else
            FMI3CSFunctor(state_value_references, output_value_references)
        end
        @parameters (functor::(typeof(_functor)))(..)[1:(length(__mtk_internal_u) + length(__mtk_internal_o))] = _functor
        # for co-simulation, we need to ensure the output buffer is solved for
        # during initialization
        for (i, x) in enumerate(collect(__mtk_internal_o))
            push!(initialization_eqs,
                x ~ functor(
                    wrapper, __mtk_internal_u, __mtk_internal_x, __mtk_internal_p, t)[i])
        end

        diffeqs = Equation[]

        # use `ImperativeAffect` for instance management here
        cb_observed = (; inputs = __mtk_internal_x, params = copy(params),
            t, wrapper, dt = communication_step_size)
        cb_modified = (;)
        # modify the outputs if present
        if symbolic_type(__mtk_internal_o) != NotSymbolic()
            cb_modified = (cb_modified..., outputs = __mtk_internal_o)
        end
        # modify the continuous state if present
        if symbolic_type(__mtk_internal_u) != NotSymbolic()
            cb_modified = (cb_modified..., states = __mtk_internal_u)
        end
        initialize_affect = MTK.ImperativeAffect(fmiCSInitialize!; observed = cb_observed,
            modified = cb_modified, ctx = _functor)
        finalize_affect = MTK.FunctionalAffect(fmiFinalize!, [], [wrapper], [])
        # the callback affect performs the stepping
        step_affect = MTK.ImperativeAffect(
            fmiCSStep!; observed = cb_observed, modified = cb_modified, ctx = _functor)
        instance_management_callback = MTK.SymbolicDiscreteCallback(
            communication_step_size, step_affect; initialize = initialize_affect,
            finalize = finalize_affect, reinitializealg = reinitializealg
        )

        # guarded in case there are no outputs/states and the variable is `[]`.
        symbolic_type(__mtk_internal_o) == NotSymbolic() || push!(params, __mtk_internal_o)
        symbolic_type(__mtk_internal_u) == NotSymbolic() || push!(params, __mtk_internal_u)

        push!(params, wrapper, functor)
    end

    eqs = [observed; diffeqs]
    return ODESystem(eqs, t, states, params; parameter_dependencies, defaults = defs,
        discrete_events = [instance_management_callback], name, initialization_eqs)
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
- `postprocess_variable`: A function applied to each created variable that should
  return the updated variable. This is useful to add metadata to variables.
"""
function fmi_variables_to_mtk_variables!(
        fmu::Union{FMI.FMU2, FMI.FMU3}, varmap::AbstractDict,
        value_references::AbstractDict, truevars, allvars,
        obseqs, defs = Dict(); parameters = false, postprocess_variable = identity)
    for (valRef, varnames) in varmap
        stateT = FMI.dataTypeForValueReference(fmu, valRef)
        snames = Symbol[]
        ders = Int[]
        for name in varnames
            sname, der = parseFMIVariableName(name)
            push!(snames, sname)
            push!(ders, der)
        end
        if parameters
            vars = [postprocess_variable(MTK.unwrap(only(@parameters $sname::stateT)))
                    for sname in snames]
        else
            vars = [postprocess_variable(MTK.unwrap(only(@variables $sname(t)::stateT)))
                    for sname in snames]
        end
        for i in eachindex(vars)
            der = ders[i]
            vars[i] = MTK.unwrap(vars[i])
            for j in 1:der
                vars[i] = D(vars[i])
            end
            vars[i] = MTK.default_toterm(vars[i])
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
        defs[vars[1]] = defval
    end
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
        idx = findfirst(',', name)
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
end

"""
    $(TYPEDSIGNATURES)

Create an `FMI2InstanceWrapper` with no instance.
"""
function FMI2InstanceWrapper(fmu, ders, states, outputs, params, inputs, tolerance)
    FMI2InstanceWrapper(fmu, ders, states, outputs, params, inputs, tolerance, nothing)
end

Base.nameof(::FMI2InstanceWrapper) = :FMI2InstanceWrapper

"""
    $(TYPEDSIGNATURES)

Common functionality for creating an instance of a v2 FMU. Does not check if
`wrapper.instance` is already present, and overwrites the existing value with
a new instance. `inputs` should be in the order of `wrapper.input_value_references`.
`params` should be in the order of `wrapper.param_value_references`. `t` is the current
time. Returns the created instance, which is also stored in `wrapper.instance`.
"""
function get_instance_common!(wrapper::FMI2InstanceWrapper, inputs, params, t)
    wrapper.instance = FMI.fmi2Instantiate!(wrapper.fmu)::FMI.FMU2Component
    if !isempty(inputs)
        @statuscheck FMI.fmi2SetReal(wrapper.instance, wrapper.input_value_references,
            Csize_t(length(wrapper.param_value_references)), inputs)
    end
    if !isempty(params)
        @statuscheck FMI.fmi2SetReal(wrapper.instance, wrapper.param_value_references,
            Csize_t(length(wrapper.param_value_references)), params)
    end
    @statuscheck FMI.fmi2SetupExperiment(
        wrapper.instance, FMI.fmi2True, wrapper.tolerance, t, FMI.fmi2False, t)
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
    if wrapper.instance === nothing
        get_instance_common!(wrapper, inputs, params, t)
        @statuscheck FMI.fmi2ExitInitializationMode(wrapper.instance)
        eventInfo = FMI.fmi2NewDiscreteStates(wrapper.instance)
        @assert eventInfo.newDiscreteStatesNeeded == FMI.fmi2False
        # TODO: Support FMU events
        @statuscheck FMI.fmi2EnterContinuousTimeMode(wrapper.instance)
    end

    return wrapper.instance
end

"""
    $(TYPEDSIGNATURES)

Create an instance of a CoSimulation FMU. Use the existing instance in `wrapper` if
present and create a new one otherwise. Return the instance.

See `get_instance_common!` for a description of the arguments.
"""
function get_instance_CS!(wrapper::FMI2InstanceWrapper, states, inputs, params, t)
    if wrapper.instance === nothing
        get_instance_common!(wrapper, inputs, params, t)
        if !isempty(states)
            @statuscheck FMI.fmi2SetReal(wrapper.instance, wrapper.state_value_references,
                Csize_t(length(wrapper.state_value_references)), states)
        end
        @statuscheck FMI.fmi2ExitInitializationMode(wrapper.instance)
    end
    return wrapper.instance
end

"""
    $(TYPEDSIGNATURES)

Call `fmiXCompletedIntegratorStep` with `noSetFMUStatePriorToCurrentPoint` as false.
"""
function partiallyCompleteIntegratorStep(wrapper::FMI2InstanceWrapper)
    @statuscheck FMI.fmi2CompletedIntegratorStep(wrapper.instance, FMI.fmi2False)
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
    wrapper.instance = nothing
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
end

"""
    $(TYPEDSIGNATURES)

Create an `FMI3InstanceWrapper` with no instance.
"""
function FMI3InstanceWrapper(fmu, ders, states, outputs, params, inputs)
    FMI3InstanceWrapper(fmu, ders, states, outputs, params, inputs, nothing)
end

Base.nameof(::FMI3InstanceWrapper) = :FMI3InstanceWrapper

"""
    $(TYPEDSIGNATURES)

Common functionality for creating an instance of a v3 FMU. Since v3 FMUs need to be
instantiated differently depending on the type, this assumes `wrapper.instance` is a
freshly instantiated FMU which needs to be initialized. `inputs` should be in the order
of `wrapper.input_value_references`. `params` should be in the order of
`wrapper.param_value_references`. `t` is the current time. Returns `wrapper.instance`.
"""
function get_instance_common!(wrapper::FMI3InstanceWrapper, inputs, params, t)
    if !isempty(params)
        @statuscheck FMI.fmi3SetFloat64(wrapper.instance, wrapper.param_value_references,
            params)
    end
    @statuscheck FMI.fmi3EnterInitializationMode(
        wrapper.instance, FMI.fmi3False, zero(FMI.fmi3Float64), t, FMI.fmi3False, t)
    if !isempty(inputs)
        @statuscheck FMI.fmi3SetFloat64(
            wrapper.instance, wrapper.input_value_references, inputs)
    end

    return wrapper.instance
end

"""
    $(TYPEDSIGNATURES)

Create an instance of a Model Exchange FMU. Use the existing instance in `wrapper` if
present and create a new one otherwise. Return the instance.

See `get_instance_common!` for a description of the arguments.
"""
function get_instance_ME!(wrapper::FMI3InstanceWrapper, inputs, params, t)
    if wrapper.instance === nothing
        wrapper.instance = FMI.fmi3InstantiateModelExchange!(wrapper.fmu)::FMI.FMU3Instance
        get_instance_common!(wrapper, inputs, params, t)
        @statuscheck FMI.fmi3ExitInitializationMode(wrapper.instance)
        eventInfo = FMI.fmi3UpdateDiscreteStates(wrapper.instance)
        @assert eventInfo[1] == FMI.fmi2False
        # TODO: Support FMU events
        @statuscheck FMI.fmi3EnterContinuousTimeMode(wrapper.instance)
    end

    return wrapper.instance
end

"""
    $(TYPEDSIGNATURES)

Create an instance of a CoSimulation FMU. Use the existing instance in `wrapper` if
present and create a new one otherwise. Return the instance.

See `get_instance_common!` for a description of the arguments.
"""
function get_instance_CS!(wrapper::FMI3InstanceWrapper, states, inputs, params, t)
    if wrapper.instance === nothing
        wrapper.instance = FMI.fmi3InstantiateCoSimulation!(
            wrapper.fmu; eventModeUsed = false)::FMI.FMU3Instance
        get_instance_common!(wrapper, inputs, params, t)
        if !isempty(states)
            @statuscheck FMI.fmi3SetFloat64(
                wrapper.instance, wrapper.state_value_references, states)
        end
        @statuscheck FMI.fmi3ExitInitializationMode(wrapper.instance)
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
        wrapper.instance, FMI.fmi3False, enterEventMode, terminateSimulation)
    @assert enterEventMode[] == FMI.fmi3False
    @assert terminateSimulation[] == FMI.fmi3False
end

"""
    $(TYPEDSIGNATURES)
"""
function reset_instance!(wrapper::FMI3InstanceWrapper)
    wrapper.instance === nothing && return
    FMI.fmi3Terminate(wrapper.instance)
    FMI.fmi3FreeInstance!(wrapper.instance)
    wrapper.instance = nothing
end

@register_array_symbolic (fn::FMI2InstanceWrapper)(
    states::Vector{<:Real}, inputs::Vector{<:Real}, params::Vector{<:Real}, t::Real) begin
    size = (length(states) + length(fn.output_value_references),)
    eltype = eltype(states)
    ndims = 1
end

@register_array_symbolic (fn::FMI3InstanceWrapper)(
    wrapper::FMI3InstanceWrapper, states::Vector{<:Real},
    inputs::Vector{<:Real}, params::Vector{<:Real}, t::Real) begin
    size = (length(states) + length(fn.output_value_references),)
    eltype = eltype(states)
    ndims = 1
end

"""
    $(TYPEDSIGNATURES)

Update the internal state of the ME FMU and return a vector of updated values
for continuous state derivatives and output variables respectively. Needs to be a
callable struct to enable symbolic registration with an inferred return size.
"""
function (wrapper::Union{FMI2InstanceWrapper, FMI3InstanceWrapper})(
        states, inputs, params, t)
    instance = get_instance_ME!(wrapper, inputs, params, t)

    # TODO: Find a way to do this without allocating. We can't pass a view to these
    # functions.
    states_buffer = zeros(length(states))
    outputs_buffer = zeros(length(wrapper.output_value_references))
    # Defined in FMIBase.jl/src/eval.jl
    # Doesn't seem to be documented, but somehow this is the only way to
    # propagate inputs to the FMU consistently. I have no idea why.
    instance(; x = states, u = inputs, u_refs = wrapper.input_value_references,
        p = params, p_refs = wrapper.param_value_references, t = t)
    # the spec requires completing the step before getting updated derivative/output values
    partiallyCompleteIntegratorStep(wrapper)
    instance(; dx = states_buffer, dx_refs = wrapper.derivative_value_references,
        y = outputs_buffer, y_refs = wrapper.output_value_references)
    return [states_buffer; outputs_buffer]
end

"""
    $(TYPEDSIGNATURES)

An affect function for use inside a `FunctionalAffect`. This should be triggered at the
end of the solve, regardless of whether it succeeded or failed. Expects `p` to be a
1-length array containing the index of the instance wrapper (`FMI2InstanceWrapper` or
`FMI3InstanceWrapper`) in the parameter object.
"""
function fmiFinalize!(integrator, u, p, ctx)
    wrapper_idx = p[1]
    wrapper = integrator.ps[wrapper_idx]
    reset_instance!(wrapper)
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
end

function (fn::FMI2CSFunctor)(wrapper::FMI2InstanceWrapper, states, inputs, params, t)
    states = states isa SubArray ? copy(states) : states
    inputs = inputs isa SubArray ? copy(inputs) : inputs
    params = params isa SubArray ? copy(params) : params
    if wrapper.instance !== nothing
        reset_instance!(wrapper)
    end
    instance = get_instance_CS!(wrapper, states, inputs, params, t)
    if isempty(fn.output_value_references)
        return eltype(states)[]
    else
        return FMI.fmi2GetReal(instance, fn.output_value_references)
    end
end

@register_array_symbolic (fn::FMI2CSFunctor)(
    wrapper::FMI2InstanceWrapper, states::Vector{<:Real},
    inputs::Vector{<:Real}, params::Vector{<:Real}, t::Real) begin
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
    if wrapper.instance !== nothing
        reset_instance!(wrapper)
    end

    instance = get_instance_CS!(wrapper, states, inputs, params, t)
    if isdefined(m, :states)
        @statuscheck FMI.fmi2GetReal!(instance, ctx.state_value_references, m.states)
    end
    if isdefined(m, :outputs)
        @statuscheck FMI.fmi2GetReal!(instance, ctx.output_value_references, m.outputs)
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
    wrapper = o.wrapper
    states = isdefined(m, :states) ? m.states : ()
    inputs = o.inputs
    params = o.params
    t = o.t
    dt = o.dt

    instance = get_instance_CS!(wrapper, states, inputs, params, integrator.t)
    if !isempty(inputs)
        FMI.fmi2SetReal(
            instance, wrapper.input_value_references, Csize_t(length(inputs)), inputs)
    end
    @statuscheck FMI.fmi2DoStep(instance, integrator.t - dt, dt, FMI.fmi2True)

    if isdefined(m, :states)
        @statuscheck FMI.fmi2GetReal!(instance, ctx.state_value_references, m.states)
    end
    if isdefined(m, :outputs)
        @statuscheck FMI.fmi2GetReal!(instance, ctx.output_value_references, m.outputs)
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
end

function (fn::FMI3CSFunctor)(wrapper::FMI3InstanceWrapper, states, inputs, params, t)
    states = states isa SubArray ? copy(states) : states
    inputs = inputs isa SubArray ? copy(inputs) : inputs
    params = params isa SubArray ? copy(params) : params
    instance = get_instance_CS!(wrapper, states, inputs, params, t)

    if isempty(fn.output_value_references)
        return eltype(states)[]
    else
        return FMI.fmi3GetFloat64(instance, fn.output_value_references)
    end
end

@register_array_symbolic (fn::FMI3CSFunctor)(
    wrapper::FMI3InstanceWrapper, states::Vector{<:Real},
    inputs::Vector{<:Real}, params::Vector{<:Real}, t::Real) begin
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
    if wrapper.instance !== nothing
        reset_instance!(wrapper)
    end
    instance = get_instance_CS!(wrapper, states, inputs, params, t)
    if isdefined(m, :states)
        @statuscheck FMI.fmi3GetFloat64!(instance, ctx.state_value_references, m.states)
    end
    if isdefined(m, :outputs)
        @statuscheck FMI.fmi3GetFloat64!(instance, ctx.output_value_references, m.outputs)
    end

    return m
end

"""
    $(TYPEDSIGNATURES)
"""
function fmiCSStep!(m, o, ctx::FMI3CSFunctor, integrator)
    wrapper = o.wrapper
    states = isdefined(m, :states) ? m.states : ()
    inputs = o.inputs
    params = o.params
    t = o.t
    dt = o.dt

    instance = get_instance_CS!(wrapper, states, inputs, params, integrator.t)
    if !isempty(inputs)
        FMI.fmi3SetFloat64(instance, wrapper.input_value_references, inputs)
    end
    eventEncountered = Ref(FMI.fmi3False)
    terminateSimulation = Ref(FMI.fmi3False)
    earlyReturn = Ref(FMI.fmi3False)
    lastSuccessfulTime = Ref(zero(FMI.fmi3Float64))
    @statuscheck FMI.fmi3DoStep!(
        instance, integrator.t - dt, dt, FMI.fmi3True, eventEncountered,
        terminateSimulation, earlyReturn, lastSuccessfulTime)
    @assert eventEncountered[] == FMI.fmi3False
    @assert terminateSimulation[] == FMI.fmi3False
    @assert earlyReturn[] == FMI.fmi3False

    if isdefined(m, :states)
        @statuscheck FMI.fmi3GetFloat64!(instance, ctx.state_value_references, m.states)
    end
    if isdefined(m, :outputs)
        @statuscheck FMI.fmi3GetFloat64!(instance, ctx.output_value_references, m.outputs)
    end

    return m
end

end # module
