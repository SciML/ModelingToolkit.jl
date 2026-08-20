module MTKFMIExt

using ModelingToolkit
using SymbolicIndexingInterface: NotSymbolic, symbolic_type, getname
using ModelingToolkitBase: t_nounits, D_nounits
using Symbolics: NAMESPACE_SEPARATOR, CallAndWrap
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ModelingToolkit as MTK
import ModelingToolkitBase as MTKBase
# for `BrownFullBasicInit`, which SciMLBase 3.49 does not define
import DiffEqBase
import SciMLBase
import SymbolicIndexingInterface as SII
import SymbolicUtils
import FMIImport as FMI

# Spelled as `const` rather than `t_nounits as t`, because ExplicitImports reports a
# renaming import as stale (it tracks the pre-rename name, which never appears in use).
const t = t_nounits
const D = D_nounits

"""
    $(TYPEDSIGNATURES)

Whether `status`, returned by the FMI version 2 function `fnname`, is a failure.
`fmi2StatusWarning` is not a failure: the standard allows the computation to continue, so it
is only logged.
"""
function fmu_status_is_error(::Val{2}, status, fnname)
    status === nothing && return false
    status == FMI.fmi2StatusOK && return false
    if status == FMI.fmi2StatusWarning
        @warn "FMU function $fnname returned `fmi2StatusWarning`." maxlog = 10
        return false
    end
    return true
end

"""
    $(TYPEDSIGNATURES)

Whether `status`, returned by the FMI version 3 function `fnname`, is a failure.
`fmi3StatusWarning` is not a failure: the standard allows the computation to continue, so it
is only logged.
"""
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
    $(TYPEDSIGNATURES)

A utility macro for FMI.jl functions that return a status. Frees the instance and errors on a
failing status, terminating it first unless the status is fatal. Must be used as
`@statuscheck FMI.fmiXFunction(...)` where `X` should be `2` or `3`. Evaluates to the value
returned by the wrapped call, so callers can consume the extra outputs of functions such as
`FMI.fmi2CompletedIntegratorStep`, which return their status as the first element of a tuple.
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
    fmiStatusFatal = is_v2 ? FMI.fmi2StatusFatal : FMI.fmi3StatusFatal
    fmiTerminate = is_v2 ? FMI.fmi2Terminate : FMI.fmi3Terminate
    fmiFreeInstance! = is_v2 ? FMI.fmi2FreeInstance! : FMI.fmi3FreeInstance!
    return quote
        result = $expr
        status = result isa Tuple ? result[1] : result
        if fmu_status_is_error($version, status, $(QuoteNode(fnname)))
            # a fatal status leaves the instance in a state where no further call is allowed
            if status != $fmiStatusFatal
                $fmiTerminate(wrapper.instance)
            end
            $fmiFreeInstance!(wrapper.instance)
            wrapper.instance = nothing
            error("FMU Error in $($(QuoteNode(fnname))): status $status")
        end
        result
    end |> esc
end

"""
The maximum number of `fmiXNewDiscreteStates`/`fmiXUpdateDiscreteStates` calls allowed in a
single event iteration before it is considered non-convergent.
"""
const MAX_EVENT_ITERATIONS = 100

"""
    $(TYPEDSIGNATURES)

Record the result of the event iteration performed while leaving initialization mode on
`wrapper`, and error if the FMU requested termination or declared a time event.
"""
function handle_initial_event_iteration!(wrapper, event_result)
    if event_result.terminate
        # the instance is in Event Mode, where the state exchange a retried solve performs is
        # illegal, so a caught error must not leave it reusable. Terminating from Event Mode is
        # legal in both FMI versions.
        reset_instance!(wrapper)
        error(
            "FMU $(FMI.getModelName(wrapper.fmu)) requested termination of the simulation \
            during the event iteration performed while leaving initialization mode."
        )
    end
    wrapper.next_event_time = event_result.next_event_time
    # this is what makes `ignore_time_events` warn once per FMU instance. A solve creates one
    # for a Model Exchange FMU, and at least two for a CoSimulation one: the functor
    # instantiates during problem initialization and `fmiCSInitialize!` re-instantiates.
    wrapper.time_event_warned = false
    handle_fmu_time_event!(wrapper, event_result.next_event_time)
    return event_result
end

"""
    $(TYPEDSIGNATURES)

Report a time event declared by the FMU at `next_event_time`, or do nothing if it declared
none. Time events are not supported, so this errors unless `ignore_time_events` was passed to
`FMIComponent`, in which case it warns once per FMU instance.
"""
function handle_fmu_time_event!(wrapper, next_event_time)
    next_event_time === nothing && return nothing
    name = FMI.getModelName(wrapper.fmu)
    if !wrapper.ignore_time_events
        # a caught error must not leave the instance behind: the instance getters would hand
        # it, mid-solve state and all, to a retried solve
        reset_instance!(wrapper)
        error(
            "FMU $name declared a time event at t = $next_event_time. Time events are not \
            supported by the FMU import. Pass `ignore_time_events = true` to \
            `FMIComponent` to ignore them and continue with a warning instead."
        )
    end
    wrapper.time_event_warned && return nothing
    wrapper.time_event_warned = true
    @warn "FMU $name declared a time event at t = $next_event_time, which the FMU import \
        does not support. Ignoring it, as requested by `ignore_time_events = true`."
    return nothing
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
    $(TYPEDEF)

Event-relevant facts parsed from an FMU's model description at import time. Attached to the
FMU's instance wrapper parameter as Symbolics metadata, with this struct as the key, so that
consumers can find it by scanning `parameters` of the compiled system.

Names are component-relative because metadata is opaque to `renamespace`; use
[`resolve_relative`](@ref) to obtain the name in the namespace of the wrapper parameter.

The fields this version's runtime does not consume (`fmi_version`, `interface`, `state_vrs`,
`can_have_time_events` and `cs_has_event_mode` outside import), along with the
`nominals_changed` result of the event iteration, are inputs for the clocks follow-up.

# Fields

$(TYPEDFIELDS)
"""
Base.@kwdef struct FMUEventMetadata
    """
    The FMI version of the FMU, `2` or `3`.
    """
    fmi_version::Int
    """
    The FMU interface in use, `:ME` or `:CS`.
    """
    interface::Symbol
    """
    The number of event indicators the FMU declares.
    """
    n_event_indicators::Int
    """
    Whether the FMU may signal time events. `true` for Model Exchange FMUs.
    """
    can_have_time_events::Bool
    """
    Whether event mode is supported by the interface in use. Always `false` for Model
    Exchange imports and for FMI2 CoSimulation.
    """
    cs_has_event_mode::Bool
    """
    Component-relative names of the continuous states of the FMU.
    """
    state_names::Vector{Symbol}
    """
    Component-relative names of the inputs of the FMU.
    """
    input_names::Vector{Symbol}
    """
    Value references of the continuous states of the FMU, ordered as `state_names`.
    """
    state_vrs::Vector
end

"""
    $(TYPEDSIGNATURES)

Return the name of the variable `name` of the FMU component that `wrapper_param` belongs to,
in the namespace of `wrapper_param`. `name` is a component-relative name from
[`FMUEventMetadata`](@ref); an un-namespaced `wrapper_param` returns it unchanged.
"""
function resolve_relative(wrapper_param, name::Symbol)
    wrapper_name = string(getname(wrapper_param))
    separator = findlast(NAMESPACE_SEPARATOR, wrapper_name)
    separator === nothing && return name
    return Symbol(wrapper_name[1:separator], name)
end

"""
    $(TYPEDSIGNATURES)

Whether the FMU declares support for CoSimulation event mode. FMI2 CoSimulation has no event
mode, so this is always `false` for v2 FMUs.
"""
cs_event_mode_supported(::FMI.FMU2) = false

function cs_event_mode_supported(fmu::FMI.FMU3)
    coSimulation = fmu.modelDescription.coSimulation
    return coSimulation !== nothing && coSimulation.hasEventMode === true
end

"""
    $(TYPEDSIGNATURES)

Build the callbacks handling the events of the Model Exchange FMU whose instance wrapper is
the parameter `wrapper_param` of the compiled system `sys`. Attached to `wrapper_param` as
`ModelingToolkit.CallbackConstructionHook` metadata, which calls this during problem
construction.

The state events of the FMU become a single `VectorContinuousCallback` over all of its event
indicators, which an FMU declaring none does not get. Every Model Exchange FMU also gets a
`DiscreteCallback` with an always-true condition, which notifies it of each accepted
integrator step and handles the step events it may report there.
"""
function build_fmu_me_callbacks(sys, wrapper_param)
    meta = SymbolicUtils.getmetadata(
        SymbolicUtils.unwrap(wrapper_param), FMUEventMetadata
    )
    get_wrapper = SII.getp(sys, wrapper_param)
    # `nothing` while the FMU's state and inputs can be exchanged with the integrator, and the
    # reason they cannot otherwise. Handling an event needs both: the continuous states to
    # write the post-event values back into, and the inputs to evaluate at the event point.
    no_event_access = nothing
    state_idxs = Int[]
    for name in meta.state_names
        resolved = resolve_relative(wrapper_param, name)
        idx = SII.variable_index(sys, resolved)
        if idx === nothing
            no_event_access = "The continuous state $resolved of the FMU wrapped by \
                $(getname(wrapper_param)) is not an unknown of the simplified system, so \
                the value it takes after an FMU event cannot be written back. Prevent it \
                from being simplified away, for example with `irreducible = true`."
            break
        end
        push!(state_idxs, idx)
    end
    input_names = [resolve_relative(wrapper_param, name) for name in meta.input_names]
    get_inputs = nothing
    if no_event_access === nothing && !isempty(input_names)
        try
            # the root search evaluates the condition at interpolated states, so the inputs
            # have to be evaluated there too rather than read off the integrator
            get_inputs = SII.observed(sys, input_names)
        catch err
            # building the getter is the only honest test of whether the inputs can be
            # evaluated: `is_observed` is true for every symbol the system does not know as
            # something else, including one it does not contain at all
            no_event_access = "The inputs $input_names of the FMU wrapped by \
                $(getname(wrapper_param)) cannot be evaluated as functions of the state of \
                the simplified system, so the FMU cannot be moved to an event point. Prevent \
                them from being simplified away, for example with `irreducible = true`. \
                Building their getter failed with: $(sprint(showerror, err))"
        end
    end
    # a state event is only meaningful if its result can be written back, so an FMU that
    # declares event indicators has to have that access; one that declares none is still told
    # about completed steps, and only a step event it requests then needs the write-back
    if no_event_access !== nothing && meta.n_event_indicators > 0
        error(no_event_access)
    end
    state_buffer = zeros(Float64, length(state_idxs))
    input_buffer = zeros(Float64, length(input_names))

    function push_fmu_event_state!(wrapper, u, p, t)
        for (i, idx) in enumerate(state_idxs)
            state_buffer[i] = u[idx]
        end
        if get_inputs !== nothing
            copyto!(input_buffer, get_inputs(u, p, t))
        end
        return force_set_fmu_state!(wrapper, state_buffer, input_buffer, t)
    end

    function fmu_event_condition!(out, u, t, integrator)
        wrapper = get_wrapper(integrator)
        # the instance is created by the first evaluation of the FMU, and unlike a completed
        # step an event indicator has no value without one
        if wrapper.instance === nothing
            error(
                "The event indicators of the FMU wrapped by $(getname(wrapper_param)) cannot \
                be evaluated before the FMU is instantiated, which the first evaluation of \
                its dynamics does."
            )
        end
        push_fmu_event_state!(wrapper, u, integrator.p, t)
        return get_fmu_event_indicators!(wrapper, out)
    end

    function fmu_event_affect!(integrator, idx)
        wrapper = get_wrapper(integrator)
        push_fmu_event_state!(wrapper, integrator.u, integrator.p, integrator.t)
        enter_fmu_event_mode!(wrapper, idx)
        event_result = do_fmu_event_iteration!(wrapper)
        if event_result.terminate
            SciMLBase.terminate!(integrator)
            # a terminating FMU is left in Event Mode, where `SetContinuousStates` and
            # `CompletedIntegratorStep` are illegal, so the interpolant rebuild that follows
            # an affect must not evaluate the RHS. Nothing is stepped after `terminate!`.
            SciMLBase.derivative_discontinuity!(integrator, false)
            return nothing
        end
        states = leave_fmu_event_mode!(wrapper, event_result.values_changed)
        if states !== nothing
            for (i, state_idx) in enumerate(state_idxs)
                integrator.u[state_idx] = states[i]
            end
        end
        wrapper.next_event_time = event_result.next_event_time
        handle_fmu_time_event!(wrapper, event_result.next_event_time)
        return nothing
    end

    # notify the FMU of the step and return whether it asks for Event Mode, which is the only
    # outcome its caller has to treat as a discontinuity
    function fmu_step_event_occurred!(integrator)
        # a terminating event leaves the FMU in Event Mode, where none of the calls below are
        # legal, and SciML applies the discrete callbacks after a continuous one terminated.
        # `check_error!` stores `Success` in the retcode of a solve that is still in flight,
        # so those are the two codes that mean "nothing has ended this solve".
        retcode = integrator.sol.retcode
        in_flight = retcode === SciMLBase.ReturnCode.Default ||
            retcode === SciMLBase.ReturnCode.Success
        in_flight || return false
        wrapper = get_wrapper(integrator)
        # the instance is created by the first evaluation of the FMU, which an FMU that only
        # feeds observed equations may not have seen by the time this initializes
        wrapper.instance === nothing && return false
        # error control and (for an implicit solver) the Jacobian leave the FMU at points
        # other than the accepted one, which is the step it has to be told about. Without
        # exchangeable states there is nothing to push, which is what FMI.jl's `stepCompleted`
        # does for every FMU.
        no_event_access === nothing &&
            push_fmu_event_state!(wrapper, integrator.u, integrator.p, integrator.t)
        step_result = partiallyCompleteIntegratorStep(wrapper)
        if step_result.terminate
            SciMLBase.terminate!(integrator)
            return false
        end
        return step_result.enter_event_mode
    end

    function fmu_step_completed!(integrator)
        if fmu_step_event_occurred!(integrator)
            if no_event_access !== nothing
                error(
                    "The FMU declares no event indicators but requested Event Mode from a \
                    completed integrator step, which needs the same access to its states as \
                    a state event. $no_event_access"
                )
            end
            # every callback ahead of this one clears the flag during the callback-initialize
            # phase (`SciMLBase.INITIALIZE_DEFAULT`), so the write-back has to declare itself.
            # On the `apply_discrete_callback!` path the flag is already true.
            SciMLBase.derivative_discontinuity!(integrator, true)
            # a step event has no triggered event indicator
            return fmu_event_affect!(integrator, nothing)
        end
        # every affect is assumed to be a discontinuity, which for a system with a singular
        # mass matrix reinitializes and thereby overwrites `u`. A step the FMU accepts as-is
        # changes nothing, so this must clear the assumption on every non-event path.
        SciMLBase.derivative_discontinuity!(integrator, false)
        return nothing
    end

    function fmu_step_initialize(cb, u, t, integrator)
        return fmu_step_completed!(integrator)
    end

    # a callback that leaves a discontinuity standing has the integrator reinitialize, which
    # for a singular mass matrix means solving for consistent initial conditions again. The
    # problem's own algorithm is MTK's `OverrideInit`, which resets `u` to `u0` and would throw
    # away exactly the post-event states the affect just wrote, so both callbacks name an
    # algorithm that keeps the differential states and re-solves only the algebraic variables.
    event_initializealg = DiffEqBase.BrownFullBasicInit()

    callbacks = SciMLBase.DECallback[]
    if meta.n_event_indicators > 0
        push!(
            callbacks,
            SciMLBase.VectorContinuousCallback(
                fmu_event_condition!, fmu_event_affect!, meta.n_event_indicators;
                # FMI defines an event as a domain change, so the FMU has to be past the
                # switch when it enters Event Mode; localizing on the pre-crossing side would
                # have it refresh its relations to the old domain and not handle the event
                rootfind = SciMLBase.RightRootFind, interp_points = 10,
                save_positions = (true, true), initializealg = event_initializealg
            )
        )
    end
    push!(
        callbacks,
        SciMLBase.DiscreteCallback(
            Returns(true), fmu_step_completed!;
            # `initialize` gives the FMU the `t0` step FMI.jl's `func_start` does
            initialize = fmu_step_initialize, save_positions = (false, false),
            initializealg = event_initializealg
        )
    )
    return callbacks
end

"""
    $(TYPEDSIGNATURES)

A component that wraps an FMU loaded via FMI.jl. The FMI version (2 or 3) should be
provided as a `Val` to the function. Supports Model Exchange and CoSimulation FMUs.
All inputs, continuous variables and outputs must be `FMI.fmi2Real` or `FMI.fmi3Float64`.
Supports the state and step events of Model Exchange FMUs, and the event mode of FMI3
CoSimulation FMUs that declare `hasEventMode`; does not support time events or discrete
variables in the FMU. Does not support automatic differentiation.
Parameters of the FMU will have defaults corresponding to their initial
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
- `ignore_time_events`: Time events are not supported, and by default an FMU declaring one
  errors. If this is `true`, a declared time event is instead warned about once per FMU
  instance and ignored, which changes the trajectory the FMU would otherwise have taken. A
  solve creates one instance for a Model Exchange FMU and at least two for a CoSimulation one.
- `name`: The name of the system.
"""
function MTK.FMIComponent(
        ::Val{Ver}; fmu = nothing, tolerance = 1.0e-6,
        communication_step_size = nothing, reinitializealg = nothing,
        stop_time = nothing, ignore_time_events = false, type, name
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

    # parse states
    fmi_variables_to_mtk_variables!(
        fmu, FMI.getStateValueReferencesAndNames(fmu),
        value_references, diffvars, states, observed
    )
    # evaluating an instance sets continuous states and reads their derivatives positionally
    # (FMIBase ignores the value references given for them), so both vectors must follow the
    # FMU's canonical state order, which the value reference dictionaries do not preserve.
    canonical_state_index = Dict(
        vr => i for (i, vr) in enumerate(fmu.modelDescription.stateValueReferences)
    )
    canonical_position = var -> canonical_state_index[value_references[var]]
    sort!(diffvars; by = canonical_position)
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
    function derivative_order_postprocess(var)
        idx = findfirst(isequal(var), states)
        idx === nothing || push!(derivative_order, states[idx])
        return var
    end
    fmi_variables_to_mtk_variables!(
        fmu, FMI.getDerivateValueReferencesAndNames(fmu), value_references, dervars,
        states, der_observed; postprocess_variable = derivative_order_postprocess
    )
    @assert length(derivative_order) == length(dervars)
    # the returned derivative vector is ordered as the continuous states, so `dervars[i]` has
    # to be the derivative of `diffvars[i]`
    der_permutation = sortperm(derivative_order; by = canonical_position)
    permute!(dervars, der_permutation)
    permute!(derivative_order, der_permutation)

    # parse the inputs to the FMU
    inputs = []
    fmi_variables_to_mtk_variables!(
        fmu, FMI.getInputValueReferencesAndNames(fmu),
        value_references, inputs, states, observed; postprocess_variable = v -> MTKBase.setinput(
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
    mark_output = v -> MTKBase.setoutput(v, true)
    output_varmap = FMI.getOutputValueReferencesAndNames(fmu)
    # An FMU may declare a continuous state as an output. The state variable is that output,
    # so it is reused instead of being created again, which would duplicate the equation
    # defining it. Such an output also needs no runtime fetch and stays out of `outputs`
    # (and thus out of `output_value_references` and the wrapper buffers).
    state_by_value_reference = Dict(value_references[var] => var for var in diffvars)
    is_state_alias = kvp -> haskey(state_by_value_reference, kvp.first)
    state_output_varmap = filter(is_state_alias, output_varmap)
    fmi_variables_to_mtk_variables!(
        fmu, filter(!is_state_alias, output_varmap),
        value_references, outputs, states, observed; postprocess_variable = mark_output
    )
    # The names of a value reference do not depend on its category, so the states pass has
    # already created and equated every name of a state-aliased value reference. Names it did
    # not create are aliases of the state and need an equation of their own.
    aliased_outputs = []
    fmi_variables_to_mtk_variables!(
        fmu, state_output_varmap, value_references, [], aliased_outputs, Equation[];
        postprocess_variable = mark_output
    )
    for var in aliased_outputs
        any(isequal(var), states) && continue
        push!(states, var)
        push!(observed, var ~ state_by_value_reference[value_references[var]])
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
    input_value_references = UInt32[value_references[var] for var in inputs]
    param_value_references = UInt32[value_references[var] for var in params]

    # create a parameter for the instance wrapper
    # this manages the creation and deallocation of FMU instances
    buffer_length = length(diffvars) + length(outputs)
    # FMI only allows setting a variable declared as an output while the FMU is in
    # initialization mode, so a CoSimulation FMU whose continuous states are outputs owns them
    # once it is initialized.
    states_settable_after_init = isempty(state_output_varmap)
    if Ver == 2
        @parameters (wrapper::FMI2InstanceWrapper)(..)[1:buffer_length] = FMI2InstanceWrapper(
            fmu, derivative_value_references, state_value_references, output_value_references,
            param_value_references, input_value_references, tolerance,
            states_settable_after_init, ignore_time_events
        )
    else
        @parameters (wrapper::FMI3InstanceWrapper)(..)[1:buffer_length] = FMI3InstanceWrapper(
            fmu, derivative_value_references, state_value_references,
            output_value_references, param_value_references, input_value_references,
            states_settable_after_init, ignore_time_events
        )
    end

    event_metadata = FMUEventMetadata(;
        fmi_version = Ver, interface = type,
        n_event_indicators = Int(something(FMI.getNumberOfEventIndicators(fmu), 0)),
        can_have_time_events = type == :ME,
        cs_has_event_mode = type == :CS && cs_event_mode_supported(fmu),
        state_names = Symbol[getname(var) for var in diffvars],
        input_names = Symbol[getname(var) for var in inputs],
        state_vrs = state_value_references
    )
    # `setmetadata` returns a new symbolic and unwraps the callable parameter, so the
    # rewrapped result is what has to be spliced into the system below.
    tagged_wrapper = SymbolicUtils.setmetadata(
        SymbolicUtils.unwrap(wrapper), FMUEventMetadata, event_metadata
    )
    # every Model Exchange FMU gets the hook: the per-accepted-step notification is
    # unconditional, whether or not the FMU declares event indicators.
    if type == :ME
        tagged_wrapper = SymbolicUtils.setmetadata(
            tagged_wrapper, MTKBase.CallbackConstructionHook, build_fmu_me_callbacks
        )
    end
    wrapper = CallAndWrap(tagged_wrapper)

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

        # instance management callback, which deallocates the instance when necessary.
        # Completed integrator steps are notified by `build_fmu_me_callbacks`.
        finalize_affect = MTKBase.ImperativeAffect(fmiFinalize!; observed = (; wrapper))
        step_affect = MTKBase.ImperativeAffect(Returns((;)))
        instance_management_callback = MTKBase.SymbolicDiscreteCallback(
            (t == t - 1), step_affect; finalize = finalize_affect, reinitializealg = SciMLBase.NoInit()
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

        # guarded in case there are no outputs/states and the variable is `[]`.
        symbolic_type(__mtk_internal_o) == NotSymbolic() || push!(params, __mtk_internal_o)
        symbolic_type(__mtk_internal_u) == NotSymbolic() || push!(params, __mtk_internal_u)

        push!(params, wrapper, functor)
    end

    eqs = [observed; diffeqs]
    bindings = [eq.lhs => eq.rhs for eq in parameter_dependencies]
    return System(
        eqs, t, states, params; bindings, initial_conditions = defs,
        discrete_events = [instance_management_callback], name, initialization_eqs
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
- `postprocess_variable`: A function applied to each created variable that should
  return the updated variable. This is useful to add metadata to variables.
"""
function fmi_variables_to_mtk_variables!(
        fmu::Union{FMI.FMU2, FMI.FMU3}, varmap::AbstractDict,
        value_references::AbstractDict, truevars, allvars,
        obseqs, defs = Dict(); parameters = false, postprocess_variable = identity
    )
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
            vars = [
                postprocess_variable(SymbolicUtils.unwrap(only(@parameters $sname::stateT)))
                    for sname in snames
            ]
        else
            vars = [
                postprocess_variable(SymbolicUtils.unwrap(only(@variables $sname(t)::stateT)))
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
        defs[vars[1]] = defval
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
    Whether the continuous states of the FMU may be set after it leaves initialization mode.
    `false` if any of them is declared as an output, since FMI forbids setting an output
    outside initialization mode.
    """
    const states_settable_after_init::Bool
    """
    Whether a time event declared by the FMU should be warned about and ignored instead of
    being an error.
    """
    const ignore_time_events::Bool
    """
    The FMU instance, if present, and `nothing` otherwise.
    """
    instance::Union{FMI.FMU2Component{FMI.FMU2}, Nothing}
    """
    The time of the next event reported by the most recent event iteration, or `nothing` if
    the FMU did not declare one.
    """
    next_event_time::Union{Nothing, Float64}
    """
    Whether the warning about an ignored time event has already been emitted for the current
    instance.
    """
    time_event_warned::Bool
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
    Buffer for the FMU's continuous states, re-read from the FMU whenever an event iteration
    reports that the values of the continuous states changed.
    """
    const event_states_buffer::Vector{FMI.fmi2Real}
end

"""
    $(TYPEDSIGNATURES)

Create an `FMI2InstanceWrapper` with no instance.
"""
function FMI2InstanceWrapper(fmu, ders, states, outputs, params, inputs, tolerance, states_settable_after_init = true, ignore_time_events = false)
    return FMI2InstanceWrapper(fmu, ders, states, outputs, params, inputs, tolerance, states_settable_after_init, ignore_time_events, nothing, nothing, false, zeros(FMI.fmi2Real, length(states)), zeros(FMI.fmi2Real, length(outputs)), zeros(FMI.fmi2Real, length(states) + length(outputs)), zeros(FMI.fmi2Real, length(states)))
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
function get_instance_common!(wrapper::FMI2InstanceWrapper, inputs, params, t, t_end = nothing)
    wrapper.instance = FMI.fmi2Instantiate!(wrapper.fmu)::FMI.FMU2Component
    if !isempty(inputs)
        @statuscheck FMI.fmi2SetReal(
            wrapper.instance, wrapper.input_value_references,
            Csize_t(length(inputs)), inputs
        )
    end
    if !isempty(params)
        @statuscheck FMI.fmi2SetReal(
            wrapper.instance, wrapper.param_value_references,
            Csize_t(length(wrapper.param_value_references)), params
        )
    end
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

Run the FMU's event iteration to convergence. Assumes the FMU is already in Event Mode, and
neither enters nor leaves any mode itself. `values_changed` and `nominals_changed` are
accumulated across iterations, as the FMI event iteration requires every super-dense time
instant to be accounted for, not only the last one. `next_event_time` keeps the value of the
last iteration that declared one.
"""
function do_fmu_event_iteration!(wrapper::FMI2InstanceWrapper)
    values_changed = false
    nominals_changed = false
    next_event_time = nothing
    terminate = false
    for _ in 1:MAX_EVENT_ITERATIONS
        eventInfo = FMI.fmi2NewDiscreteStates(wrapper.instance)
        values_changed |= eventInfo.valuesOfContinuousStatesChanged == FMI.fmi2True
        nominals_changed |= eventInfo.nominalsOfContinuousStatesChanged == FMI.fmi2True
        if eventInfo.nextEventTimeDefined == FMI.fmi2True
            next_event_time = Float64(eventInfo.nextEventTime)
        end
        terminate = eventInfo.terminateSimulation == FMI.fmi2True
        if terminate || eventInfo.newDiscreteStatesNeeded == FMI.fmi2False
            return (; values_changed, nominals_changed, next_event_time, terminate)
        end
    end
    # the iteration is abandoned with the instance in Event Mode, which no retried solve can
    # exchange state from, so a caught error must not leave it reusable
    reset_instance!(wrapper)
    return error(
        "FMU $(FMI.getModelName(wrapper.fmu)) did not converge in $MAX_EVENT_ITERATIONS \
        event iterations."
    )
end

"""
    $(TYPEDSIGNATURES)

Push the continuous states, time and input values an event occurs at into the FMU. The v2
setters memoize on the values the instance was last set to, which the derivative evaluations
and the root search preceding an event leave at a different point, so they are forced.
"""
function force_set_fmu_state!(wrapper::FMI2InstanceWrapper, states, inputs, t)
    if !isempty(states)
        @statuscheck FMI.fmi2SetContinuousStates(wrapper.instance, states; force = true)
    end
    @statuscheck FMI.fmi2SetTime(wrapper.instance, t; force = true)
    if !isempty(inputs)
        @statuscheck FMI.fmi2SetReal(
            wrapper.instance, wrapper.input_value_references,
            Csize_t(length(inputs)), inputs
        )
    end
    return wrapper.instance
end

"""
    $(TYPEDSIGNATURES)

Read the FMU's event indicators into `out`, which must have one element per indicator.
"""
function get_fmu_event_indicators!(wrapper::FMI2InstanceWrapper, out)
    @statuscheck FMI.fmi2GetEventIndicators!(
        wrapper.instance, out, Csize_t(length(out))
    )
    return out
end

"""
    $(TYPEDSIGNATURES)

Read the FMU's outputs into `out`, which must have one element per output.
"""
function get_fmu_outputs!(wrapper::FMI2InstanceWrapper, out)
    isempty(out) && return out
    # the count has to be passed: the method without it returns `nothing` instead of a status
    @statuscheck FMI.fmi2GetReal!(
        wrapper.instance, wrapper.output_value_references, Csize_t(length(out)), out
    )
    return out
end

"""
    $(TYPEDSIGNATURES)

Enter Event Mode for the event described by `idx`: a step event if it is `nothing`, and the
state event on event indicator `idx` otherwise. The v2 function takes no arguments describing
the cause of the event.
"""
function enter_fmu_event_mode!(wrapper::FMI2InstanceWrapper, idx)
    @statuscheck FMI.fmi2EnterEventMode(wrapper.instance)
    return wrapper.instance
end

"""
    $(TYPEDSIGNATURES)

Leave Event Mode for Continuous-Time Mode. If `values_changed`, the event iteration changed
the continuous states of the FMU, which FMI then requires to be re-read; they are returned in
`wrapper.event_states_buffer`. Returns `nothing` otherwise.
"""
function leave_fmu_event_mode!(wrapper::FMI2InstanceWrapper, values_changed)
    @statuscheck FMI.fmi2EnterContinuousTimeMode(wrapper.instance)
    values_changed || return nothing
    buffer = wrapper.event_states_buffer
    @statuscheck FMI.fmi2GetContinuousStates!(
        wrapper.instance, buffer, Csize_t(length(buffer))
    )
    return buffer
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
        # leaving initialization mode puts the FMU in Event Mode, so the initial event
        # iteration has to converge before Continuous-Time Mode may be entered
        @statuscheck FMI.fmi2ExitInitializationMode(wrapper.instance)
        event_result = do_fmu_event_iteration!(wrapper)
        handle_initial_event_iteration!(wrapper, event_result)
        leave_fmu_event_mode!(wrapper, event_result.values_changed)
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
        if !isempty(states) && wrapper.states_settable_after_init
            @statuscheck FMI.fmi2SetReal(
                wrapper.instance, wrapper.state_value_references,
                Csize_t(length(wrapper.state_value_references)), states
            )
        end
        if !isempty(inputs)
            @statuscheck FMI.fmi2SetReal(
                wrapper.instance, wrapper.input_value_references,
                Csize_t(length(inputs)), inputs
            )
        end
    end
    return wrapper.instance
end

"""
    $(TYPEDSIGNATURES)

Call `fmiXCompletedIntegratorStep` with `noSetFMUStatePriorToCurrentPoint` as false. Returns
the FMU's requests to enter event mode and to terminate the simulation.
"""
function partiallyCompleteIntegratorStep(wrapper::FMI2InstanceWrapper)
    _, enter_event_mode, terminate = @statuscheck FMI.fmi2CompletedIntegratorStep(
        wrapper.instance, FMI.fmi2False
    )
    return (;
        enter_event_mode = enter_event_mode == FMI.fmi2True,
        terminate = terminate == FMI.fmi2True,
    )
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
    Whether the continuous states of the FMU may be set after it leaves initialization mode.
    `false` if any of them is declared as an output, since FMI forbids setting an output
    outside initialization mode.
    """
    const states_settable_after_init::Bool
    """
    Whether a time event declared by the FMU should be warned about and ignored instead of
    being an error.
    """
    const ignore_time_events::Bool
    """
    The FMU instance, if present, and `nothing` otherwise.
    """
    instance::Union{FMI.FMU3Instance{FMI.FMU3}, Nothing}
    """
    The time of the next event reported by the most recent event iteration, or `nothing` if
    the FMU did not declare one.
    """
    next_event_time::Union{Nothing, Float64}
    """
    Whether the warning about an ignored time event has already been emitted for the current
    instance.
    """
    time_event_warned::Bool
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
    Buffer for the FMU's continuous states, re-read from the FMU whenever an event iteration
    reports that the values of the continuous states changed.
    """
    const event_states_buffer::Vector{FMI.fmi3Float64}
end

"""
    $(TYPEDSIGNATURES)

Create an `FMI3InstanceWrapper` with no instance.
"""
function FMI3InstanceWrapper(fmu, ders, states, outputs, params, inputs, states_settable_after_init = true, ignore_time_events = false)
    return FMI3InstanceWrapper(fmu, ders, states, outputs, params, inputs, states_settable_after_init, ignore_time_events, nothing, nothing, false, zeros(FMI.fmi3Float64, length(states)), zeros(FMI.fmi3Float64, length(outputs)), zeros(FMI.fmi3Float64, length(states) + length(outputs)), zeros(FMI.fmi3Float64, length(states)))
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
function get_instance_common!(wrapper::FMI3InstanceWrapper, inputs, params, t, t_end = nothing)
    if !isempty(params)
        @statuscheck FMI.fmi3SetFloat64(
            wrapper.instance, wrapper.param_value_references,
            params
        )
    end
    stop_time_defined = t_end === nothing ? FMI.fmi3False : FMI.fmi3True
    stop_time = something(t_end, t)
    @statuscheck FMI.fmi3EnterInitializationMode(
        wrapper.instance, FMI.fmi3False, zero(FMI.fmi3Float64), t, stop_time_defined, stop_time
    )
    if !isempty(inputs)
        @statuscheck FMI.fmi3SetFloat64(
            wrapper.instance, wrapper.input_value_references, inputs
        )
    end

    return wrapper.instance
end

"""
    $(TYPEDSIGNATURES)

Run the FMU's event iteration to convergence. Assumes the FMU is already in Event Mode, and
neither enters nor leaves any mode itself. `values_changed` and `nominals_changed` are
accumulated across iterations, as the FMI event iteration requires every super-dense time
instant to be accounted for, not only the last one. `next_event_time` keeps the value of the
last iteration that declared one.
"""
function do_fmu_event_iteration!(wrapper::FMI3InstanceWrapper)
    values_changed = false
    nominals_changed = false
    next_event_time = nothing
    terminate = false
    for _ in 1:MAX_EVENT_ITERATIONS
        needs_update, terminate_simulation, states_nominals_changed,
            states_values_changed, event_time_defined,
            event_time = FMI.fmi3UpdateDiscreteStates(wrapper.instance)
        values_changed |= states_values_changed == FMI.fmi3True
        nominals_changed |= states_nominals_changed == FMI.fmi3True
        if event_time_defined == FMI.fmi3True
            next_event_time = Float64(event_time)
        end
        terminate = terminate_simulation == FMI.fmi3True
        if terminate || needs_update == FMI.fmi3False
            return (; values_changed, nominals_changed, next_event_time, terminate)
        end
    end
    # the iteration is abandoned with the instance in Event Mode, which no retried solve can
    # exchange state from, so a caught error must not leave it reusable
    reset_instance!(wrapper)
    return error(
        "FMU $(FMI.getModelName(wrapper.fmu)) did not converge in $MAX_EVENT_ITERATIONS \
        event iterations."
    )
end

"""
    $(TYPEDSIGNATURES)

Push the continuous states, time and input values an event occurs at into the FMU. The v3
setters do not memoize, so nothing has to be forced for them.
"""
function force_set_fmu_state!(wrapper::FMI3InstanceWrapper, states, inputs, t)
    if !isempty(states)
        @statuscheck FMI.fmi3SetContinuousStates(wrapper.instance, states)
    end
    @statuscheck FMI.fmi3SetTime(wrapper.instance, t)
    if !isempty(inputs)
        @statuscheck FMI.fmi3SetFloat64(
            wrapper.instance, wrapper.input_value_references, inputs
        )
    end
    return wrapper.instance
end

"""
    $(TYPEDSIGNATURES)

Read the FMU's event indicators into `out`, which must have one element per indicator.
"""
function get_fmu_event_indicators!(wrapper::FMI3InstanceWrapper, out)
    @statuscheck FMI.fmi3GetEventIndicators!(
        wrapper.instance, out, Csize_t(length(out))
    )
    return out
end

"""
    $(TYPEDSIGNATURES)

Read the FMU's outputs into `out`, which must have one element per output.
"""
function get_fmu_outputs!(wrapper::FMI3InstanceWrapper, out)
    isempty(out) && return out
    # the count has to be passed: the method without it returns `nothing` instead of a status
    n_outputs = Csize_t(length(out))
    @statuscheck FMI.fmi3GetFloat64!(
        wrapper.instance, wrapper.output_value_references, n_outputs, out, n_outputs
    )
    return out
end

"""
    $(TYPEDSIGNATURES)

Enter Event Mode for the event described by `idx`: a step event if it is `nothing`, which has
no triggered event indicator, and a state event otherwise. DiffEqBase 7.16 always passes one
crossing direction per event indicator, which is what `rootsFound` wants; the scalar branch
handles the index older versions pass, as FMIBase's two `setEventFlags!` methods do.
"""
function enter_fmu_event_mode!(wrapper::FMI3InstanceWrapper, idx)
    n_event_indicators = wrapper.fmu.modelDescription.numberOfEventIndicators
    roots_found = zeros(FMI.fmi3Int32, n_event_indicators)
    step_event = idx === nothing
    if idx isa Integer
        roots_found[idx] = one(FMI.fmi3Int32)
    elseif !step_event
        copyto!(roots_found, idx)
    end
    # FMIImport 1.3.1 wraps the pre-release `fmi3EnterEventMode` signature, which takes the
    # cause of the event; final FMI3 takes none and the extra arguments are ignored
    @statuscheck FMI.fmi3EnterEventMode(
        wrapper.instance, step_event, !step_event, roots_found,
        Csize_t(n_event_indicators), false
    )
    return wrapper.instance
end

"""
    $(TYPEDSIGNATURES)

Enter Event Mode from Step Mode, to handle the event a CoSimulation `fmi3DoStep` reported.
Such an event has no cause the importer could report, so the arguments of the pre-release
`fmi3EnterEventMode` signature FMIImport 1.3.1 wraps are all empty; final FMI3 has none.
"""
function enter_fmu_cs_event_mode!(wrapper::FMI3InstanceWrapper)
    @statuscheck FMI.fmi3EnterEventMode(
        wrapper.instance, false, false, FMI.fmi3Int32[], Csize_t(0), false
    )
    return wrapper.instance
end

"""
    $(TYPEDSIGNATURES)

Leave Event Mode for Continuous-Time Mode. If `values_changed`, the event iteration changed
the continuous states of the FMU, which FMI then requires to be re-read; they are returned in
`wrapper.event_states_buffer`. Returns `nothing` otherwise.
"""
function leave_fmu_event_mode!(wrapper::FMI3InstanceWrapper, values_changed)
    @statuscheck FMI.fmi3EnterContinuousTimeMode(wrapper.instance)
    values_changed || return nothing
    buffer = wrapper.event_states_buffer
    @statuscheck FMI.fmi3GetContinuousStates!(
        wrapper.instance, buffer, Csize_t(length(buffer))
    )
    return buffer
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
        # leaving initialization mode puts the FMU in Event Mode, so the initial event
        # iteration has to converge before Continuous-Time Mode may be entered
        @statuscheck FMI.fmi3ExitInitializationMode(wrapper.instance)
        event_result = do_fmu_event_iteration!(wrapper)
        handle_initial_event_iteration!(wrapper, event_result)
        leave_fmu_event_mode!(wrapper, event_result.values_changed)
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
        # `earlyReturnAllowed` is not ours to choose: FMIImport derives it from the model
        # description, where it is only ever true for an FMU declaring
        # `canReturnEarlyAfterIntermediateUpdate` — which FMIImport 1.3.1 never parses.
        event_mode = cs_event_mode_supported(wrapper.fmu)
        wrapper.instance = FMI.fmi3InstantiateCoSimulation!(
            wrapper.fmu; eventModeUsed = event_mode
        )::FMI.FMU3Instance
        get_instance_common!(wrapper, inputs, params, t, t_end)
        if !isempty(states)
            @statuscheck FMI.fmi3SetFloat64(
                wrapper.instance, wrapper.state_value_references, states
            )
        end
        @statuscheck FMI.fmi3ExitInitializationMode(wrapper.instance)
        if event_mode
            # an instance using event mode leaves initialization mode in Event Mode, so the
            # initial event iteration has to converge before Step Mode may be entered
            event_result = do_fmu_event_iteration!(wrapper)
            handle_initial_event_iteration!(wrapper, event_result)
            @statuscheck FMI.fmi3EnterStepMode(wrapper.instance)
        end
    else
        if !isempty(inputs)
            @statuscheck FMI.fmi3SetFloat64(
                wrapper.instance, wrapper.input_value_references, inputs
            )
        end
        if !isempty(states) && wrapper.states_settable_after_init
            @statuscheck FMI.fmi3SetFloat64(
                wrapper.instance, wrapper.state_value_references, states
            )
        end
    end
    return wrapper.instance
end

"""
    $(TYPEDSIGNATURES)

Call `fmiXCompletedIntegratorStep` with `noSetFMUStatePriorToCurrentPoint` as false. Returns
the FMU's requests to enter event mode and to terminate the simulation.
"""
function partiallyCompleteIntegratorStep(wrapper::FMI3InstanceWrapper)
    enterEventMode = Ref(FMI.fmi3False)
    terminateSimulation = Ref(FMI.fmi3False)
    @statuscheck FMI.fmi3CompletedIntegratorStep!(
        wrapper.instance, FMI.fmi3False, enterEventMode, terminateSimulation
    )
    return (;
        enter_event_mode = enterEventMode[] == FMI.fmi3True,
        terminate = terminateSimulation[] == FMI.fmi3True,
    )
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
    size = (length(states) + length(fn.output_value_references),)
    eltype = eltype(states)
    ndims = 1
end

@register_array_symbolic (fn::FMI3InstanceWrapper)(
    wrapper::FMI3InstanceWrapper, states::Vector{<:Real},
    inputs::Vector{<:Real}, params::Vector{<:Real}, t::Real
) begin
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
        states, inputs, params, t
    )
    instance = get_instance_ME!(wrapper, inputs, params, t)

    states_buffer = wrapper.states_buffer
    outputs_buffer = wrapper.outputs_buffer
    #buffer size matches
    @assert length(states_buffer) == length(states)
    # Defined in FMIBase.jl/src/eval.jl
    # Doesn't seem to be documented, but somehow this is the only way to
    # propagate inputs to the FMU consistently. I have no idea why.
    instance(;
        x = states, u = inputs, u_refs = wrapper.input_value_references,
        p = params, p_refs = wrapper.param_value_references, t = t
    )
    # the outputs have to be read before the derivatives. An OpenModelica FMU evaluates its
    # algebraic and output equations when a variable is read, but reading the derivatives
    # evaluates only the ODE block and still clears the flag saying an evaluation is due
    # (`sources/fmi-export/fmu2_model_interface.c.inc`: `updateIfNeeded` against
    # `internalGetDerivatives`), which leaves a purely-output variable one evaluation stale.
    get_fmu_outputs!(wrapper, outputs_buffer)
    instance(; dx = states_buffer, dx_refs = wrapper.derivative_value_references)
    wrapper.res_buffer[1:length(states_buffer)] .= states_buffer
    wrapper.res_buffer[(length(states_buffer) + 1):end] .= outputs_buffer
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
end

function (fn::FMI2CSFunctor)(wrapper::FMI2InstanceWrapper, states, inputs, params, t)
    states = states isa SubArray ? copy(states) : states
    inputs = inputs isa SubArray ? copy(inputs) : inputs
    params = params isa SubArray ? copy(params) : params
    instance = get_instance_CS!(wrapper, states, inputs, params, t, fn.t_end[])
    if isempty(fn.output_value_references)
        return eltype(states)[]
    else
        return FMI.fmi2GetReal(instance, fn.output_value_references)
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
    if isnan(ctx.t_end[])
        ctx.t_end[] = integrator.sol.prob.tspan[2]
    end
    if wrapper.instance !== nothing
        reset_instance!(wrapper)
    end

    instance = get_instance_CS!(wrapper, states, inputs, params, t, ctx.t_end[])
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

    instance = get_instance_CS!(wrapper, states, inputs, params, integrator.t, ctx.t_end[])
    if !isempty(inputs)
        FMI.fmi2SetReal(
            instance, wrapper.input_value_references, Csize_t(length(inputs)), inputs
        )
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
    """
    `tspan[2]`, set during callback initialization
    """
    t_end::Base.RefValue{FMI.fmi3Float64}
end

function (fn::FMI3CSFunctor)(wrapper::FMI3InstanceWrapper, states, inputs, params, t)
    states = states isa SubArray ? copy(states) : states
    inputs = inputs isa SubArray ? copy(inputs) : inputs
    params = params isa SubArray ? copy(params) : params
    instance = get_instance_CS!(wrapper, states, inputs, params, t, fn.t_end[])

    if isempty(fn.output_value_references)
        return eltype(states)[]
    else
        return FMI.fmi3GetFloat64(instance, fn.output_value_references)
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
    if isnan(ctx.t_end[])
        ctx.t_end[] = integrator.sol.prob.tspan[2]
    end
    if wrapper.instance !== nothing
        reset_instance!(wrapper)
    end
    instance = get_instance_CS!(wrapper, states, inputs, params, t, ctx.t_end[])
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

An FMU instantiated with event mode reports the events it encounters during a communication
step, which the importer has to handle for the FMU to apply them.
"""
function fmiCSStep!(m, o, ctx::FMI3CSFunctor, integrator)
    wrapper = o.wrapper
    states = isdefined(m, :states) ? m.states : ()
    inputs = o.inputs
    params = o.params
    t = o.t
    dt = o.dt

    instance = get_instance_CS!(wrapper, states, inputs, params, integrator.t, ctx.t_end[])
    if !isempty(inputs)
        FMI.fmi3SetFloat64(instance, wrapper.input_value_references, inputs)
    end
    eventEncountered = Ref(FMI.fmi3False)
    terminateSimulation = Ref(FMI.fmi3False)
    earlyReturn = Ref(FMI.fmi3False)
    lastSuccessfulTime = Ref(zero(FMI.fmi3Float64))
    @statuscheck FMI.fmi3DoStep!(
        instance, integrator.t - dt, dt, FMI.fmi3True, eventEncountered,
        terminateSimulation, earlyReturn, lastSuccessfulTime
    )
    @assert earlyReturn[] == FMI.fmi3False
    if terminateSimulation[] == FMI.fmi3True
        SciMLBase.terminate!(integrator)
        return m
    end
    if eventEncountered[] == FMI.fmi3True
        enter_fmu_cs_event_mode!(wrapper)
        event_result = do_fmu_event_iteration!(wrapper)
        if event_result.terminate
            # the instance stays in Event Mode: nothing is stepped after `terminate!`
            SciMLBase.terminate!(integrator)
            return m
        end
        @statuscheck FMI.fmi3EnterStepMode(instance)
        wrapper.next_event_time = event_result.next_event_time
        handle_fmu_time_event!(wrapper, event_result.next_event_time)
    end

    # read after the event iteration, so the outputs carry the post-event values
    if isdefined(m, :states)
        @statuscheck FMI.fmi3GetFloat64!(instance, ctx.state_value_references, m.states)
    end
    if isdefined(m, :outputs)
        @statuscheck FMI.fmi3GetFloat64!(instance, ctx.output_value_references, m.outputs)
    end

    return m
end

end # module
