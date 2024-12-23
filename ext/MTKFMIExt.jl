module MTKFMIExt

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit as MTK
import FMI

macro statuscheck(expr)
    @assert Meta.isexpr(expr, :call)
    fn = expr.args[1]
    @assert Meta.isexpr(fn, :.)
    @assert fn.args[1] == :FMI
    fnname = fn.args[2]

    instance = expr.args[2]

    return quote
        status = $expr
        fnname = $fnname
        if (status isa Tuple && status[1] == FMI.fmi2True) ||
           (!(status isa Tuple) && status != FMI.fmi2StatusOK &&
            status != FMI.fmi2StatusWarning)
            if status != FMI.fmi2StatusFatal
                FMI.fmi2Terminate(wrapper.instance)
            end
            FMI.fmi2FreeInstance!(wrapper.instance)
            wrapper.instance = nothing
            error("FMU Error: status $status")
        end
    end |> esc
end

function MTK.FMIComponent(::Val{2}, ::Val{:ME}; fmu = nothing, tolerance = 1e-6, name)
    value_references = Dict()
    defs = Dict()
    states = []
    diffvars = []
    observed = Equation[]
    stateT = Float64
    for (valRef, snames) in FMI.getStateValueReferencesAndNames(fmu)
        stateT = FMI.dataTypeForValueReference(fmu, valRef)
        snames = map(parseFMIVariableName, snames)
        vars = [MTK.unwrap(only(@variables $sname(t)::stateT)) for sname in snames]
        for i in eachindex(vars)
            if i == 1
                push!(diffvars, vars[i])
            else
                push!(observed, vars[i] ~ vars[1])
            end
            value_references[vars[i]] = valRef
        end
        append!(states, vars)
    end

    inputs = []
    for (valRef, snames) in FMI.getInputValueReferencesAndNames(fmu)
        snames = map(parseFMIVariableName, snames)
        vars = [MTK.unwrap(only(@variables $sname(t)::stateT)) for sname in snames]
        for i in eachindex(vars)
            if i == 1
                push!(inputs, vars[i])
            else
                push!(observed, vars[i] ~ vars[1])
            end
            value_references[vars[i]] = valRef
        end
        append!(states, vars)
    end

    outputs = []
    for (valRef, snames) in FMI.getOutputValueReferencesAndNames(fmu)
        snames = map(parseFMIVariableName, snames)
        vars = [MTK.unwrap(only(@variables $sname(t)::stateT)) for sname in snames]
        for i in eachindex(vars)
            if i == 1
                push!(outputs, vars[i])
            else
                push!(observed, vars[i] ~ vars[1])
            end
            value_references[vars[i]] = valRef
        end
        append!(states, vars)
    end

    params = []
    parameter_dependencies = Equation[]
    for (valRef, pnames) in FMI.getParameterValueReferencesAndNames(fmu)
        defval = FMI.getStartValue(fmu, valRef)
        T = FMI.dataTypeForValueReference(fmu, valRef)
        pnames = map(parseFMIVariableName, pnames)
        vars = [MTK.unwrap(only(@parameters $pname::T)) for pname in pnames]
        for i in eachindex(vars)
            if i == 1
                push!(params, vars[i])
            else
                push!(parameter_dependencies, vars[i] ~ vars[1])
            end
            value_references[vars[i]] = valRef
        end
        defs[vars[1]] = defval
    end

    input_value_references = UInt32[value_references[var] for var in inputs]
    param_value_references = UInt32[value_references[var] for var in params]
    @parameters wrapper::FMI2InstanceWrapper = FMI2InstanceWrapper(
        fmu, param_value_references, input_value_references, tolerance)

    output_value_references = UInt32[value_references[var] for var in outputs]
    buffer_length = length(diffvars) + length(outputs)
    _functor = FMI2MEFunctor(zeros(buffer_length), output_value_references)
    @parameters (functor::(typeof(_functor)))(..)[1:buffer_length] = _functor
    call_expr = functor(wrapper, copy(diffvars), copy(inputs), copy(params), t)

    diffeqs = Equation[]
    for (i, var) in enumerate([D.(diffvars); outputs])
        push!(diffeqs, var ~ call_expr[i])
    end

    finalize_affect = MTK.FunctionalAffect(fmi2MEFinalize!, [], [wrapper], [])
    step_affect = MTK.FunctionalAffect(fmi2MEStep!, [], [wrapper], [])
    instance_management_callback = MTK.SymbolicDiscreteCallback(
        (t != t - 1), step_affect; finalize = finalize_affect)

    push!(params, wrapper, functor)
    eqs = [observed; diffeqs]
    return ODESystem(eqs, t, states, params; parameter_dependencies, defaults = defs,
        discrete_events = [instance_management_callback], name)
end

function parseFMIVariableName(name::AbstractString)
    return Symbol(replace(name, "." => "__"))
end

mutable struct FMI2InstanceWrapper
    const fmu::FMI.FMU2
    const param_value_references::Vector{UInt32}
    const input_value_references::Vector{UInt32}
    const tolerance::FMI.fmi2Real
    instance::Union{FMI.FMU2Component, Nothing}
end

function FMI2InstanceWrapper(fmu, params, inputs, tolerance)
    FMI2InstanceWrapper(fmu, params, inputs, tolerance, nothing)
end

function get_instance!(wrapper::FMI2InstanceWrapper, states, inputs, params, t)
    if wrapper.instance === nothing
        wrapper.instance = FMI.fmi2Instantiate!(wrapper.fmu)::FMI.FMU2Component
        if !isempty(params)
            @statuscheck FMI.fmi2SetReal(wrapper.instance, wrapper.param_value_references,
                Csize_t(length(wrapper.param_value_references)), params)
        end
        @statuscheck FMI.fmi2SetupExperiment(
            wrapper.instance, FMI.fmi2True, wrapper.tolerance, t, FMI.fmi2False, t)
        @statuscheck FMI.fmi2EnterInitializationMode(wrapper.instance)
        if !isempty(inputs)
            @statuscheck FMI.fmi2SetReal(wrapper.instance, wrapper.input_value_references,
                Csize_t(length(wrapper.param_value_references)), inputs)
        end
        @statuscheck FMI.fmi2ExitInitializationMode(wrapper.instance)
        eventInfo = FMI.fmi2NewDiscreteStates(wrapper.instance)
        @assert eventInfo.newDiscreteStatesNeeded == FMI.fmi2False
        # TODO: Support FMU events
        @statuscheck FMI.fmi2EnterContinuousTimeMode(wrapper.instance)
    end
    instance = wrapper.instance
    @statuscheck FMI.fmi2SetTime(instance, t)
    @statuscheck FMI.fmi2SetContinuousStates(instance, states)
    if !isempty(inputs)
        @statuscheck FMI.fmi2SetReal(instance, wrapper.input_value_references,
            Csize_t(length(wrapper.param_value_references)), inputs)
    end

    return instance
end

function complete_step!(wrapper::FMI2InstanceWrapper)
    wrapper.instance === nothing && return
    @statuscheck FMI.fmi2CompletedIntegratorStep(wrapper.instance, FMI.fmi2True)
end

function reset_instance!(wrapper::FMI2InstanceWrapper)
    wrapper.instance === nothing && return
    FMI.fmi2Terminate(wrapper.instance)
    FMI.fmi2FreeInstance!(wrapper.instance)
    wrapper.instance = nothing
end

struct FMI2MEFunctor{T}
    return_buffer::Vector{T}
    output_value_references::Vector{UInt32}
end

@register_array_symbolic (fn::FMI2MEFunctor)(
    wrapper::FMI2InstanceWrapper, states::Vector{<:Real},
    inputs::Vector{<:Real}, params::Vector{<:Real}, t::Real) begin
    size = (length(states) + length(fn.output_value_references),)
    eltype = eltype(states)
    ndims = 1
end

function (fn::FMI2MEFunctor)(wrapper::FMI2InstanceWrapper, states, inputs, params, t)
    instance = get_instance!(wrapper, states, inputs, params, t)

    states_buffer = zeros(length(states))
    @statuscheck FMI.fmi2GetDerivatives!(instance, states_buffer)
    outputs_buffer = zeros(length(fn.output_value_references))
    FMI.fmi2GetReal!(instance, fn.output_value_references, outputs_buffer)
    return [states_buffer; outputs_buffer]
end

function fmi2MEStep!(integrator, u, p, ctx)
    wrapper_idx = p[1]
    wrapper = integrator.ps[wrapper_idx]
    complete_step!(wrapper)
end

function fmi2MEFinalize!(integrator, u, p, ctx)
    wrapper_idx = p[1]
    wrapper = integrator.ps[wrapper_idx]
    reset_instance!(wrapper)
end

end # module
