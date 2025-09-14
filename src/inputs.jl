using SymbolicIndexingInterface
using Setfield
using StaticArrays

struct Input
    var::Num
    data::SVector
    time::SVector
end

function Input(var, data::Vector{<:Real}, time::Vector{<:Real}) 
    n = length(data)
    return Input(var, SVector{n}(data), SVector{n}(time))
end

struct InputFunctions
    events::Tuple
    vars::Tuple
    setters::Tuple
end

InputFunctions(events::Vector, vars::Vector, setters::Vector) = InputFunctions(Tuple(events), Tuple(vars), Tuple(setters))

function set_input!(input_funs::InputFunctions, integrator, var, value::Real) 
    i = findfirst(isequal(var), input_funs.vars)
    setter = input_funs.setters[i]
    event = input_funs.events[i]
    
    setter(integrator, value)           
    save_callback_discretes!(integrator, event)
    u_modified!(integrator, true)
    return nothing
end

function finalize!(input_funs::InputFunctions, integrator) 

    for i in eachindex(input_funs.vars)
        save_callback_discretes!(integrator, input_funs.events[i])
    end

    return nothing
end

(input_funs::InputFunctions)(integrator, var, value::Real) = set_input!(input_funs, integrator, var, value)
(input_funs::InputFunctions)(integrator) = finalize!(input_funs, integrator)

function setup_inputs(sys, inputs = unbound_inputs(sys))
    
    vars = SymbolicUtils.BasicSymbolic[isparameter(x) ? x : toparam(x) for x in unwrap.(inputs)]
    setters = []
    events = SymbolicDiscreteCallback[]
    if !isempty(vars)
        
        for x in vars
            affect = ImperativeAffect((m, o, c, i)->m, modified=(;x))
            sdc = SymbolicDiscreteCallback(Inf, affect)

            push!(events, sdc)
        end

        @set! sys.discrete_events = events
        @set! sys.index_cache = ModelingToolkit.IndexCache(sys)

        setters = [SymbolicIndexingInterface.setsym(sys, x) for x in vars]
        
    end

    return sys, InputFunctions(events, vars, setters)
end




function DiffEqBase.solve(prob::SciMLBase.AbstractDEProblem, inputs::Vector{Input}, args...; input_funs::InputFunctions, kwargs...)

    tstops = Float64[]
    callbacks = DiscreteCallback[]

    for input::Input in inputs

        tstops = union(tstops, input.time)
        condition = (u,t,integrator) -> any(t .== input.time)
        affect! = function (integrator) 
            @inbounds begin
                i = findfirst(integrator.t .== input.time)
                input_funs(integrator, input.var, input.data[i])
            end 
        end
        push!(callbacks, DiscreteCallback(condition, affect!))
    
    end

    # finalize!
    t_end = prob.tspan[2]
    condition = (u,t,integrator) -> (t == t_end)
    affect! = (integrator) -> input_funs(integrator)
    push!(callbacks, DiscreteCallback(condition, affect!))
    push!(tstops, t_end)

    return solve(prob, args...; tstops, callback=CallbackSet(callbacks...), kwargs...)
end
