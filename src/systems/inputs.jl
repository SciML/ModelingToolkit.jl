using SymbolicIndexingInterface
using Setfield
using StaticArrays
using OrdinaryDiffEqCore

"""
    Input(var, data::Vector{<:Real}, time::Vector{<:Real})

Create an `Input` object that specifies predetermined input values for a variable at specific time points.

# Arguments
- `var`: The symbolic variable (marked with `[input=true]` metadata) to be used as an input.
- `data`: A vector of real values that the input variable should take at the corresponding time points.
- `time`: A vector of time points at which the input values should be applied. Must be the same length as `data`.

# Description
The `Input` struct is used with the extended `solve` method to provide time-varying inputs to a system
during simulation. When passed to `solve(prob, [input1, input2, ...], alg)`, the solver will automatically
set the input variable to the specified values at the specified times using discrete callbacks.

This provides a "determinate form" of input handling where all input values are known a priori,
as opposed to setting inputs manually during integration with [`set_input!`](@ref).

See also [`set_input!`](@ref), [`finalize!`](@ref)
"""
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
    events::Tuple{SymbolicDiscreteCallback}
    vars::Tuple{SymbolicUtils.BasicSymbolic{Real}}
    setters::Tuple{SymbolicIndexingInterface.ParameterHookWrapper}
end

InputFunctions(events::Vector, vars::Vector, setters::Vector) = InputFunctions(Tuple(events), Tuple(vars), Tuple(setters))

"""
    set_input!(integrator, var, value::Real)

Set the value of an input variable during integration.

# Arguments
- `integrator`: An ODE integrator object (from `init(prob, alg)` or available in callbacks).
- `var`: The symbolic input variable to set (must be marked with `[input=true]` metadata and included in the `inputs` keyword of `@mtkcompile`).
- `value`: The new real-valued input to assign to the variable.
- `input_funs` (optional): The `InputFunctions` object associated with the system. If not provided, it will be retrieved from `integrator.f.sys`.

# Description
This function allows you to manually set input values during integration, providing an "indeterminate form"
of input handling where inputs can be computed on-the-fly. This is useful when input values depend on
runtime conditions, external data sources, or interactive user input.

After setting input values with `set_input!`, you must call [`finalize!`](@ref) at the end of integration
to ensure all discrete callbacks are properly saved.

# Example
```julia
@variables x(t) [input=true]
@variables y(t) = 0

eqs = [D(y) ~ x]
@mtkcompile sys = System(eqs, t, [x, y], []) inputs=[x]

prob = ODEProblem(sys, [], (0, 4))
integrator = init(prob, Tsit5())

# Set input and step forward
set_input!(integrator, sys.x, 1.0)
step!(integrator, 1.0, true)

set_input!(integrator, sys.x, 2.0)
step!(integrator, 1.0, true)

# Must call finalize! at the end
finalize!(integrator)
```

See also [`finalize!`](@ref), [`Input`](@ref)
"""
function set_input!(input_funs::InputFunctions, integrator::OrdinaryDiffEqCore.ODEIntegrator, var, value::Real)
    i = findfirst(isequal(var), input_funs.vars)
    setter = input_funs.setters[i]
    event = input_funs.events[i]

    setter(integrator, value)
    save_callback_discretes!(integrator, event)
    u_modified!(integrator, true)
    return nothing
end
set_input!(integrator, var, value::Real) = set_input!(get_input_functions(integrator.f.sys), integrator, var, value)

"""
    finalize!(integrator)

Finalize all input callbacks at the end of integration.

# Arguments
- `integrator`: An ODE integrator object (from `init(prob, alg)` or available in callbacks).
- `input_funs` (optional): The `InputFunctions` object associated with the system. If not provided, it will be retrieved from `integrator.f.sys`.

# Description
This function must be called after using [`set_input!`](@ref) to manually set input values during integration.
It ensures that all discrete callbacks associated with input variables are properly saved in the solution,
making the input values accessible when querying the solution at specific time points.

Without calling `finalize!`, input values set with `set_input!` may not be correctly recorded in the
final solution object, leading to incorrect results when indexing the solution.

See also [`set_input!`](@ref), [`Input`](@ref)
"""
function finalize!(input_funs::InputFunctions, integrator)

    for i in eachindex(input_funs.vars)
        save_callback_discretes!(integrator, input_funs.events[i])
    end

    return nothing
end
finalize!(integrator) = finalize!(get_input_functions(integrator.f.sys), integrator)

(input_funs::InputFunctions)(integrator, var, value::Real) = set_input!(input_funs, integrator, var, value)
(input_funs::InputFunctions)(integrator) = finalize!(input_funs, integrator)

function build_input_functions(sys, inputs)
    
    # Here we ensure the inputs have metadata marking the discrete variables as parameters.  In some
    # cases the inputs can be fed to this function before they are converted to parameters by mtkcompile.
    vars = SymbolicUtils.BasicSymbolic[isparameter(x) ? x : toparam(x) for x in unwrap.(inputs)] 
    setters = []
    events = SymbolicDiscreteCallback[]
    defaults = get_defaults(sys)
    if !isempty(vars)
        
        for x in vars
            affect = ImperativeAffect((m, o, c, i)->m, modified=(;x))
            sdc = SymbolicDiscreteCallback(Inf, affect)

            push!(events, sdc)

            # ensure that the ODEProblem does not complain about missing parameter map
            if !haskey(defaults, x)
                push!(defaults, x => 0.0)
            end

        end

        @set! sys.discrete_events = events
        @set! sys.index_cache = ModelingToolkit.IndexCache(sys)
        @set! sys.defaults = defaults

        setters = [SymbolicIndexingInterface.setsym(sys, x) for x in vars]
    
        @set! sys.input_functions = InputFunctions(events, vars, setters)

    end
    

    return sys 
end




function DiffEqBase.solve(prob::SciMLBase.AbstractDEProblem, inputs::Vector{Input}, args...; kwargs...)

    tstops = Float64[]
    callbacks = DiscreteCallback[]

    # set_input!
    for input::Input in inputs

        tstops = union(tstops, input.time)
        condition = (u,t,integrator) -> any(t .== input.time)
        affect! = function (integrator) 
            @inbounds begin
                i = findfirst(integrator.t .== input.time)
                set_input!(integrator, input.var, input.data[i])
            end 
        end
        push!(callbacks, DiscreteCallback(condition, affect!))

        # DiscreteCallback doesn't hit on t==0, workaround...
        if input.time[1] == 0
            prob.ps[input.var] = input.data[1]
        end
    
    end

    # finalize!
    t_end = prob.tspan[2]
    condition = (u,t,integrator) -> (t == t_end)
    affect! = (integrator) -> finalize!(integrator)
    push!(callbacks, DiscreteCallback(condition, affect!))
    push!(tstops, t_end)

    return solve(prob, args...; tstops, callback=CallbackSet(callbacks...), kwargs...)
end
