using SymbolicIndexingInterface
using Setfield
using StaticArrays

struct Input
    var::Num
    data::SVector
    time::SVector
end

function DiffEqBase.solve(prob::SciMLBase.AbstractDEProblem, inputs::Union{Input, Vector{Input}}, args...; input_funs, kwargs...)

    set_input!, finalize! = input_funs

    tstops = Float64[]
    callbacks = DiscreteCallback[]
    if !isa(inputs, Vector)
        inputs = [inputs]
    end

    for input::Input in inputs
        tstops = union(tstops, input.time)
        
        condition = (u,t,integrator) -> any(t .== input.time)
        affect! = function (integrator) 
            i = findfirst(integrator.t .== input.time)
            set_input!(integrator, input.var, input.data[i])
        end
        callback = DiscreteCallback(condition, affect!)

        push!(callbacks, callback)
    end

    # finalize!
    t_end = prob.tspan[2]
    condition = (u,t,integrator) -> (t == t_end)
    affect! = (integrator) -> finalize!(integrator)
    callback = DiscreteCallback(condition, affect!)

    push!(callbacks, callback)
    push!(tstops, t_end)

    return solve(prob, args...; tstops, callback=CallbackSet(callbacks...), kwargs...)
end

function setup_inputs(sys)

    inputs = ModelingToolkit.unbound_inputs(sys)
    setters = Dict{Num, Function}()

    if !isempty(inputs)
        sdcs = ModelingToolkit.SymbolicDiscreteCallback[]
        for x in inputs
            affect = ModelingToolkit.ImperativeAffect((m, o, c, i)->m, modified=(;x))
            sdc = ModelingToolkit.SymbolicDiscreteCallback(Inf, affect)

            push!(sdcs, sdc)
        end

        @set! sys.discrete_events = sdcs
        sys = complete(sys)  # @set! sys.index_cache = ModelingToolkit.IndexCache(sys)

        for (i,x) in enumerate(inputs)
            setter = SymbolicIndexingInterface.setsym(sys, x)
            sdc = sdcs[i]

            setval = function (integrator, set, val=NaN)
                if set
                    println("::setting $x to $val @ $(integrator.t)s")
                    setter(integrator, val)
                else
                    println("::saving $x @ $(integrator.t)s")
                end
                ModelingToolkit.save_callback_discretes!(integrator, sdc)
            end

            setters[x] = setval
        end
    end

    set_input! = function (integrator, var, value) 
        setters[var](integrator, true, value)
        u_modified!(integrator, true)
    end

    finalize! = function (integrator) 
        for ky in keys(setters)
            setters[ky](integrator, false)
        end
    end

    return sys, set_input!, finalize!
end

