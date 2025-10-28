# Building Models with Discrete Data, Interpolations, and Lookup Tables

There are 4 ways to include data as part of a model.

 1. using `ModelingToolkitStandardLibrary.Blocks.Interpolation`
 2. using `ModelingToolkitStandardLibrary.Blocks.ParametrizedInterpolation`
 3. using a custom component with external data (not recommended)
 4. using `ModelingToolkitStandardLibrary.Blocks.SampledData` (legacy)

This tutorial demonstrate each case and explain the pros and cons of each.

## `Interpolation` Block

The `ModelingToolkitStandardLibrary.Blocks.Interpolation` component is easy to use and is performant.
It is similar to using callable parameters, but it provides a block interface with `RealInput` and `RealOutput` connectors.
The `Interpolation` is compatible with interpolation types from `DataInterpolation`.

```@docs
ModelingToolkitStandardLibrary.Blocks.Interpolation
```

Here is an example on how to use it. Let's consider a mass-spring-damper system, where
we have an external force as an input. We then generate some example data in a `DataFrame`
that would represent a measurement of the input. In a more realistic case, this `DataFrame`
would be read from a file.

```@example interpolation_block
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Blocks
using DataInterpolations
using OrdinaryDiffEq
using DataFrames
using Plots

function MassSpringDamper(; name)
    @named input = RealInput()
    @variables f(t) x(t)=0 dx(t)=0 ddx(t)
    @parameters m=10 k=1000 d=1

    eqs = [f ~ input.u
           ddx * 10 ~ k * x + d * dx + f
           D(x) ~ dx
           D(dx) ~ ddx]

    System(eqs, t; name, systems = [input])
end

function MassSpringDamperSystem(data, time; name)
    @named src = Interpolation(LinearInterpolation, data, time)
    @named clk = ContinuousClock()
    @named model = MassSpringDamper()

    eqs = [connect(clk.output, src.input)
           connect(src.output, model.input)]

    System(eqs, t, [], []; name, systems = [src, clk, model])
end

function generate_data()
    dt = 4e-4
    time = 0:dt:0.1
    data = sin.(2 * pi * time * 100)

    return DataFrame(; time, data)
end

df = generate_data() # example data

@named system = MassSpringDamperSystem(df.data, df.time)
sys = mtkcompile(system)
prob = ODEProblem(sys, [], (0, df.time[end]))
sol = solve(prob)
plot(sol)
```

Note that in the case of the `Interpolation` block, the `data` and the `time` act like
structural parameters.

As such, we can also build the interpolation object outside of the model

```@example interpolation_block
my_interpolation = LinearInterpolation(df.data, df.time)

@mtkmodel MassSpringDamperSystem2 begin
    @components begin
        src = Interpolation(itp = my_interpolation)
        clk = ContinuousClock()
        model = MassSpringDamper()
    end
    @equations begin
        connect(src.input, clk.output)
        connect(src.output, model.input)
    end
end;
@mtkcompile sys = MassSpringDamperSystem2()

prob = ODEProblem(sys, [], (0, df.time[end]))
sol = solve(prob, Tsit5())
plot(sol)
```

Note that the interpolation is constructed outside of the model, so we cannot use `remake` to change the
data. For that usecase, see the `ParametrizedInterpolation`.

## `ParametrizedInterpolation` Block

The `ModelingToolkitStandardLibrary.Blocks.ParametrizedInterpolation` component is similar to `Interpolation`, but as the name suggests, it is parametrized by the data, allowing one to change the underlying data without rebuilding the model as the data is represented via vector parameters.
The main advantage of this block over the [`Interpolation`](@ref) one is that one can use it for optimization problems. Currently, this supports forward mode AD via ForwardDiff, but due to the increased flexibility of the types in the component, this is not as fast as the `Interpolation` block,
so it is recommended to use only when the added flexibility is required.

```@docs
ModelingToolkitStandardLibrary.Blocks.ParametrizedInterpolation
```

Here is an example on how to use it

```@example parametrized_interpolation
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Blocks
using DataInterpolations
using OrdinaryDiffEq
using DataFrames
using Plots

function MassSpringDamper(; name)
    @named input = RealInput()
    vars = @variables f(t) x(t)=0 dx(t) [guess = 0] ddx(t)
    pars = @parameters m=10 k=1000 d=1

    eqs = [f ~ input.u
           ddx * 10 ~ k * x + d * dx + f
           D(x) ~ dx
           D(dx) ~ ddx]

    System(eqs, t, vars, pars; name, systems = [input])
end

function MassSpringDamperSystem(data, time; name)
    @named src = ParametrizedInterpolation(LinearInterpolation, data, time)
    @named clk = ContinuousClock()
    @named model = MassSpringDamper()

    eqs = [connect(model.input, src.output)
           connect(clk.output, src.input)]

    System(eqs, t; name, systems = [src, clk, model])
end

function generate_data()
    dt = 4e-4
    time = 0:dt:0.1
    data = sin.(2 * pi * time * 100)

    return DataFrame(; time, data)
end

df = generate_data() # example data

@named system = MassSpringDamperSystem(df.data, df.time)
sys = mtkcompile(system)
prob = ODEProblem(sys, [], (0, df.time[end]))
sol = solve(prob)
plot(sol)
```

If we want to run a new data set, this requires only remaking the problem and solving again

```@example parametrized_interpolation
prob2 = remake(prob, p = [sys.src.data => ones(length(df.data))])
sol2 = solve(prob2)
plot(sol2)
```

!!! note
    
    Note that when changing the data, the length of the new data must be the same as the length of the original data.

## Custom Component with External Data

The below code shows how to include data using a `Ref` and registered `get_sampled_data` function.  This example uses a very basic function which requires non-adaptive solving and sampled data.  As can be seen, the data can easily be set and changed before solving.

```@example custom_component_external_data
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq

const rdata = Ref{Vector{Float64}}()

dt = 4e-4
time = 0:dt:0.1
# Data Sets
data1 = sin.(2 * pi * time * 100)
data2 = cos.(2 * pi * time * 50)

function get_sampled_data(t)
    i = floor(Int, t / dt) + 1
    x = rdata[][i]

    return x
end

Symbolics.@register_symbolic get_sampled_data(t)

function CustomSystem(; name)
    vars = @variables f(t)=0 x(t)=0 dx(t)=0 ddx(t)=0
    pars = @parameters m=10 k=1000 d=1

    eqs = [f ~ get_sampled_data(t)
           ddx * 10 ~ k * x + d * dx + f
           D(x) ~ dx
           D(dx) ~ ddx]

    System(eqs, t, vars, pars; name)
end

@named system = CustomSystem()
sys = mtkcompile(system)
prob = ODEProblem(sys, [], (0, time[end]))

rdata[] = data1
sol1 = solve(prob, ImplicitEuler(); dt, adaptive = false)
ddx1 = sol1[sys.ddx]

rdata[] = data2
sol2 = solve(prob, ImplicitEuler(); dt, adaptive = false)
ddx2 = sol2[sys.ddx]
```

The drawback of this method is that the solution observables can be linked to the data `Ref`, which means that if the data changes then the observables are no longer valid.  In this case `ddx` is an observable that is derived directly from the data.  Therefore, `sol1[sys.ddx]` is no longer correct after the data is changed for `sol2`.

```julia
# the following test will fail
@test all(ddx1 .== sol1[sys.ddx]) #returns false
```

Additional code could be added to resolve this issue, for example by using a `Ref{Dict}` that could link a parameter of the model to the data source.  This would also be necessary for parallel processing.

## `SampledData` Component

To resolve the issues presented above, the `ModelingToolkitStandardLibrary.Blocks.SampledData` component can be used which allows for a resusable `System` and self contained data which ensures a solution which remains valid for it's lifetime.  Now it's possible to also parallelize the call to `solve()`.

```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq

function SampledDataSystem(; name)
    src = SampledData(Float64, name = :src)

    vars = @variables f(t)=0 x(t)=0 dx(t)=0 ddx(t)=0
    pars = @parameters m=10 k=1000 d=1

    eqs = [f ~ src.output.u
           ddx * 10 ~ k * x + d * dx + f
           D(x) ~ dx
           D(dx) ~ ddx]

    System(eqs, t, vars, pars; systems = [src], name)
end

@named system = SampledDataSystem()
sys = mtkcompile(system, split = false)
s = complete(system)

dt = 4e-4
time = 0:dt:0.1
data1 = sin.(2 * pi * time * 100)
data2 = cos.(2 * pi * time * 50)

prob = ODEProblem(sys, [], (0, time[end]); split = false, tofloat = false, use_union = true)
defs = ModelingToolkit.defaults(sys)

function get_prob(data)
    defs[s.src.buffer] = Parameter(data, dt)
    # ensure p is a uniform type of Vector{Parameter{Float64}} (converting from Vector{Any})
    p = Parameter.(ModelingToolkit.varmap_to_vars(defs, parameters(sys); tofloat = false))
    remake(prob; p, build_initializeprob = false)
end

prob1 = get_prob(data1)
prob2 = get_prob(data2)

sol1 = Ref{ODESolution}()
sol2 = Ref{ODESolution}()
@sync begin
    @async sol1[] = solve(prob1, ImplicitEuler())
    @async sol2[] = solve(prob2, ImplicitEuler())
end
```

Note, in the above example, we can build the system with an empty `SampledData` component, only setting the expected data type: `@named src = SampledData(Float64)`.  It's also possible to initialize the component with real sampled data: `@named src = SampledData(data, dt)`.  Additionally note that before running an `ODEProblem` using the `SampledData` component, one must be careful about the parameter vector Type.  The `SampledData` component contains a `buffer` parameter of type `Parameter`, therefore we must generate the problem using `tofloat=false`.  This will initially give a parameter vector of type `Vector{Any}` with a mix of numbers and `Parameter` type.  We can convert the vector to a uniform `Parameter` type by running `p = Parameter.(p)`.  This will wrap all the single values in a `Parameter` which will be mathematically equivalent to a `Number`.
