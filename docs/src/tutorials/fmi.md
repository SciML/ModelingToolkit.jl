```@meta
Draft = true
```

# Importing FMUs

ModelingToolkit is able to import FMUs following the [FMI Standard](https://fmi-standard.org/) versions 2 and 3.
This integration is done through [FMI.jl](https://github.com/ThummeTo/FMI.jl) and requires importing it to
enable the relevant functionality. Currently Model Exchange (ME) and CoSimulation (CS) FMUs are supported.
Model Exchange supports time and state events, accepted-step events, and discrete input changes at accepted steps.
Numeric scalar ports include floating-point, Boolean, enumeration and integer types.
Enumerations must be exactly representable as `Float64`. FMI 3 `Int64` and `UInt64` ports and
parameters support the inclusive ranges `[-2^53, 2^53]` and `[0, 2^53]`, respectively. Supply
integer values so the range check can inspect the value before floating-point conversion.
Array variables are unsupported.

!!! danger "Experimental"

    This functionality is currently experimental and subject to change without a breaking release of
    ModelingToolkit.jl.

## FMUs of full models

Here, we will demonstrate the usage of an FMU of an entire model (as opposed to a single component).
First, the required libraries must be imported and the FMU loaded using FMI.jl.

```@example fmi
using ModelingToolkit, FMI, FMIZoo, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit

# This is a spring-pendulum FMU from FMIZoo.jl. It is a v2 FMU
# and we are importing it in ModelExchange format.
fmu = loadFMU("SpringPendulum1D", "Dymola", "2022x"; type = :ME)
```

Following are the variables in the FMU (both states and parameters):

```@example fmi
fmu.modelDescription.modelVariables
```

Next, [`FMIComponent`](@ref ModelingToolkit.FMIComponent) is used to import the FMU as an MTK component. We provide the FMI
major version as a `Val` to the constructor, along with the loaded FMU and the type as keyword
arguments.

```@example fmi
@named model = ModelingToolkit.FMIComponent(Val(2); fmu, type = :ME)
```

Note how hierarchical names in the FMU (e.g. `mass.m` or `spring.f`) are turned into flattened
names, with `__` being the namespace separator (`mass__m` and `spring__f`).

!!! note

    Eventually we plan to reconstruct a hierarchical system structure mirroring the one indicated
    by the variables in the FMU. This would allow accessing the above mentioned variables as
    `model.mass.m` and `model.spring.f` instead of `model.mass__m` and `model.spring__f` respectively.

Derivative variables such as `der(mass.v)` use the dummy derivative notation, and are hence transformed
into a form similar to `mass__vˍt`. However, they can still be referred to as `D(model.mass__v)`.

```@example fmi
equations(model)
```

Since the FMI spec allows multiple names to alias the same quantity, ModelingToolkit.jl creates
equations to alias them. For example, it can be seen above that `der(mass.v)` and `mass.a` have the
same reference, and hence refer to the same quantity. Correspondingly, there is an equation
`mass__vˍt(t) ~ mass__a(t)` in the system.

!!! note

    Any variables and/or parameters that are not part of the FMU should be ignored, as ModelingToolkit
    creates them to manage the FMU. Unexpected usage of these variables/parameters can lead to errors.

```@example fmi
initial_conditions(model)
```

FMU parameters default to their start values. Supply continuous initial state values in the problem map or
determine them with `initialization_eqs` on the enclosing system or problem. Alternatively, set
`fmu_initial_states = true` on the component to obtain initial states from the FMU. Problem-map values
take precedence over these defaults. Use the default `fmu_initial_states = false` when explicit
initialization equations determine FMU states. A Co-Simulation FMU accepts writes only for states with
a start value and computes states declared with `initial = "calculated"` itself.

We can simulate this model like any other ModelingToolkit system.

```@repl fmi
sys = mtkcompile(model)
prob = ODEProblem(sys, [sys.mass__s => 0.5, sys.mass__v => 0.0], (0.0, 5.0))
sol = solve(prob, Tsit5(); saveat = 0.0:0.1:5.0)
```

FMU state variables support interpolation within the solved interval. When the FMU stays out of
Event Mode during integration and no supplied solver callback affect runs, continuous outputs can
also be evaluated at arbitrary times. A Model Exchange event or an affect from an `ODEProblem`
`callback` activates replay of continuous-output snapshots. The solution records both saved sides
of the affect, including input changes from ordinary parameters, even when the FMU stays out of
Event Mode. Use `continuity = :left` or `continuity = :right` to select the value at an event time.
A continuous-output query without a snapshot raises an error because the historical FMI state or
inputs may be unavailable. Pass `saveat` to `solve` to record outputs at additional times.

```@repl fmi
sol(0.0:0.1:1.0; idxs = sys.mass__a)
```

FMUs following version 3 of the specification can be simulated with almost the same process. This time,
we will create a model from a CoSimulation FMU.

```@example fmi
fmu = loadFMU("SpringPendulum1D", "Dymola", "2023x", "3.0"; type = :CS)
@named inner = ModelingToolkit.FMIComponent(
    Val(3); fmu, communication_step_size = 0.001, type = :CS,
    reinitializealg = BrownFullBasicInit())
nothing # hide
```

This FMU has fewer equations, partly due to missing aliasing variables and partly due to being a CS FMU.
CoSimulation FMUs are bundled with an integrator. As such, they do not function like ME FMUs. Instead,
a callback steps the FMU at periodic intervals in time and obtains the updated state. This state is held
constant until the next time the callback triggers. The periodic interval must be specified through the
`communication_step_size` keyword argument. A smaller step size typically leads to less error but is
more computationally expensive. Variable-step Co-Simulation FMUs communicate at the solve
endpoint. Fixed-step FMUs communicate full steps through the last communication-grid point
at or before the endpoint; a later endpoint uses the held sample from that grid point,
including when it differs only by roundoff. The
`stop_time` keyword sets the FMU experiment limit. The default covers the solve endpoint
and, for fixed-step FMUs, also covers the last native request if adding the configured step
size to the preceding communication point rounds past that endpoint. An explicit
`stop_time` is a hard limit. It must cover the solve endpoint for variable-step FMUs or
the last native request for fixed-step FMUs. If rounding the final variable communication
interval would cross the limit, the native step is rounded down within floating-point
tolerance and the host still reports the scheduled endpoint. A fixed-step solve can
continue past the limit on the held sample. If a scheduled fixed step would exceed an
explicit limit, initialization reports the minimum usable `stop_time`.

The FMU holds the inputs sampled at the left communication point for the entire interval.
States and outputs, including direct-feedthrough outputs, are read at the right communication
point with that same input. The host then samples inputs for the next interval after any
coupled reinitialization. FMI 3 Co-Simulation checks that the FMU reached each requested time.

This model alone does not have any differential variables, and calling `mtkcompile` will lead
to an `System` with no unknowns.

```@example fmi
mtkcompile(inner)
```

Simulating this model will cause the OrdinaryDiffEq integrator to immediately finish, and will not
trigger the callback. Thus, we wrap this system in a trivial system with a differential variable.

```@example fmi
@variables x(t) = 1.0
@mtkcompile sys = System([D(x) ~ x], t; systems = [inner])
```

We can now simulate `sys`.

```@example fmi
prob = ODEProblem(sys, [sys.inner.mass__s => 0.5, sys.inner.mass__v => 0.0], (0.0, 5.0))
sol = solve(prob, Tsit5())
```

The variables of the FMU are discrete, and their timeseries can be obtained at intervals of
`communication_step_size`.

```@example fmi
sol[sys.inner.mass__s]
```

## FMUs of components

FMUs can also be imported as individual components. For this example, we will use custom FMUs used
in the test suite of ModelingToolkit.jl.

```@example fmi
fmu = loadFMU(
    joinpath(@__DIR__, "..", "..", "..", "test", "fmi", "fmus", "SimpleAdder.fmu");
    type = :ME)
fmu.modelDescription.modelVariables
```

For this FMU, `fmu_initial_states = true` uses its built-in initial state `c = 1`.

```@example fmi
@named initial_adder = ModelingToolkit.FMIComponent(Val(2); fmu, type = :ME, fmu_initial_states = true)
@mtkcompile initialized = System(
    [initial_adder.a ~ 1.0, initial_adder.b ~ 2.0], t; systems = [initial_adder])
initialized_sol = solve(ODEProblem(initialized, [], (0.0, 0.1)), Tsit5())
@assert initialized_sol[initialized.initial_adder.c][1] == 1.0
initialized_sol
```

This FMU is equivalent to the following model:

```julia
@mtkmodel SimpleAdder begin
    @variables begin
        a(t)
        b(t)
        c(t)
        out(t)
        out2(t)
    end
    @parameters begin
        value = 1.0
    end
    @equations begin
        out ~ a + b + value
        D(c) ~ out
        out2 ~ 2c
    end
end
```

`a` and `b` are inputs, `c` is a state, and `out` and `out2` are outputs of the component.

```@repl fmi
@named adder = ModelingToolkit.FMIComponent(
    Val(2); fmu, type = :ME, reinitializealg = BrownFullBasicInit());
isinput(adder.a)
isinput(adder.b)
isoutput(adder.out)
isoutput(adder.out2)
```

ModelingToolkit recognizes input and output variables of the component, and attaches the appropriate
metadata. We can now use this component as a subcomponent of a larger system.

```@repl fmi
@variables a(t) b(t) c(t) [guess = 1.0];
@mtkcompile sys = System(
    [adder.a ~ a, adder.b ~ b, D(a) ~ t,
        D(b) ~ adder.out + adder.c, c^2 ~ adder.out + adder.value],
    t;
    systems = [adder])
equations(sys)
```

Note how the output `adder.out` is used in an algebraic equation of the system. We have also given
`sys.c` a guess, expecting it to be solved for by initialization. ModelingToolkit is able to use
FMUs in initialization to solve for initial states. For solvers requiring derivatives through a Model Exchange FMU, explicitly select finite differences.

```@example fmi
prob = ODEProblem(
    sys, [sys.adder.c => 2.0, sys.a => 1.0, sys.b => 1.0, sys.adder.value => 2.0],
    (0.0, 1.0))
solve(prob, Rodas5P(autodiff = AutoFiniteDiff()))
```

CoSimulation FMUs follow a nearly identical process. Since CoSimulation FMUs operate using callbacks,
after triggering the callbacks and altering the discrete state the algebraic equations may no longer
be satisfied. To resolve for the values of algebraic variables, we use the `reinitializealg` keyword
of `FMIComponent`. This is a DAE initialization algorithm to use at the end of every callback. Since
CoSimulation FMUs are not directly involved in the RHS of the system - instead operating through
callbacks - we can use a solver with automatic differentiation.

```@example fmi
fmu = loadFMU(
    joinpath(@__DIR__, "..", "..", "..", "test", "fmi", "fmus", "SimpleAdder.fmu");
    type = :CS)
@named adder = ModelingToolkit.FMIComponent(
    Val(2); fmu, type = :CS, communication_step_size = 1e-3,
    reinitializealg = BrownFullBasicInit())
@mtkcompile sys = System(
    [adder.a ~ a, adder.b ~ b, D(a) ~ t,
        D(b) ~ adder.out + adder.c, c^2 ~ adder.out + adder.value],
    t;
    systems = [adder])
prob = ODEProblem(
    sys, [sys.adder.c => 2.0, sys.a => 1.0, sys.b => 1.0, sys.adder.value => 2.0],
    (0.0, 1.0))
solve(prob, Rodas5P())
```
