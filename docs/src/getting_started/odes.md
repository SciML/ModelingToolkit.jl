# [Building ODEs and DAEs with ModelingToolkit.jl](@id getting_started_ode)

This is an introductory tutorial for ModelingToolkit.jl (MTK). We will demonstrate the
basics of the package by demontrating how to build systems of Ordinary Differential
Equations (ODEs) and Differential-Algebraic Equations (DAEs).

## Installing ModelingToolkit

To install ModelingToolkit, use the Julia package manager. This can be done as follows:

```julia
using Pkg
Pkg.add("ModelingToolkit")
```

## The end goal

TODO

## Basics of MTK

ModelingToolkit.jl is a symbolic-numeric system. This means it allows specifying a model
(such as an ODE) in a similar way to how it would be written on paper. Let's start with a
simple example. The system to be modeled is a first-order lag element:

```math
\dot{x} = \frac{f(t) - x(t)}{\tau}
```

Here, ``t`` is the independent variable (time), ``x(t)`` is the (scalar) unknown
variable, ``f(t)`` is an external forcing function, and ``\tau`` is a
parameter.

For simplicity, we will start off by setting the forcing function to a constant `1`. Every
ODE has a single independent variable. MTK has a common definition for time `t` and the
derivative with respect to it.

```@example ode2
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
```

Next, we declare the (dependent) variables and the parameters of our model:

```@example ode2
@variables x(t)
@parameters τ
```

Note the syntax `x(t)`. We must declare that the variable `x` is a function of the independent
variable `t`. Next, we define the equations of the system:

```@example ode2
eqs = [D(x) ~ (1 - x) / τ]
```

Since `=` is reserved as the assignment operator, MTK uses `~` to denote equality between
expressions. Now we must consolidate all of this information about our system of ODEs into
ModelingToolkit's `System` type.

```@example ode2
sys = System(eqs, t, [x], [τ]; name = :sys)
```

The `System` constructor accepts a `Vector{Equation}` as the first argument, followed by the
independent variable, a list of dependent variables, and a list of parameters. Every system
must be given a name via the `name` keyword argument. Most of the time, we want to name our
system the same as the variable it is assigned to. The `@named` macro helps with this:

```@example ode2
@named sys = System(eqs, t, [x], [τ])
```

Additionally, it may become inconvenient to specify all variables and parameters every time
a system is created. MTK allows omitting these arguments, and will automatically infer them
from the equations.

```@example ode2
@named sys = System(eqs, t)
```

Our system is not quite ready for simulation yet. First, we must use the `mtkcompile`
function which transforms the system into a form that MTK can handle. For our trivial
system, this does not do much.

```@example ode2
simp_sys = mtkcompile(sys)
```

Since building and simplifying a system is a common workflow, MTK provides the `@mtkcompile`
macro for convenience.

```@example ode2
@mtkcompile sys = System(eqs, t)
```

We can now build an `ODEProblem` from the system. ModelingToolkit generates the necessary
code for numerical ODE solvers to solve this system. We need to provide an initial
value for the variable `x` and a value for the parameter `p`, as well as the time span
for which to simulate the system.

```@example ode2
prob = ODEProblem(sys, [x => 0.0, τ => 3.0], (0.0, 10.0))
```

Here, we are saying that `x` should start at `0.0`, `τ` should be `3.0` and the system
should be simulated from `t = 0.0` to `t = 10.0`. To solve the system, we must import a
solver.

```@example ode2
using OrdinaryDiffEq

sol = solve(prob)
```

[OrdinaryDiffEq.jl](https://docs.sciml.ai/DiffEqDocs/stable/) contains a large number of
numerical solvers. It also comes with a default solver which is used when calling
`solve(prob)` and is capable of handling a large variety of systems.

We can obtain the timeseries of `x` by indexing the solution with the symbolic variable:

```@example ode2
sol[x]
```

We can even obtain timeseries of complicated expressions involving the symbolic variables in
the model

```@example ode2
sol[(1 - x) / τ]
```

Perhaps more interesting is a plot of the solution. This can easily be achieved using Plots.jl.

```@example ode2
using Plots

plot(sol)
```

Similarly, we can plot different expressions:

```@example ode2
plot(sol; idxs = (1 - x) / τ)
```
