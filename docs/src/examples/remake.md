# Optimizing through an ODE solve and re-creating MTK Problems

Solving an ODE as part of an `OptimizationProblem`'s loss function is a common scenario.
In this example, we will go through an efficient way to model such scenarios using
ModelingToolkit.jl.

First, we build the ODE to be solved. For this example, we will use a Lotka-Volterra model:

```@example Remake
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters α β γ δ
@variables x(t) y(t)
eqs = [D(x) ~ (α - β * y) * x
       D(y) ~ (δ * x - γ) * y]
@mtkbuild odesys = ODESystem(eqs, t)
```

To create the "data" for optimization, we will solve the system with a known set of
parameters.

```@example Remake
using OrdinaryDiffEq

odeprob = ODEProblem(
    odesys, [x => 1.0, y => 1.0], (0.0, 10.0), [α => 1.5, β => 1.0, γ => 3.0, δ => 1.0])
timesteps = 0.0:0.1:10.0
sol = solve(odeprob, Tsit5(); saveat = timesteps)
data = Array(sol)
# add some random noise
data = data + 0.01 * randn(size(data))
```

Now we will create the loss function for the Optimization solve. This will require creating
an `ODEProblem` with the parameter values passed to the loss function. Creating a new
`ODEProblem` is expensive and requires differentiating through the code generation process.
This can be bug-prone and is unnecessary. Instead, we will leverage the `remake` function.
This allows creating a copy of an existing problem with updating state/parameter values. It
should be noted that the types of the values passed to the loss function may not agree with
the types stored in the existing `ODEProblem`. Thus, we cannot use `setp` to modify the
problem in-place. Here, we will use the `replace` function from SciMLStructures.jl since
it allows updating the entire `Tunable` portion of the parameter object which contains the
parameters to optimize.

```@example Remake
using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!

function loss(x, p)
    odeprob = p[1] # ODEProblem stored as parameters to avoid using global variables
    ps = parameter_values(odeprob) # obtain the parameter object from the problem
    ps = replace(Tunable(), ps, x) # create a copy with the values passed to the loss function
    # remake the problem, passing in our new parameter object
    newprob = remake(odeprob; p = ps)
    timesteps = p[2]
    sol = solve(newprob, AutoTsit5(Rosenbrock23()); saveat = timesteps)
    truth = p[3]
    data = Array(sol)
    return sum((truth .- data) .^ 2) / length(truth)
end
```

Note how the problem, timesteps and true data are stored as model parameters. This helps
avoid referencing global variables in the function, which would slow it down significantly.

We could have done the same thing by passing `remake` a map of parameter values. For example,
let us enforce that the order of ODE parameters in `x` is `[α β γ δ]`. Then, we could have
done:

```julia
remake(odeprob; p = [α => x[1], β => x[2], γ => x[3], δ => x[4]])
```

However, passing a symbolic map to `remake` is significantly slower than passing it a
parameter object directly. Thus, we use `replace` to speed up the process. In general,
`remake` is the most flexible method, but the flexibility comes at a cost of performance.

We can perform the optimization as below:

```@example Remake
using Optimization
using OptimizationOptimJL

# manually create an OptimizationFunction to ensure usage of `ForwardDiff`, which will
# require changing the types of parameters from `Float64` to `ForwardDiff.Dual`
optfn = OptimizationFunction(loss, Optimization.AutoForwardDiff())
# parameter object is a tuple, to store differently typed objects together
optprob = OptimizationProblem(
    optfn, rand(4), (odeprob, timesteps, data), lb = 0.1zeros(4), ub = 3ones(4))
sol = solve(optprob, BFGS())
```

To identify which values correspond to which parameters, we can `replace!` them into the
`ODEProblem`:

```@example Remake
replace!(Tunable(), parameter_values(odeprob), sol.u)
odeprob.ps[[α, β, γ, δ]]
```

`replace!` operates in-place, so the values being replaced must be of the same type as those
stored in the parameter object, or convertible to that type. For demonstration purposes, we
can construct a loss function that uses `replace!`, and calculate gradients using
`AutoFiniteDiff` rather than `AutoForwardDiff`.

```@example Remake
function loss2(x, p)
    odeprob = p[1] # ODEProblem stored as parameters to avoid using global variables
    newprob = remake(odeprob) # copy the problem with `remake`
    # update the parameter values in-place
    replace!(Tunable(), parameter_values(newprob), x)
    timesteps = p[2]
    sol = solve(newprob, AutoTsit5(Rosenbrock23()); saveat = timesteps)
    truth = p[3]
    data = Array(sol)
    return sum((truth .- data) .^ 2) / length(truth)
end

# use finite-differencing to calculate derivatives
optfn2 = OptimizationFunction(loss2, Optimization.AutoFiniteDiff())
optprob2 = OptimizationProblem(
    optfn2, rand(4), (odeprob, timesteps, data), lb = 0.1zeros(4), ub = 3ones(4))
sol = solve(optprob2, BFGS())
```

# Re-creating the problem

There are multiple ways to re-create a problem with new state/parameter values. We will go
over the various methods, listing their use cases.

## Pure `remake`

This method is the most generic. It can handle symbolic maps, initializations of
parameters/states dependent on each other and partial updates. However, this comes at the
cost of performance. `remake` is also not always inferable.

## `remake` and `setp`/`setu`

Calling `remake(prob)` creates a copy of the existing problem. This new problem has the
exact same types as the original one, and the `remake` call is fully inferred.
State/parameter values can be modified after the copy by using `setp` and/or `setu`. This
is most appropriate when the types of state/parameter values does not need to be changed,
only their values.

## `replace` and `remake`

`replace` returns a copy of a parameter object, with the appropriate portion replaced by new
values. This is useful for changing the type of an entire portion, such as during the
optimization process described above. `remake` is used in this case to create a copy of the
problem with updated state/unknown values.

## `remake` and `replace!`

`replace!` is similar to `replace`, except that it operates in-place. This means that the
parameter values must be of the same types. This is useful for cases where bulk parameter
replacement is required without needing to change types. For example, optimization methods
where the gradient is not computed using dual numbers (as demonstrated above).
