# Initialization of ODESystems

While for simple numerical ODEs choosing an initial condition can be an easy
affair, with ModelingToolkit's more general differential-algebraic equation
(DAE) system there is more care needed due to the flexibility of the solver
state. In this tutorial we will walk through the functionality involved in
initialization of ODESystem and the diagnostics to better understand and
debug the initialization problem.

## Primer on Initialization of Differential-Algebraic Equations

Before getting started, let's do a brief walkthrough of the mathematical
principles of initialization of DAE systems. Take a DAE written in semi-explicit
form:

```math
\begin{aligned}
    x^\prime &= f(x,y,t) \\
    0 &= g(x,y,t)
\end{aligned}
```

where ``x`` are the differential variables and ``y`` are the algebraic variables.
An initial condition ``u0 = [x(t_0) y(t_0)]`` is said to be consistent if
``g(x(t_0),y(t_0),t_0) = 0``.

For ODEs, this is trivially satisfied. However, for more complicated systems it may
not be easy to know how to choose the variables such that all of the conditions
are satisfied. This is even more complicated when taking into account ModelingToolkit's
simplification engine, given that variables can be eliminated and equations can be
changed. If this happens, how do you know how to initialize the system?

## Initialization By Example: The Cartesian Pendulum

To illustrate how to perform the initialization, let's take a look at the Cartesian
pendulum:

```@example init
using ModelingToolkit, OrdinaryDiffEq, Plots
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters g
@variables x(t) y(t) [state_priority = 10] λ(t)
eqs = [D(D(x)) ~ λ * x
       D(D(y)) ~ λ * y - g
       x^2 + y^2 ~ 1]
@mtkbuild pend = ODESystem(eqs, t)
```

While we defined the system using second derivatives and a length constraint,
the structural simplification system improved the numerics of the system to
be solvable using the dummy derivative technique, which results in 3 algebraic
equations and 2 differential equations. In this case, the differential equations
with respect to `y` and `D(y)`, though it could have just as easily have been
`x` and `D(x)`. How do you initialize such a system if you don't know in advance
what variables may defined the equation's state?

To see how the system works, let's start the pendulum in the far right position,
i.e. `x(0) = 1` and `y(0) = 0`. We can do this by:

```@example init
prob = ODEProblem(pend, [x => 1, y => 0], (0.0, 1.5), [g => 1], guesses = [λ => 1])
```

This solves via:

```@example init
sol = solve(prob, Rodas5P())
plot(sol, idxs = (x, y))
```

and we can check it satisfies our conditions via:

```@example init
conditions = getfield.(equations(pend)[3:end], :rhs)
```

```@example init
[sol[conditions][1]; sol[x][1] - 1; sol[y][1]]
```

Notice that we set `[x => 1, y => 0]` as our initial conditions and `[λ => 1]` as our guess.
The difference is that the initial conditions are **required to be satisfied**, while the
guesses are simply a guess for what the initial value might be. Every variable must have
either an initial condition or a guess, and thus since we did not know what `λ` would be
we set it to 1 and let the initialization scheme find the correct value for λ. Indeed,
the value for `λ` at the initial time is not 1:

```@example init
sol[λ][1]
```

We can similarly choose `λ = 0` and solve for `y` to start the system:

```@example init
prob = ODEProblem(pend, [x => 1, λ => 0], (0.0, 1.5), [g => 1], guesses = [y => 1])
sol = solve(prob, Rodas5P())
plot(sol, idxs = (x, y))
```

or choose to satisfy derivative conditions:

```@example init
prob = ODEProblem(
    pend, [x => 1, D(y) => 0], (0.0, 1.5), [g => 1], guesses = [λ => 0, y => 1])
sol = solve(prob, Rodas5P())
plot(sol, idxs = (x, y))
```

Notice that since a derivative condition is given, we are required to give a
guess for `y`.

We can also directly give equations to be satisfied at the initial point by using
the `initialization_eqs` keyword argument, for example:

```@example init
prob = ODEProblem(pend, [x => 1], (0.0, 1.5), [g => 1], guesses = [λ => 0, y => 1],
    initialization_eqs = [y ~ 0])
sol = solve(prob, Rodas5P())
plot(sol, idxs = (x, y))
```

Additionally, note that the initial conditions are allowed to be functions of other
variables and parameters:

```@example init
prob = ODEProblem(
    pend, [x => 1, D(y) => g], (0.0, 3.0), [g => 1], guesses = [λ => 0, y => 1])
sol = solve(prob, Rodas5P())
plot(sol, idxs = (x, y))
```

## Determinability: Underdetermined and Overdetermined Systems

For this system we have 3 conditions to satisfy:

```@example init
conditions = getfield.(equations(pend)[3:end], :rhs)
```

when we initialize with

```@example init
prob = ODEProblem(pend, [x => 1, y => 0], (0.0, 1.5), [g => 1], guesses = [y => 0, λ => 1])
```

we have two extra conditions to satisfy, `x ~ 1` and `y ~ 0` at the initial point. That gives
5 equations for 5 variables and thus the system is well-formed. What happens if that's not the
case?

```@example init
prob = ODEProblem(pend, [x => 1], (0.0, 1.5), [g => 1], guesses = [y => 0, λ => 1])
```

Here we have 4 equations for 5 unknowns (note: the warning is post-simplification of the
nonlinear system, which solves the trivial `x ~ 1` equation analytical and thus says
3 equations for 4 unknowns). This warning thus lets you know the system is underdetermined
and thus the solution is not necessarily unique. It can still be solved:

```@example init
sol = solve(prob, Rodas5P())
plot(sol, idxs = (x, y))
```

and the found initial condition satisfies all constraints which were given. In the opposite
direction, we may have an overdetermined system:

```@example init
prob = ODEProblem(
    pend, [x => 1, y => 0.0, D(y) => 0], (0.0, 1.5), [g => 1], guesses = [λ => 1])
```

Can that be solved?

```@example init
sol = solve(prob, Rodas5P())
plot(sol, idxs = (x, y))
```

Indeed since we saw `D(y) = 0` at the initial point above, it turns out that this solution
is solvable with the chosen initial conditions. However, for overdetermined systems we often
aren't that lucky. If the set of initial conditions cannot be satisfied, then you will get
a `SciMLBase.ReturnCode.InitialFailure`:

```@example init
prob = ODEProblem(
    pend, [x => 1, y => 0.0, D(y) => 2.0, λ => 1], (0.0, 1.5), [g => 1], guesses = [λ => 1])
sol = solve(prob, Rodas5P())
```

What this means is that the initial condition finder failed to find an initial condition.
While this can be sometimes due to numerical error (which is then helped by picking guesses closer
to the correct value), most circumstances of this come from ill-formed models. Especially
**if your system is overdetermined and you receive an InitialFailure, the initial conditions
may not be analytically satisfiable!**. In our case here, if you sit down with a pen and paper
long enough you will see that `λ = 0` is required for this equation, but since we chose
`λ = 1` we end up with a set of equations that are impossible to satisfy.

!!! note
    
    If you would prefer to have an error instead of a warning in the context of non-fully
    determined systems, pass the keyword argument `fully_determined = true` into the
    problem constructor. Additionally, any warning about not being fully determined can
    be suppressed via passing `warn_initialize_determined = false`.

## Initialization of parameters

Parameters may also be treated as unknowns in the initialization system. Doing so works
almost identically to the standard case. For a parameter to be an initialization unknown
(henceforth referred to as "solved parameter") it must represent a floating point number
(have a `symtype` of `Real` or `<:AbstractFloat`) or an array of such numbers. Additionally,
it must have a guess and one of the following conditions must be satisfied:

 1. The value of the parameter as passed to `ODEProblem` is an expression involving other
    variables/parameters. For example, if `[p => 2q + x]` is passed to `ODEProblem`. In
    this case, `p ~ 2q + x` is used as an equation during initialization.
 2. The parameter has a default (and no value for it is given to `ODEProblem`, since
    that is condition 1). The default will be used as an equation during initialization.
 3. The parameter has a default of `missing`. If `ODEProblem` is given a value for this
    parameter, it is used as an equation during initialization (whether the value is an
    expression or not).
 4. `ODEProblem` is given a value of `missing` for the parameter. If the parameter has a
    default, it will be used as an equation during initialization.

All parameter dependencies (where the dependent parameter is a floating point number or
array thereof) also become equations during initialization, and the dependent parameters
become unknowns.

`remake` will reconstruct the initialization system and problem, given the new
constraints provided to it. The new values will be combined with the original
variable-value mapping provided to `ODEProblem` and used to construct the initialization
problem.

### Parameter initialization by example

Consider the following system, where the sum of two unknowns is a constant parameter
`total`.

```@example paraminit
using ModelingToolkit, OrdinaryDiffEq # hidden
using ModelingToolkit: t_nounits as t, D_nounits as D # hidden

@variables x(t) y(t)
@parameters total
@mtkbuild sys = ODESystem([D(x) ~ -x, total ~ x + y], t;
    defaults = [total => missing], guesses = [total => 1.0])
```

Given any two of `x`, `y` and `total` we can determine the remaining variable.

```@example paraminit
prob = ODEProblem(sys, [x => 1.0, y => 2.0], (0.0, 1.0))
integ = init(prob, Tsit5())
@assert integ.ps[total] ≈ 3.0 # hide
integ.ps[total]
```

Suppose we want to re-create this problem, but now solve for `x` given `total` and `y`:

```@example paraminit
prob2 = remake(prob; u0 = [y => 1.0], p = [total => 4.0])
initsys = prob2.f.initializeprob.f.sys
```

The system is now overdetermined. In fact:

```@example paraminit
[equations(initsys); observed(initsys)]
```

The system can never be satisfied and will always lead to an `InitialFailure`. This is
due to the aforementioned behavior of retaining the original variable-value mapping
provided to `ODEProblem`. To fix this, we pass `x => nothing` to `remake` to remove its
retained value.

```@example paraminit
prob2 = remake(prob; u0 = [y => 1.0, x => nothing], p = [total => 4.0])
initsys = prob2.f.initializeprob.f.sys
```

The system is fully determined, and the equations are solvable.

```@example paraminit
[equations(initsys); observed(initsys)]
```

## Diving Deeper: Constructing the Initialization System

To get a better sense of the initialization system and to help debug it, you can construct
the initialization system directly. The initialization system is a NonlinearSystem
which requires the system-level information and the additional nonlinear equations being
tagged to the system.

```@example init
isys = generate_initializesystem(pend, u0map = [x => 1.0, y => 0.0], guesses = [λ => 1])
```

We can inspect what its equations and unknown values are:

```@example init
equations(isys)
```

```@example init
unknowns(isys)
```

Notice that all initial conditions are treated as initial equations. Additionally, for systems
with observables, those observables are too treated as initial equations. We can see the
resulting simplified system via the command:

```@example init
isys = structural_simplify(isys; fully_determined = false)
```

Note `fully_determined=false` allows for the simplification to occur when the number of equations
does not match the number of unknowns, which we can use to investigate our overdetermined system:

```@example init
isys = ModelingToolkit.generate_initializesystem(
    pend, u0map = [x => 1, y => 0.0, D(y) => 2.0, λ => 1], guesses = [λ => 1])
```

```@example init
isys = structural_simplify(isys; fully_determined = false)
```

```@example init
equations(isys)
```

```@example init
unknowns(isys)
```

```@example init
observed(isys)
```

After simplification we see that we have 5 equatinos to solve with 3 variables, and the
system that is given is not solvable.

## Numerical Isolation: InitializationProblem

To inspect the numerics of the initialization problem, we can use the `InitializationProblem`
constructor which acts just like an `ODEProblem` or `NonlinearProblem` constructor, but
creates the special initialization system for a given `sys`. This is done as follows:

```@example init
iprob = ModelingToolkit.InitializationProblem(pend, 0.0,
    [x => 1, y => 0.0, D(y) => 2.0, λ => 1], [g => 1], guesses = [λ => 1])
```

We can see that because the system is overdetermined we receive a NonlinearLeastSquaresProblem,
solvable by [NonlinearSolve.jl](https://docs.sciml.ai/NonlinearSolve/stable/). Using NonlinearSolve
we can recreate the initialization solve directly:

```@example init
using NonlinearSolve
sol = solve(iprob)
```

!!! note
    
    For more information on solving NonlinearProblems and NonlinearLeastSquaresProblems,
    check out the [NonlinearSolve.jl tutorials!](https://docs.sciml.ai/NonlinearSolve/stable/tutorials/getting_started/).

We can see that the default solver stalls

```@example init
sol.stats
```

after doing many iterations, showing that it tried to compute but could not find a valid solution.
Trying other solvers:

```@example init
sol = solve(iprob, GaussNewton())
```

gives the same issue, indicating that the chosen initialization system is unsatisfiable. We can
check the residuals:

```@example init
sol.resid
```

to see the problem is not equation 2 but other equations in the system. Meanwhile, changing
some of the conditions:

```@example init
iprob = ModelingToolkit.InitializationProblem(pend, 0.0,
    [x => 1, y => 0.0, D(y) => 0.0, λ => 0], [g => 1], guesses = [λ => 1])
```

gives a NonlinearLeastSquaresProblem which can be solved:

```@example init
sol = solve(iprob)
```

```@example init
sol.resid
```

In comparison, if we have a well-conditioned system:

```@example init
iprob = ModelingToolkit.InitializationProblem(pend, 0.0,
    [x => 1, y => 0.0], [g => 1], guesses = [λ => 1])
```

notice that we instead obtained a NonlinearSystem. In this case we have to use
different solvers which can take advantage of the fact that the Jacobian is square.

```@example init
sol = solve(iprob)
```

```@example init
sol = solve(iprob, TrustRegion())
```

## More Features of the Initialization System: Steady-State and Observable Initialization

Let's take a Lotka-Volterra system:

```@example init
@variables x(t) y(t) z(t)
@parameters α=1.5 β=1.0 γ=3.0 δ=1.0

eqs = [D(x) ~ α * x - β * x * y
       D(y) ~ -γ * y + δ * x * y
       z ~ x + y]

@named sys = ODESystem(eqs, t)
simpsys = structural_simplify(sys)
tspan = (0.0, 10.0)
```

Using the derivative initializations, we can set the ODE to start at the steady state
by initializing the derivatives to zero:

```@example init
prob = ODEProblem(simpsys, [D(x) => 0.0, D(y) => 0.0], tspan, guesses = [x => 1, y => 1])
sol = solve(prob, Tsit5(), abstol = 1e-16)
```

Notice that this is a "numerical zero", not an exact zero, and thus the solution will leave the
steady state in this instance because it's an unstable steady state.

Additionally, notice that in this setup we have an observable `z ~ x + y`. If we instead know the
initial condition for the observable we can use that directly:

```@example init
prob = ODEProblem(simpsys, [D(x) => 0.0, z => 2.0], tspan, guesses = [x => 1, y => 1])
sol = solve(prob, Tsit5())
```

We can check that indeed the solution does satisfy that D(x) = 0 at the start:

```@example init
sol[α * x - β * x * y]
```

```@example init
plot(sol)
```
