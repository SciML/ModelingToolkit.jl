# Modelingtoolkitize: Automatically Translating Numerical to Symbolic Code

## What is `modelingtoolkitize`?

From the other tutorials you will have learned that ModelingToolkit is a symbolic library
with all kinds of goodies, such as the ability to derive analytical expressions for things
like Jacobians, determine the sparsity of a set of equations, perform index reduction,
tearing, and other transformations to improve both stability and performance. All of these
are good things, but all of these require that one has defined the problem symbolically.

**But what happens if one wants to use ModelingToolkit functionality on code that is already
written for DifferentialEquations.jl, NonlinearSolve.jl, Optimization.jl, or beyond?**

`modelingtoolktize` is a function in ModelingToolkit which takes a numerically-defined
`SciMLProblem` and transforms it into its symbolic ModelingToolkit equivalent. By doing
so, ModelingToolkit analysis passes and transformations can be run as intermediate steps
to improve a simulation code before it's passed to the solver.

!!! note
    
    `modelingtoolkitize` does have some limitations, i.e. not all codes that work with the
    numerical solvers will work with `modelingtoolkitize`. Namely, it requires the ability
    to trace the equations with Symbolics.jl `Num` types. Generally, a code which is
    compatible with forward-mode automatic differentiation is compatible with
    `modelingtoolkitize`.

!!! warn
    
    `modelingtoolkitize` expressions cannot keep control flow structures (loops), and thus
    equations with long loops will be translated into large expressions, which can increase
    the compile time of the equations and reduce the SIMD vectorization achieved by LLVM.

## Example Usage: Generating an Analytical Jacobian Expression for an ODE Code

Take, for example, the Robertson ODE
defined as an `ODEProblem` for OrdinaryDiffEq.jl:

```@example mtkize
using OrdinaryDiffEq, ModelingToolkit
function rober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃
    du[3] = k₂ * y₂^2
    nothing
end
prob = ODEProblem(rober, [1.0, 0.0, 0.0], (0.0, 1e5), (0.04, 3e7, 1e4))
```

If we want to get a symbolic representation, we can simply call `modelingtoolkitize`
on the `prob`, which will return an `ODESystem`:

```@example mtkize
@mtkbuild sys = modelingtoolkitize(prob)
```

Using this, we can symbolically build the Jacobian and then rebuild the ODEProblem:

```@example mtkize
prob_jac = ODEProblem(sys, [], (0.0, 1e5), jac = true)
```
