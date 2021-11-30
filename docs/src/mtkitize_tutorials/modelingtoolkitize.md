# Automatically Accelerating ODEProblem Code

For some `DEProblem` types, automatic tracing functionality is already included
via the `modelingtoolkitize` function. Take, for example, the Robertson ODE
defined as an `ODEProblem` for DifferentialEquations.jl:

```julia
using DifferentialEquations
function rober(du,u,p,t)
  y₁,y₂,y₃ = u
  k₁,k₂,k₃ = p
  du[1] = -k₁*y₁+k₃*y₂*y₃
  du[2] =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  du[3] =  k₂*y₂^2
  nothing
end
prob = ODEProblem(rober,[1.0,0.0,0.0],(0.0,1e5),(0.04,3e7,1e4))
```

If we want to get a symbolic representation, we can simply call `modelingtoolkitize`
on the `prob`, which will return an `ODESystem`:

```julia
sys = modelingtoolkitize(prob)
```

Using this, we can symbolically build the Jacobian and then rebuild the ODEProblem:

```julia
prob_jac = ODEProblem(sys,[],(0.0,1e5),jac=true)
```
