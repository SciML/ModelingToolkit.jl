# Modeling Optimization Problems

```julia
using ModelingToolkit, GalacticOptim

@variables x y
@parameters a b
loss = (a - x)^2 + b * (y - x^2)^2
sys = OptimizationSystem(loss,[x,y],[a,b])

u0 = [
    x=>1.0
    y=>2.0
]
p = [
    a => 6.0
    b => 7.0
]

prob = OptimizationProblem(sys,u0,p,grad=true,hess=true)
solve(prob,Newton())
```

Needs more text but it's super cool and auto-parallelizes and sparsifies too.
Plus you can hierarchically nest systems to have it generate huge
optimization problems.
