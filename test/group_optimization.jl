include("shared/mtktestset.jl")

@mtktestset("OptimizationSystem Test", "optimization/optimizationsystem.jl")
@mtktestset("InfiniteOpt Extension Test", "optimization/test_infiniteopt.jl")
@mtktestset(
    "Dynamic Optimization Collocation Solvers",
    "optimization/dynamic_optimization.jl"
)
