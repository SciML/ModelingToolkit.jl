include("shared/mtktestset.jl")

@mtktestset("Guess Propagation", "guess_propagation.jl")
@safetestset "Hierarchical Initialization Equations" include("hierarchical_initialization_eqs.jl")
@mtktestset("InitializationSystem Test", "initializationsystem.jl")
@mtktestset("Initial Values Test", "initial_values.jl")
