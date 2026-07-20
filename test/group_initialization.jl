include("shared/mtktestset.jl")

@mtktestset("Guess Propagation", "guess_propagation.jl")
@safetestset "Hierarchical Initialization Equations" include("hierarchical_initialization_eqs.jl")
@safetestset "Initialization Jacobian analysis" include("initialization_jacobian_analysis.jl")
@mtktestset("InitializationSystem Test", "initializationsystem.jl")
@mtktestset("Initial Values Test", "initial_values.jl")
@safetestset "Homotopy initialization (SCC bypass)" include("homotopy_initialization_scc.jl")
