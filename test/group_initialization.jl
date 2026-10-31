include("shared/mtktestset.jl")

@mtktestset("Guess Propagation", "guess_propagation.jl")
@safetestset "Hierarchical Initialization Equations" include("hierarchical_initialization_eqs.jl")
@safetestset "Initialization Jacobian analysis" include("initialization_jacobian_analysis.jl")
@safetestset "Residual slice initialization with algebraic unknown" include("residual_slice_initialization.jl")
@mtktestset("InitializationSystem Test", "initializationsystem.jl")
@mtktestset("Initial Values Test", "initial_values.jl")
@safetestset "strip_initialization_symbolics, nonlinear SCC" include("strip_initialization_symbolics_scc.jl")
@safetestset "Homotopy initialization (SCC bypass)" include("homotopy_initialization_scc.jl")

@safetestset "Model-independent initialization data" include("model_independent_initialization.jl")
