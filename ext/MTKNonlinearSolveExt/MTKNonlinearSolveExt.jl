module MTKNonlinearSolveExt

using ModelingToolkit
using NonlinearSolve
using SCCNonlinearSolve

# Export the TrustRegion algorithm for use in linearization
ModelingToolkit._get_default_nlsolve_alg() = TrustRegion()

end