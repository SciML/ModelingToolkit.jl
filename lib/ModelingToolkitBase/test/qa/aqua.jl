using ModelingToolkitBase
using Aqua

# This is loaded only by MTKBifurcationKitExt, while remaining a hard dependency so
# loading BifurcationKit alone still activates the extension.
Aqua.test_all(ModelingToolkitBase; stale_deps = (; ignore = [:SimpleNonlinearSolve]))
