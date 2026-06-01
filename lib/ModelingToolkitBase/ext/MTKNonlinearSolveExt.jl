module MTKNonlinearSolveExt

using ModelingToolkitBase: ModelingToolkitBase
using NonlinearSolve: NewtonRaphson

# Make `NonlinearSolve.NewtonRaphson` the default inner solver for the homotopy
# continuation algorithms (`HomotopySweep`/`TrivialHomotopy`/`TrivialThenSweep`)
# whenever NonlinearSolve is loaded. NonlinearSolve is not a runtime dependency
# of ModelingToolkitBase, so `_default_inner` cannot reference it directly; this
# extension populates the factory Ref at load time. Without NonlinearSolve loaded
# ModelingToolkitBase falls back to `SimpleNonlinearSolve.SimpleNewtonRaphson`
# (a hard `[deps]`), so a `homotopy(...)` system is still buildable out of the box.
function __init__()
    ModelingToolkitBase._DEFAULT_INNER_FACTORY[] = () -> NewtonRaphson()
    return nothing
end

end
