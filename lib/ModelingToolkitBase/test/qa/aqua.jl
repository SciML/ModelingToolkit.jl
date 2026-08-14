using ModelingToolkitBase
using Aqua
using JET
using SciMLTesting

# This is loaded only by MTKBifurcationKitExt, while remaining a hard dependency so
# loading BifurcationKit alone still activates the extension.
const STALE_EXTENSION_DEPENDENCIES = [:SimpleNonlinearSolve]

# ModelingToolkitBase supplies Symbolics expression-protocol methods for equations and
# JumpProcesses jump types. These methods intentionally extend the corresponding
# SymbolicUtils generic functions, while all other piracies remain checked.
const INTENTIONAL_EXTERNAL_GENERIC_EXTENSIONS = (
    ModelingToolkitBase.SymbolicUtils.search_variables!,
    ModelingToolkitBase.toexpr,
)

run_qa(
    ModelingToolkitBase;
    Aqua,
    JET,
    jet = true,
    aqua_kwargs = (
        ;
        stale_deps = (; ignore = STALE_EXTENSION_DEPENDENCIES),
        piracies = (; treat_as_own = INTENTIONAL_EXTERNAL_GENERIC_EXTENSIONS),
    ),
)
