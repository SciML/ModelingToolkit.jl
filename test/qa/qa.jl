using ModelingToolkit
using Aqua
using JET
using SciMLTesting
import ModelingToolkitBase as MTKBase

const PIRACY_TREAT_AS_OWN = Union{Function, Type}[
    MTKBase.System,
    MTKBase.Sample,
    MTKBase.Hold,
    MTKBase.ShiftIndex,
    ModelingToolkit.StateSelection.DiffGraph,
    ModelingToolkit.MTKTearing.TearingState,
]

run_qa(
    ModelingToolkit;
    Aqua = Aqua,
    JET = JET,
    jet = true,
    aqua_kwargs = (; piracies = (; treat_as_own = PIRACY_TREAT_AS_OWN)),
    jet_kwargs = (; target_modules = (ModelingToolkit,), mode = :typo),
    api_docs_kwargs = (; rendered = true),
)
