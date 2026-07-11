using ModelingToolkit
using Aqua
using JET
using SciMLTesting

run_qa(
    ModelingToolkit;
    Aqua = Aqua,
    JET = JET,
    jet = true,
    jet_kwargs = (; target_defined_modules = true),
    api_docs_kwargs = (; rendered = true),
)
