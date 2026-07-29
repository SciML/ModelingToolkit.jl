using ModelingToolkit
using Aqua
using JET
using SciMLTesting

# Public names that reach ModelingToolkit's API surface only through
# `@reexport using Symbolics` in ModelingToolkitBase. They are owned (and undocumented) by
# Symbolics/SymbolicUtils, so ModelingToolkit is not the right place to document them.
const SYMBOLICS_OWNED_REEXPORTS = (
    Symbol("@symbolic_wrap"),
    Symbol("@wrapped"),
    :RuleSet,
    :get_canonical_expr,
    :infimum,
    :is_derivative,
    :istree,
    :solve_for,
    :supremum,
)

run_qa(
    ModelingToolkit;
    Aqua = Aqua,
    JET = JET,
    jet = true,
    jet_kwargs = (; target_defined_modules = true),
    api_docs_kwargs = (; ignore = SYMBOLICS_OWNED_REEXPORTS),
)
