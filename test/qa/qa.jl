using ModelingToolkit
using Aqua
using JET
using SciMLTesting
using Test

# ExplicitImports only sees an extension module once its trigger package is loaded, so
# every weakdep is loaded here to bring the extensions into the checked module set.
using FMIImport
using OrdinaryDiffEqBDF
using OrdinaryDiffEqDefault
using OrdinaryDiffEqRosenbrock

const MTK_EXTENSIONS = (
    :MTKFMIExt,
    :MTKOrdinaryDiffEqBDFExt,
    :MTKOrdinaryDiffEqDefaultExt,
    :MTKOrdinaryDiffEqRosenbrockExt,
)

# ExplicitImports silently skips an extension that fails to load, so assert the
# extension modules actually exist rather than trusting a green run_qa.
@testset "Extensions loaded" begin
    for ext in MTK_EXTENSIONS
        @test Base.get_extension(ModelingToolkit, ext) !== nothing
    end
end

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

# Names the extensions must reach for which no public spelling exists.
const NONPUBLIC_QUALIFIED_ACCESSES = (
    # ModelingToolkitBase: the FMI extension's variable metadata helpers.
    :default_toterm,
    :setinput,
    :setoutput,
    # ModelingToolkit's own stub that MTKFMIExt implements; it is documented API but has
    # not been declared `public`, and an extension has no way to add a method unqualified.
    :FMIComponent,
    # Base: `Base.RefValue` has no public spelling.
    :RefValue,
)

run_qa(
    ModelingToolkit;
    Aqua = Aqua,
    JET = JET,
    jet = true,
    jet_kwargs = (; target_defined_modules = true),
    ei_kwargs = (;
        all_qualified_accesses_are_public = (; ignore = NONPUBLIC_QUALIFIED_ACCESSES),
    ),
    api_docs_kwargs = (; ignore = SYMBOLICS_OWNED_REEXPORTS),
)
