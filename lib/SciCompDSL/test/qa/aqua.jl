using SciCompDSL
using ModelingToolkitBase
using Aqua
using SciMLTesting
using Test

# ExplicitImports only sees an extension module once its trigger package is loaded.
using DynamicQuantities

@testset "Extensions loaded" begin
    @test Base.get_extension(SciCompDSL, :SciCompDSLDynamicQuantitiesExt) !== nothing
end

# SciCompDSL implements ModelingToolkitBase's `@connector` hook
# (`__mtkmodel_connector`). Aqua flags the method as piracy because the generic
# lives in ModelingToolkitBase; treat it as owned by this package.
const INTENTIONAL_EXTERNAL_GENERIC_EXTENSIONS = (
    ModelingToolkitBase.__mtkmodel_connector,
)

# Extension / parent seams that have no public spelling yet. Declaring them
# public would advertise internal extension contracts; ignore until the owners
# promote them.
#
#   SciCompDSL. Extension methods for DynamicQuantities unit conversion and the
#   unit-aware variable constructor. These are the package's own extension hooks,
#   not end-user API.
#
#   DynamicQuantities.SymbolicUnits. `as_quantity` is the conversion entry point
#   the extension needs; DynamicQuantities does not declare it public.
#
#   ModelingToolkitBase. `__mtkmodel_connector` is the `@connector` hook this
#   package implements; it is intentionally internal to the MTKBase/SciCompDSL pair.
#
#   Base / Core. `@nospecializeinfer` and `eval` have no public spelling.
const NONPUBLIC_EXPLICIT_IMPORTS = (
    # SciCompDSL
    :NO_VALUE, :NoValue, :convert_units,
)

const NONPUBLIC_QUALIFIED_ACCESSES = (
    # SciCompDSL
    :__generate_variable_with_unit, :convert_units,
    # DynamicQuantities.SymbolicUnits
    :as_quantity,
    # ModelingToolkitBase
    :__mtkmodel_connector,
    # Base / Core
    Symbol("@nospecializeinfer"), :eval,
)

# `@mtkmodel` is the deprecated DSL surface left alone in this PR; it is covered
# narratively in docs/src/basics/MTKLanguage.md rather than a `@docs` block.
const API_DOCS_IGNORE = (
    Symbol("@mtkmodel"),
)

run_qa(
    SciCompDSL;
    Aqua,
    aqua_kwargs = (;
        piracies = (; treat_as_own = INTENTIONAL_EXTERNAL_GENERIC_EXTENSIONS),
    ),
    ei_kwargs = (;
        all_explicit_imports_are_public = (; ignore = NONPUBLIC_EXPLICIT_IMPORTS),
        all_qualified_accesses_are_public = (; ignore = NONPUBLIC_QUALIFIED_ACCESSES),
    ),
    api_docs_kwargs = (;
        ignore = API_DOCS_IGNORE,
        rendered_ignore = API_DOCS_IGNORE,
    ),
)
