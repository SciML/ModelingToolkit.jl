using ModelingToolkit
using ModelingToolkitBase
using ModelingToolkitTearing
using StateSelection
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

# ModelingToolkit is the upper half of ModelingToolkitBase: the two are one library split
# across a monorepo boundary, with the compilation and problem-construction layer living
# here and dispatching on the system and operator types defined below it. Methods added to
# those types are extensions of the pair's own API, so Aqua should not count them as piracy.
const MTKBASE_OWNED_TYPES = (
    ModelingToolkitBase.Hold,
    ModelingToolkitBase.Sample,
    ModelingToolkitBase.SampleTime,
    ModelingToolkitBase.ShiftIndex,
    ModelingToolkitBase.System,
)

# ModelingToolkit is also where the structural-transformation stack is wired to
# ModelingToolkitBase's hooks, which means forwarding a handful of MTKBase functions to the
# StateSelection/ModelingToolkitTearing implementations on their own types. `TearingState`
# is not yet public in ModelingToolkitTearing, so it has to be reached by qualified access.
const STRUCTURAL_TYPES = (
    ModelingToolkitTearing.TearingState,
    StateSelection.DiffGraph,
)

run_qa(
    ModelingToolkit;
    Aqua = Aqua,
    JET = JET,
    jet = true,
    aqua_kwargs = (;
        piracies = (; treat_as_own = (MTKBASE_OWNED_TYPES..., STRUCTURAL_TYPES...)),
    ),
    ei_kwargs = (;
        all_qualified_accesses_are_public = (; ignore = NONPUBLIC_QUALIFIED_ACCESSES),
    ),
)
