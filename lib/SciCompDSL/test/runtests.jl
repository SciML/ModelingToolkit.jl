using SafeTestsets, Pkg, Test

# The centralized SciML sublibrary CI (SublibraryCI.yml -> sublibrary-tests.yml@v1)
# emits GROUP="SciCompDSL" for the [Core] section and
# GROUP="SciCompDSL_<Section>" for every other section in test_groups.toml.
# Strip that "<pkg>_" prefix so the bare section names below ("Core", "QA", …)
# drive the dispatch. A bare "<pkg>" maps to "Core"; anything else (e.g. "All" for
# local runs) is passed through unchanged.
const _G = get(ENV, "GROUP", "All")
const _SUB = "SciCompDSL"
const GROUP = _G == _SUB ? "Core" :
    (startswith(_G, _SUB * "_") ? _G[(length(_SUB) + 2):end] : _G)

const MTKBasePath = joinpath(dirname(dirname(@__DIR__)), "ModelingToolkitBase")
const MTKBasePkgSpec = PackageSpec(; path = MTKBasePath)

const MTKPath = dirname(dirname(dirname(@__DIR__)))
const MTKPkgSpec = PackageSpec(; path = MTKPath)

function activate_qa_env()
    Pkg.activate("qa")
    Pkg.develop([PackageSpec(path = dirname(@__DIR__))])
    return Pkg.instantiate()
end

@time begin
    if GROUP == "All" || GROUP == "Core"
        Pkg.develop([MTKBasePkgSpec, MTKPkgSpec]; preserve = Pkg.PRESERVE_ALL)

        @safetestset "Model macro hygiene" include("macro_hygiene.jl")
        @safetestset "Model parsing - MTKBase" include("model_parsing.jl")
        @safetestset "Model parsing - MTK" begin
            using ModelingToolkit
            import ModelingToolkitBase
            include("model_parsing.jl")
        end
    end

    if GROUP == "All" || GROUP == "QA"
        activate_qa_env()
        @safetestset "Aqua Tests" include("qa/aqua.jl")
    end
end
