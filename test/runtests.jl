using SafeTestsets, Pkg, Test

const GROUP = get(ENV, "GROUP", "All")

function activate_extensions_env()
    Pkg.activate("extensions")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
end

function activate_downstream_env()
    Pkg.activate("downstream")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
end

@time begin
    if GROUP == "All" || GROUP == "InterfaceI"
        @testset "InterfaceI" begin
            @safetestset "Equation Type Accessors Test" include("equation_type_accessors.jl")           
        end
    end

    if GROUP == "All" || GROUP == "InterfaceI" || GROUP == "SymbolicIndexingInterface"
        @safetestset "SymbolicIndexingInterface test" include("symbolic_indexing_interface.jl")
        @safetestset "MTKParameters Test" include("mtkparameters.jl")
    end

    if GROUP == "All" || GROUP == "InterfaceII"
        println("C compilation test requires gcc available in the path!")
        @safetestset "C Compilation Test" include("ccompile.jl")
        @testset "Distributed Test" include("distributed.jl")
        @testset "Serialization" include("serialization.jl")
    end

    if GROUP == "All" || GROUP == "RegressionI"
        @safetestset "Latexify recipes Test" include("latexify.jl")
    end

    if GROUP == "All" || GROUP == "Downstream"
        activate_downstream_env()
        @safetestset "Linearization Tests" include("downstream/linearize.jl")
        @safetestset "Inverse Models Test" include("downstream/inversemodel.jl")
    end

    if GROUP == "All" || GROUP == "Extensions"
        activate_extensions_env()
        @safetestset "BifurcationKit Extension Test" include("extensions/bifurcationkit.jl")
    end
end

@safetestset "Model Parsing Test" include("model_parsing.jl")
