using SafeTestsets

@safetestset "Utilities" begin
    include("utils.jl")
end
@safetestset "Index Reduction & SCC" begin
    include("index_reduction.jl")
end
@safetestset "Tearing" begin
    include("tearing.jl")
end
@safetestset "Array equations" begin
    include("array_equations.jl")
end
