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
@safetestset "Selection at the initial point" begin
    include("initial_point_selection.jl")
end
