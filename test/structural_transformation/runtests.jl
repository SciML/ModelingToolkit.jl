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
# Array-equation groups are provided by ModelingToolkitTearing; the tests run
# against releases that define them.
import ModelingToolkitTearing
if isdefined(ModelingToolkitTearing, :row_to_equation_indices)
    @safetestset "Array equations" begin
        include("array_equations.jl")
    end
end
@safetestset "Stationary linear tearing" begin
    include("stationary_linear_tearing.jl")
end
