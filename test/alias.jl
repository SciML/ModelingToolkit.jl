using Test
using ModelingToolkit: AliasGraph

ag = AliasGraph(10)
ag[1] = 1 => 2
ag[2] = -1 => 3
ag[4] = -1 => 1
ag[5] = -1 => 4
for _ in 1:5 # check ag is robust
    @test ag[1] == (-1, 3)
    @test ag[2] == (-1, 3)
    @test ag[4] == (1, 3)
    @test ag[5] == (-1, 3)
end

@test 1 in keys(ag)
@test 2 in keys(ag)
@test !(3 in keys(ag))
@test 4 in keys(ag)
@test 5 in keys(ag)

# Import breaking example from DAECompiler when applied to MSL's rotational testcase
using ModelingToolkit
using ModelingToolkit: AbstractSystem, SparseMatrixCLIL, StructuralTransformations
using ModelingToolkit.SystemStructures: TransformationState, DiffGraph, SystemStructure
using ModelingToolkit.BipartiteGraphs: BipartiteGraph

# Create a dummy system and state so that we can override `linear_subsys_adjmat!`
# for our test case to avoid having to perform a bunch of difficult computation.
struct DummySystem <: AbstractSystem
end
struct DummyState <: TransformationState{DummySystem}
    structure::Any
end

# Embed the test case data
fadjlist = [
    [1], [3, 6], [3, 8], [4, 31], [5, 32], [5, 7, 9],
    [10, 13], [10, 15], [11, 33], [12, 34], [12, 14, 16],
    [17, 19, 21], [18, 22], [18, 20], [17, 18], [23, 27, 29],
    [24, 35], [25, 36], [26, 30], [26, 28], [24, 26],
    [1, 27], [1, 8], [1, 19], [2, 9, 20, 28], [21, 29],
    [13, 21], [14, 22, 30], [7], [16],
]
badjlist = [
    [1, 22, 23, 24], [25], [2, 3], [4], [5, 6], [2], [6, 29],
    [3, 23], [6, 25], [7, 8], [9], [10, 11], [7, 27],
    [11, 28], [8], [11, 30], [12, 15], [13, 14, 15],
    [12, 24], [14, 25], [12, 26, 27], [13, 28], [16],
    [17, 21], [18], [19, 20, 21], [16, 22], [20, 25],
    [16, 26], [19, 28], [4], [5], [9], [10], [17], [18],
]

# Some variables were differentiated
var_to_diff = DiffGraph(36, false)
for (val, idx) in enumerate([3, 4, 10, 11, 23, 24])
    var_to_diff[idx] = 30 + val
end
eq_to_diff = DiffGraph(30, false)

function ModelingToolkit.linear_subsys_adjmat!(::DummyState)
    # Just return known value instead of calculating it
    return ModelingToolkit.SparseMatrixCLIL(30,
                                            36,
                                            collect(1:30),
                                            deepcopy(fadjlist),
                                            [[-1], [-1, 1], [-1, 1], [1, -1], [1, -1],
                                                [-2, 1, 1],
                                                [-1, 1], [-1, 1], [1, -1], [1, -1],
                                                [-2, 1, 1],
                                                [-1, -1, 1], [1, -1], [-1, -1], [10000, -1],
                                                [-1, -1, 1], [1, -1], [1, -1], [1, -1],
                                                [-1, -1],
                                                [10, -1], [-1, 1], [-1, 1], [-1, 1],
                                                [1, 1, 1, 1],
                                                [-1, 1], [1, -1], [1, 1, 1], [1], [1]])
end

function StructuralTransformations.var_derivative!(state::DummyState, var)
    return StructuralTransformations.var_derivative_graph!(state.structure, var)
end

graph = BipartiteGraph(deepcopy(fadjlist), deepcopy(badjlist))
solvable_graph = BipartiteGraph(deepcopy(fadjlist), deepcopy(badjlist))
structure = SystemStructure(complete(var_to_diff),
                            complete(eq_to_diff),
                            graph,
                            solvable_graph,
                            nothing,
                            false)

# Ensure that both `mm` results from `alias_eliminate_graph!()` do not contain
# single-entry rows, as this represents a trivially-zero coefficient in the
# homogenous equation that `mm` represents, and should have been reduced to zero.
_, mm1, ag, mm2 = ModelingToolkit.alias_eliminate_graph!(DummyState(structure))
@test all(count(!iszero, mm1, dims = 2) .!= 1)
@test all(count(!iszero, mm2, dims = 2) .!= 1)
