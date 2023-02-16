using Test
using ModelingToolkit.StructuralTransformations: computed_highest_diff_variables
using ModelingToolkit: SystemStructure, AliasGraph, BipartiteGraph, DiffGraph, complete,
                       𝑑neighbors

### Test `computed_highest_diff_variables`, which is used in the Pantelides algorithm. ###

begin
    """
       Vars: x, y
       Eqs: 0 = f(x)
       Alias: ẋ = ẏ
    """
    n_vars = 4
    ag = AliasGraph(n_vars)

    # Alias: ẋ = 1 * ẏ
    ag[4] = 1 => 2

    # 0 = f(x)
    graph = complete(BipartiteGraph([Int[1]], n_vars))

    # [x, ẋ, y, ẏ]
    var_to_diff = DiffGraph([2, nothing, 4, nothing], # primal_to_diff
                            [nothing, 1, nothing, 3]) # diff_to_primal

    # [f(x)]
    eq_to_diff = DiffGraph([nothing], # primal_to_diff
                           [nothing]) # diff_to_primal
    structure = SystemStructure(var_to_diff, eq_to_diff, graph, nothing, nothing, false)
    varwhitelist = computed_highest_diff_variables(structure, ag)

    # Correct answer is: ẋ
    @assert varwhitelist == Bool[0, 1, 0, 0]
end

begin
    """
       Vars: x, y
       Eqs: 0 = f(x)
       Alias: ẋ = ẏ, ̈x = ̈y
    """
    n_vars = 6
    ag = AliasGraph(n_vars)

    # Alias: ẋ = 1 * ẏ
    ag[5] = 1 => 2
    # Alias: ẍ = 1 * ̈y
    ag[6] = 1 => 3

    # 0 = f(x)
    graph = complete(BipartiteGraph([Int[1]], n_vars))

    # [x, ẋ, ̈x, y, ẏ, ̈x]
    var_to_diff = DiffGraph([2, 3, nothing, 5, 6, nothing], # primal_to_diff
                            [nothing, 1, 2, nothing, 4, 5]) # diff_to_primal

    # [f(x)]
    eq_to_diff = DiffGraph([nothing], # primal_to_diff
                           [nothing]) # diff_to_primal

    structure = SystemStructure(var_to_diff, eq_to_diff, graph, nothing, nothing, false)
    varwhitelist = computed_highest_diff_variables(structure, ag)

    # Correct answer is: ẋ
    @assert varwhitelist == Bool[0, 1, 0, 0, 0, 0]
end
