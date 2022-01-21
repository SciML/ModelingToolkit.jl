# This code is from the Modia project and is licensed as follows:
# https://github.com/ModiaSim/Modia.jl/blob/b61daad643ef7edd0c1ccce6bf462c6acfb4ad1a/LICENSE

################################################
#
# Functions to tear systems of equations
#
# Author: Martin Otter, DLR-SR (first version: Jan. 14, 2017)
#
# Details are described in the paper:
#   Otter, Elmqvist (2017): Transformation of Differential Algebraic Array Equations to
#                           Index One Form. Modelica'2017 Conference.
#
################################################

"""
    (eSolved, vSolved, eResidue, vTear) = tearEquations!(td, Gsolvable, es, vs)

Equations es shall be solved with respect to variables vs. The function returns
the teared equation so that if vTear is given, vSolved can be computed from eSolved
in a forward sequence (so solving eSolved[1] for vSolved[1], eSolved[2] for vSolved[2],
and so on). vTear must be selected, so that the equations eResidues are fulfilled.
Equations es are the union of eSolved and eResidue.
Variables vs are the union of vSolved and vTear.

Gsolvable defines the variables that can be explicitly solved in every equation without influencing the solution space
(= rank preserving operation).
"""
function tearEquations!(ict::IncrementalCycleTracker, Gsolvable, es::Vector{Int}, vs::Vector{Int})
    G = ict.graph
    vActive = BitSet(vs)

    for eq in es  # iterate only over equations that are not in eSolvedFixed
        for vj in Gsolvable[eq]
            if G.matching[vj] === unassigned && (vj in vActive)
                r = add_edge_checked!(ict, Iterators.filter(!=(vj), 𝑠neighbors(G.graph, eq)), vj) do G
                    G.matching[vj] = eq
                    G.ne += length(𝑠neighbors(G.graph, eq)) - 1
                end
                r && break
            end
        end
    end

    return ict
end

"""
    tear_graph_modia(sys) -> sys

Tear the bipartite graph in a system. End users are encouraged to call [`structural_simplify`](@ref)
instead, which calls this function internally.
"""
function tear_graph_modia(graph::BipartiteGraph, solvable_graph::BipartiteGraph; varfilter=v->true, eqfilter=eq->true)
    var_eq_matching = complete(maximal_matching(graph, eqfilter, varfilter))
    var_sccs::Vector{Union{Vector{Int}, Int}} = find_var_sccs(graph, var_eq_matching)

    for vars in var_sccs
        filtered_vars = filter(varfilter, vars)
        ieqs = Int[var_eq_matching[v] for v in filtered_vars if var_eq_matching[v] !== unassigned]

        ict = IncrementalCycleTracker(DiCMOBiGraph{true}(graph); dir=:in)
        tearEquations!(ict, solvable_graph.fadjlist, ieqs, filtered_vars)

        for var in vars
            var_eq_matching[var] = ict.graph.matching[var]
        end
    end

    return var_eq_matching
end
