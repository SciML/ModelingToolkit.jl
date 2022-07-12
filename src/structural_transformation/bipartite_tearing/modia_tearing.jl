# This code is derived from the Modia project and is licensed as follows:
# https://github.com/ModiaSim/Modia.jl/blob/b61daad643ef7edd0c1ccce6bf462c6acfb4ad1a/LICENSE

function try_assign_eq!(ict::IncrementalCycleTracker, vj::Integer, eq::Integer)
    G = ict.graph
    add_edge_checked!(ict, Iterators.filter(!=(vj), 𝑠neighbors(G.graph, eq)), vj) do G
        G.matching[vj] = eq
        G.ne += length(𝑠neighbors(G.graph, eq)) - 1
    end
end

function tearEquations!(ict::IncrementalCycleTracker, Gsolvable, es::Vector{Int},
                        vs::Vector{Int})
    G = ict.graph
    vActive = BitSet(vs)

    for eq in es  # iterate only over equations that are not in eSolvedFixed
        for vj in Gsolvable[eq]
            if G.matching[vj] === unassigned && (vj in vActive)
                r = try_assign_eq!(ict, vj, eq)
                r && break
            end
        end
    end

    return ict
end

function tear_graph_block_modia!(var_eq_matching, graph, solvable_graph, eqs, vars)
    ict = IncrementalCycleTracker(DiCMOBiGraph{true}(graph); dir = :in)
    tearEquations!(ict, solvable_graph.fadjlist, eqs, vars)
    for var in vars
        var_eq_matching[var] = ict.graph.matching[var]
    end
    return nothing
end

function tear_graph_modia(structure::SystemStructure, ::Type{U} = Unassigned;
                          varfilter = v -> true, eqfilter = eq -> true) where {U}
    # It would be possible here to simply iterate over all variables and attempt to
    # use tearEquations! to produce a matching that greedily selects the minimal
    # number of torn variables. However, we can do this process faster if we first
    # compute the strongly connected components. In the absence of cycles and
    # non-solvability, a maximal matching on the original graph will give us an
    # optimal assignment. However, even with cycles, we can use the maximal matching
    # to give us a good starting point for a good matching and then proceed to
    # reverse edges in each scc to improve the solution. Note that it is possible
    # to have optimal solutions that cannot be found by this process. We will not
    # find them here [TODO: It would be good to have an explicit example of this.]

    @unpack graph, solvable_graph = structure
    var_eq_matching = complete(maximal_matching(graph, eqfilter, varfilter, U))
    var_sccs::Vector{Union{Vector{Int}, Int}} = find_var_sccs(graph, var_eq_matching)

    for vars in var_sccs
        filtered_vars = filter(varfilter, vars)
        ieqs = Int[var_eq_matching[v]
                   for v in filtered_vars if var_eq_matching[v] !== unassigned]
        for var in vars
            var_eq_matching[var] = unassigned
        end
        tear_graph_block_modia!(var_eq_matching, graph, solvable_graph, ieqs, filtered_vars)
    end

    return var_eq_matching
end
