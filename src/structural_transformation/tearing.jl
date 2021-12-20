struct EquationSolveError
    eq
    var
    rhs
end

function Base.showerror(io::IO, ese::EquationSolveError)
    print(io, "EquationSolveError: While solving\n\n\t")
    print(io, ese.eq)
    print(io, "\nfor ")
    printstyled(io, var, bold=true)
    print(io, ", obtained RHS\n\n\tt")
    println(io, rhs)
end

function masked_cumsum!(A::Vector)
    acc = zero(eltype(A))
    for i in eachindex(A)
        iszero(A[i]) && continue
        A[i] = (acc += A[i])
    end
end

function contract_variables(graph::BipartiteGraph, var_eq_matching::Matching, eliminated_variables)
    var_rename = ones(Int64, ndsts(graph))
    eq_rename = ones(Int64, nsrcs(graph))
    for v in eliminated_variables
        eq_rename[var_eq_matching[v]] = 0
        var_rename[v] = 0
    end
    masked_cumsum!(var_rename)
    masked_cumsum!(eq_rename)

    dig = DiCMOBiGraph{true}(graph, var_eq_matching)

    # Update bipartite graph
    var_deps = map(1:ndsts(graph)) do v
        [var_rename[v′] for v′ in neighborhood(dig, v, Inf; dir=:in) if var_rename[v′] != 0]
    end

    nelim = length(eliminated_variables)
    newgraph = BipartiteGraph(nsrcs(graph) - nelim, ndsts(graph) - nelim)
    for e in 𝑠vertices(graph)
        ne = eq_rename[e]
        ne == 0 && continue
        for v in 𝑠neighbors(graph, e)
            newvar = var_rename[v]
            if newvar != 0
                add_edge!(newgraph, ne, newvar)
            else
                for nv in var_deps[v]
                    add_edge!(newgraph, ne, nv)
                end
            end
        end
    end

    return newgraph
end

"""
    algebraic_variables_scc(sys)

Find strongly connected components of algebraic variables in a system.
"""
function algebraic_variables_scc(state::TearingState)
    graph = state.structure.graph
    # skip over differential equations
    algvars = BitSet(findall(v->isalgvar(state.structure, v), 1:ndsts(graph)))
    algeqs = BitSet(findall(map(1:nsrcs(graph)) do eq
        all(v->!isdervar(state.structure, v), 𝑠neighbors(graph, eq))
    end))
    var_eq_matching = complete(maximal_matching(graph, e->e in algeqs, v->v in algvars))
    var_sccs = find_var_sccs(complete(graph), var_eq_matching)

    return var_eq_matching, var_sccs
end
