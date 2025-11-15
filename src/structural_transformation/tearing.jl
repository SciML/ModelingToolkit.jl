struct EquationSolveError
    eq::Any
    var::Any
    rhs::Any
end

function Base.showerror(io::IO, ese::EquationSolveError)
    print(io, "EquationSolveError: While solving\n\n\t")
    print(io, ese.eq)
    print(io, "\nfor ")
    printstyled(io, var, bold = true)
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

function contract_variables(graph::BipartiteGraph, var_eq_matching::Matching,
        var_rename, eq_rename, nelim_eq, nelim_var)
    dig = DiCMOBiGraph{true}(graph, var_eq_matching)

    # Update bipartite graph
    var_deps = map(1:ndsts(graph)) do v
        [var_rename[v′]
         for v′ in neighborhood(dig, v, Inf; dir = :in) if var_rename[v′] != 0]
    end

    newgraph = BipartiteGraph(nsrcs(graph) - nelim_eq, ndsts(graph) - nelim_var)
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
    algvars = BitSet(findall(v -> isalgvar(state.structure, v), 1:ndsts(graph)))
    algeqs = BitSet(findall(map(1:nsrcs(graph)) do eq
        all(v -> !isdervar(state.structure, v),
            𝑠neighbors(graph, eq))
    end))
    var_eq_matching = complete(
        maximal_matching(graph, e -> e in algeqs, v -> v in algvars), ndsts(graph))
    var_sccs = find_var_sccs(complete(graph), var_eq_matching)

    return var_eq_matching, var_sccs
end

function free_equations(graph, vars_scc, var_eq_matching, varfilter::F) where {F}
    ne = nsrcs(graph)
    seen_eqs = falses(ne)
    for vars in vars_scc, var in vars

        varfilter(var) || continue
        ieq = var_eq_matching[var]
        if ieq isa Int
            seen_eqs[ieq] = true
        end
    end
    findall(!, seen_eqs)
end

struct SelectedState end
const MatchingT{T} = Matching{T, Vector{Union{T, Int}}}
const MatchedVarT = Union{Unassigned, SelectedState}
const VarEqMatchingT = MatchingT{MatchedVarT}

"""
    $TYPEDEF

A struct containing the results of tearing.

# Fields

$TYPEDFIELDS
"""
struct TearingResult
    """
    The variable-equation matching. Differential variables are matched to `SelectedState`.
    The derivative of a differential variable is matched to the corresponding differential
    equation. Solved variables are matched to the equation they are solved from. Algebraic
    variables are matched to `unassigned`.
    """
    var_eq_matching::VarEqMatchingT
    """
    The variable-equation matching prior to tearing. This is the maximal matching used to
    compute `var_sccs` (see below). For generating the torn system, `var_eq_matching` is
    the source of truth. This should only be used to identify algebraic equations in each
    SCC.
    """
    full_var_eq_matching::VarEqMatchingT
    """
    The partitioning of variables into strongly connected components (SCCs). The SCCs are
    sorted in dependency order, so each SCC depends on variables in previous SCCs.
    """
    var_sccs::Vector{Vector{Int}}
end

"""
    $TYPEDEF

Supertype for all tearing algorithms. A tearing algorithm takes as input the
`SystemStructure` along with any other necessary arguments.

The output of a tearing algorithm must be a `TearingResult` and a `NamedTuple` of
any additional data computed in the process that may be useful for further processing.
"""
abstract type TearingAlgorithm end

"""
    $TYPEDEF

Supertype for all reassembling algorithms. A reassembling algorithm takes as input the
`TearingState`, `TearingResult` and integer incidence matrix `mm::SparseMatrixCLIL`. The
matrix `mm` may be `nothing`. The algorithm must also accept arbitrary keyword arguments.
The following keyword arguments will always be provided:
- `fully_determined::Bool`: flag indicating whether the system is fully determined.

The output of a reassembling algorithm must be the torn system.

A reassemble algorithm must also implement `with_fully_determined`
"""
abstract type ReassembleAlgorithm end
