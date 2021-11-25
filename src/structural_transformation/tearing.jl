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

    rg = ResidualCMOGraph(graph, var_eq_matching)

    # Update bipartite graph
    var_deps = Union{Vector{Int}, Nothing}[nothing for v in eliminated_variables]
    var_idxs = Dict(v => i for (i,v) in enumerate(eliminated_variables))
    for (i, v) in enumerate(eliminated_variables)
        isa(var_deps[i], Vector{Int}) && continue
        var_deps[i] = deps = Vector{Int}()
        for v′ in neighborhood(rg, v)
            if var_rename[v′] != 0
                push!(deps, var_rename[v′])
            else
                var_deps[var_idxs[v′]] = deps
            end
        end
    end

    new_fadjlist = Vector{Int}[
        let new_list = Vector{Int}()
            for v in graph.fadjlist[i]
                if var_rename[v] != 0
                    push!(new_list, var_rename[v])
                else
                    append!(new_list, var_deps[var_idxs[v]])
                end
            end
            new_list
        end for i = 1:nsrcs(graph) if eq_rename[i] != 0]

    return BipartiteGraph(new_fadjlist, ndsts(graph) - length(eliminated_variables))
end

"""
    algebraic_variables_scc(sys)

Find strongly connected components of algebraic variables in a system.
"""
function algebraic_variables_scc(sys)
    s = structure(sys)
    if !(s isa SystemStructure)
        sys = initialize_system_structure(sys)
        s = structure(sys)
    end

    # skip over differential equations
    algvars = isalgvar.(Ref(s), 1:ndsts(s.graph))

    var_eq_matching = complete(maximal_matching(s, e->s.algeqs[e], v->algvars[v]))
    var_sccs = find_var_sccs(complete(s.graph), var_eq_matching)

    return var_eq_matching, var_sccs
end
