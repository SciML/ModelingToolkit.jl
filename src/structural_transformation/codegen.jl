using LinearAlgebra

using ModelingToolkit: process_events

const MAX_INLINE_NLSOLVE_SIZE = 8

function torn_system_with_nlsolve_jacobian_sparsity(state, var_eq_matching, var_sccs,
        nlsolve_scc_idxs, eqs_idxs, states_idxs)
    graph = state.structure.graph

    # The sparsity pattern of `nlsolve(f, u, p)` w.r.t `p` is difficult to
    # determine in general. Consider the "simplest" case, a linear system. We
    # have
    #                   A u = p.
    # Clearly, the sparsity of `u` depends on the sparsity of both `p` and `A`
    # in a non-trivial way. However, in the generic case, `u` is dense even when
    # `A` and `p` are sparse. For instance
    #
    # ```julia
    # julia> using Random, SparseArrays
    #
    # julia> A = sprand(MersenneTwister(1234), 100, 100, 0.1);
    #
    # julia> p = sprand(MersenneTwister(12345), 100, 0.05);
    #
    # julia> count(x->abs(x) < 1e-5, A \ Vector(p))
    # 0
    # ```
    #
    # Let 𝑇 be the set of tearing variables and 𝑉 be the set of all *unknowns* in
    # the residual equations. In the following code, we are going to assume the
    # connection between 𝑇 (the `u` in from above) and 𝑉 ∖ 𝑇 (the `p` in from
    # above) has full incidence.
    #
    # Note that as we are reducing algebraic equations numerically, it could be
    # the case that a later partition (a BLT block) contains tearing variables
    # from other partitions.
    #
    # We know that partitions are BLT ordered. Hence, the tearing variables in
    # each partition is unique, and all unknowns in a partition must be
    # either differential variables or algebraic tearing variables that are
    # from previous partitions. Hence, we can build the dependency chain as we
    # traverse the partitions.

    var_rename = ones(Int64, ndsts(graph))
    nlsolve_vars = Int[]
    for i in nlsolve_scc_idxs, c in var_sccs[i]

        append!(nlsolve_vars, c)
        for v in c
            var_rename[v] = 0
        end
    end
    masked_cumsum!(var_rename)

    dig = DiCMOBiGraph{true}(graph, var_eq_matching)

    fused_var_deps = map(1:ndsts(graph)) do v
        BitSet(v′ for v′ in neighborhood(dig, v, Inf; dir = :in) if var_rename[v′] != 0)
    end

    for scc in var_sccs[nlsolve_scc_idxs]
        if length(scc) >= 2
            deps = fused_var_deps[scc[1]]
            for c in 2:length(scc)
                union!(deps, fused_var_deps[c])
                fused_var_deps[c] = deps
            end
        end
    end

    var2idx = Dict{Int, Int}(v => i for (i, v) in enumerate(states_idxs))
    eqs2idx = Dict{Int, Int}(v => i for (i, v) in enumerate(eqs_idxs))

    I = Int[]
    J = Int[]
    s = state.structure
    for ieq in 𝑠vertices(graph)
        nieq = get(eqs2idx, ieq, 0)
        nieq == 0 && continue
        for ivar in 𝑠neighbors(graph, ieq)
            isdervar(s, ivar) && continue
            if var_rename[ivar] != 0
                push!(I, nieq)
                push!(J, var2idx[ivar])
            else
                for dvar in fused_var_deps[ivar]
                    isdervar(s, dvar) && continue
                    niv = get(var2idx, dvar, 0)
                    niv == 0 && continue
                    push!(I, nieq)
                    push!(J, niv)
                end
            end
        end
    end
    sparse(I, J, true, length(eqs_idxs), length(states_idxs))
end

"""
    find_solve_sequence(sccs, vars)

given a set of `vars`, find the groups of equations we need to solve for
to obtain the solution to `vars`
"""
function find_solve_sequence(sccs, vars)
    subset = filter(i -> !isdisjoint(sccs[i], vars), 1:length(sccs))
    isempty(subset) && return []
    vars′ = mapreduce(i -> sccs[i], union, subset)
    if vars′ == vars
        return subset
    else
        return find_solve_sequence(sccs, vars′)
    end
end
