function prepare_state_selection!(sys)
    s = get_structure(sys)
    if !(s isa SystemStructure)
        sys = initialize_system_structure(sys)
        s = structure(sys)
    end
    sys, assign, eqassoc = pantelides!(sys; kwargs...)
    s = get_structure(sys)

    # TODO: use eqmask
    # When writing updating a graph we often have two choices:
    #   1. Construct a new graph with updated vertices and edges;
    #   2. Make all the functions take masks for equations and variables.
    #
    # We want to use the second choice as much as possible as it minimizes
    # memory footprint while potentially strengthen temporal locality.
    #
    # N.B.: assign changes meaning here. In Pantelides, `assign` is w.r.t. to
    # all the equations and variables with the highest derivative. Here,
    # `assign` is w.r.t. equations and variables with the highest derivative.
    highest_order_graph = BipartiteGraph(0, length(assign), Val(false))
    eq_reidx = Int[]
    for (i, es) in enumerate(eqassoc); es == 0 || continue
        vars = graph.fadjlist[i]
        # N.B.: this is an alias to the original graph
        push!(highest_order_graph.fadjlist, vars)
        highest_order_graph.ne += length(vars)
        push!(eq_reidx, i)
    end
    # matching is read-only, so aliasing is fine
    assign = matching(highest_order_graph, s.varmask)

    # Compute SCC and map the indices back to the original system
    scc = Vector{Int}[]
    for component in find_scc(highest_order_graph, assign)
        push!(scc, eq_reidx[component])
    end
    for (j, a) in enumerate(assign); a == 0 && continue
        assign[j] = eq_reidx[a]
    end

    # Otter and Elmqvist (2017) 4.3.1:
    # All variables $v_{j}$ that appear differentiated and their derivatives
    # with exception of their highest derivatives are potential states, i.e.,
    # with slight abuse of notation:
    # $v_{j} ∈ u$, where $u' = f(u, p, t)$.
    #
    # The `varmask` will be used to sort the algebraic equations, so being a
    # $v_{j}$ is a state => `varmask[j]` is `false`.
    #
    # Since $a = varassoc[j]$ when $a > 0 => vⱼ = v̇ₐ$, it follows that
    # $a > 0 => vⱼ$ is a potential state.
    for (j, a) in enumerate(s.varassoc)
        is_potential_state = a > 0
        s.varmask[j] = !is_potential_state
    end

    e_constraint_set = Vector{Vector{Int}}[]
    v_constraint_set = Vector{Vector{Int}}[]

    tearing_graph = BipartiteGraph(nsrcs(s.graph), 0, Val(false))

    function check_component(c, eqassoc) # return is_highest_order
        # Fast paths
        isempty(c) && return true
        length(c) == 1 && return eqassoc[c[1]] == 0

        # If there are more than one equation in the component, it must contain
        # highest order equations only.
        for eq in c; eqassoc[eq] == 0 && continue
            error("Internal error: please file a bug report with a full reproducer.")
        end
        return true
    end

    for c in scc
        check_component(c, eqassoc) || break # skip lower order equations
        # TODO: constraint set
    end

    @set! sys.structure.assign = assign
    @set! sys.structure.inv_assign = inverse_mapping(assign)
    @set! sys.inv_varassoc = inverse_mapping(s.varassoc)
    @set! sys.inv_eqassoc = inverse_mapping(s.eqassoc)
    @set! sys.structure.scc = scc
    return sys
end
