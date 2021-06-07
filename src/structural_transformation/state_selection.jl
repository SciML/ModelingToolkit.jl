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

    @set! s.inv_varassoc = inverse_mapping(s.varassoc)
    @set! s.inv_eqassoc = inverse_mapping(s.eqassoc)
    @set! s.inv_assign = inverse_mapping(assign)
    @set! s.assign = assign
    @set! s.scc = scc

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

    eq_constraint_set = Vector{Vector{Int}}[]
    var_constraint_set = Vector{Vector{Int}}[]

    @unpack graph = s
    unknow_subgraph = BipartiteGraph(nsrcs(graph), ndsts(graph), Val(false))
    original_vars = trues(ndsts(graph))

    for c in scc
        check_component(c, eqassoc) || break # skip lower order equations
        eq_constraints, var_constraints = sorted_constraints(c, s)
        push!(eq_constraint_set, eq_constraints)
        push!(var_constraint_set, var_constraints)

        var_set = Set(var_constraints)
        for (i, eqs) in enumerate(eq_constraints), eq in eqs
            for var in 𝑠neighbors(graph, eq); var in var_set || continue
                add_edge!(unknow_subgraph, eq, var)
            end

            if inv_eqassoc[eq] != 0 # differentiated equations
                for var in 𝑠neighbors(graph, eq)
                    original_vars[var] = false
                end
            end
        end
    end

    for constraints in Iterators.reverse(zip(eq_constraints, var_constraints))
        highest_order = length(first(constraints))
        for (order, (eq_constraint, var_constraint)) in enumerate(zip(constraints...))
            if length(eq_constraint) == 1 == length(var_constraints)
                if order < highest_order
                    # There is only one equation and one unknown in this local
                    # contraint set, while the equation is lower order. This
                    # implies that the unknown must be a dummy state.
                    v = var_constraint[1]
                    @assert !s.varmask[v]
                    s.varmask[v] = true
                else

                end
            end
        end
    end

    @set! sys.structure = s
    return sys
end

@noinline scc_throw() = error("Internal error: " *
                              "Strongly connected component invariant is violated!\n" *
                              "Please file a bug report with full stack trace " *
                              "and a reproducer.")
function check_component(c, eqassoc) # return `is_highest_order`
    # Fast paths
    isempty(c) && return true
    length(c) == 1 && return eqassoc[c[1]] == 0

    # If there are more than one equation in the component, it must contain
    # highest order equations only.
    for eq in c; eqassoc[eq] == 0 && continue
        scc_throw()
    end
    return true
end

"""
$(TYPEDSIGNATURES)

The constraints are sorted in descending order.
"""
function sorted_constraints(c::Vector{Int}, s::SystemStructure)
    @unpack inv_varassoc, inv_eqassoc, inv_assign = s
    # The set of variables that this compoent solves for
    unknowns = [inv_assign[eq] for eq in c]
    eq_constraints = [c]
    var_constraints = [unknowns]
    while true
        lower_eqs = Int[]
        for eq in eq_constraints[end]
            lower_eq = inv_eqassoc[eq]
            lower_eq > 0 && push!(lower_eqs, lower_eq)
        end
        isempty(lower_eqs) && break

        push!(eq_constraints, lower_eqs)

        vars = Int[]
        for var in var_constraints[end]
            lower_var = inv_varassoc[var]
            lower_var && push!(vars, lower_var)
        end
        push!(var_constraints, vars)

        length(vars) < length(lower_eqs) && scc_throw()
    end
    return eq_constraints, var_constraints
end
