function try_assign_eq!(ict::IncrementalCycleTracker, vj::Integer, eq::Integer)
    G = ict.graph
    add_edge_checked!(ict, Iterators.filter(!=(vj), 𝑠neighbors(G.graph, eq)), vj) do G
        G.matching[vj] = eq
        G.ne += length(𝑠neighbors(G.graph, eq)) - 1
    end
end
function try_assign_eq!(ict::IncrementalCycleTracker, vars, v_active, eq::Integer,
        condition::F = _ -> true) where {F}
    G = ict.graph
    for vj in vars
        (vj in v_active && G.matching[vj] === unassigned && condition(vj)) || continue
        try_assign_eq!(ict, vj, eq) && return true
    end
    return false
end
function tearEquations!(ict::IncrementalCycleTracker, Gsolvable, es::Vector{Int},
        v_active::BitSet, isder′::F) where {F}
    check_der = isder′ !== nothing
    if check_der
        has_der = Ref(false)
        isder = let has_der = has_der, isder′ = isder′
            v -> begin
                r = isder′(v)
                has_der[] |= r
                r
            end
        end
    end
    for only_single_solvable in (true, false)
        for eq in es
            vs = Gsolvable[eq]
            ((length(vs) == 1) ⊻ only_single_solvable) && continue
            if check_der
                try_assign_eq!(ict, vs, v_active, eq, isder)
                if has_der[]
                    has_der[] = false
                    continue
                end
            end
            try_assign_eq!(ict, vs, v_active, eq)
        end
    end
    return ict
end
function tear_graph_block_modia!(var_eq_matching, ict, solvable_graph, eqs, vars,
        isder::F) where {F}
    tearEquations!(ict, solvable_graph.fadjlist, eqs, vars, isder)
    for var in vars
        var_eq_matching[var] = ict.graph.matching[var]
    end
    return nothing
end
function build_var_eq_matching(structure::SystemStructure, ::Type{U} = Unassigned;
        varfilter::F2 = v -> true, eqfilter::F3 = eq -> true) where {U, F2, F3}
    @unpack graph, solvable_graph = structure
    var_eq_matching = maximal_matching(graph, eqfilter, varfilter, U)
    matching_len = max(length(var_eq_matching),
        maximum(x -> x isa Int ? x : 0, var_eq_matching, init = 0))
    return complete(var_eq_matching, matching_len), matching_len
end
function tear_graph_modia(structure::SystemStructure, isder::F = nothing,
        ::Type{U} = Unassigned;
        varfilter::F2 = v -> true,
        eqfilter::F3 = eq -> true) where {F, U, F2, F3}
    @unpack graph, solvable_graph = structure
    var_eq_matching, matching_len = build_var_eq_matching(structure, U; varfilter, eqfilter)
    full_var_eq_matching = copy(var_eq_matching)
    var_sccs = find_var_sccs(graph, var_eq_matching)
    vargraph = DiCMOBiGraph{true}(graph, 0, Matching(matching_len))
    ict = IncrementalCycleTracker(vargraph; dir = :in)
    ieqs = Int[]
    filtered_vars = BitSet()
    free_eqs = free_equations(graph, var_sccs, var_eq_matching, varfilter)
    is_overdetemined = !isempty(free_eqs)
    for vars in var_sccs
        for var in vars
            if varfilter(var)
                push!(filtered_vars, var)
                if var_eq_matching[var] !== unassigned
                    ieq = var_eq_matching[var]
                    push!(ieqs, ieq)
                end
            end
            var_eq_matching[var] = unassigned
        end
        tear_graph_block_modia!(var_eq_matching, ict, solvable_graph, ieqs,
            filtered_vars, isder)
        if !is_overdetemined
            vargraph.ne = 0
            for var in vars
                vargraph.matching[var] = unassigned
            end
        end
        empty!(ieqs)
        empty!(filtered_vars)
    end
    if is_overdetemined
        free_vars = findall(x -> !(x isa Int), var_eq_matching)
        tear_graph_block_modia!(var_eq_matching, ict, solvable_graph, free_eqs,
            BitSet(free_vars), isder)
    end
    return var_eq_matching, full_var_eq_matching, var_sccs
end
