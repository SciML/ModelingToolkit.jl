function partial_state_selection_graph!(state::TransformationState)
    find_solvables!(state; allow_symbolic=true)
    var_eq_matching = complete(pantelides!(state))
    complete!(state.structure)
    partial_state_selection_graph!(state.structure, var_eq_matching)
end

function ascend_dg(xs, dg, level)
    while level > 0
        xs = Int[dg[x] for x in xs]
        level -= 1
    end
    return xs
end

function ascend_dg_all(xs, dg, level, maxlevel)
    r = Int[]
    while true
        if level <= 0
            append!(r, xs)
        end
        maxlevel <= 0 && break
        xs = Int[dg[x] for x in xs if dg[x] !== nothing]
        level -= 1
        maxlevel -= 1
    end
    return r
end

function pss_graph_modia!(structure::SystemStructure, var_eq_matching, varlevel, inv_varlevel, inv_eqlevel)
    @unpack eq_to_diff, var_to_diff, graph, solvable_graph = structure

    # var_eq_matching is a maximal matching on the top-differentiated variables.
    # Find Strongly connected components. Note that after pantelides, we expect
    # a balanced system, so a maximal matching should be possible.
    var_sccs::Vector{Union{Vector{Int}, Int}} = find_var_sccs(graph, var_eq_matching)
    var_eq_matching = Matching{Union{Unassigned, SelectedState}}(var_eq_matching)
    for vars in var_sccs
        # TODO: We should have a way to not have the scc code look at unassigned vars.
        if length(vars) == 1 && varlevel[vars[1]] != 0
            continue
        end

        # Now proceed level by level from lowest to highest and tear the graph.
        eqs = [var_eq_matching[var] for var in vars if var_eq_matching[var] !== unassigned]
        isempty(eqs) && continue
        maxlevel = level = maximum(map(x->inv_eqlevel[x], eqs))
        old_level_vars = ()
        ict = IncrementalCycleTracker(DiCMOBiGraph{true}(graph, complete(Matching(ndsts(graph)))); dir=:in)
        while level >= 0
            to_tear_eqs_toplevel = filter(eq->inv_eqlevel[eq] >= level, eqs)
            to_tear_eqs = ascend_dg(to_tear_eqs_toplevel, invview(eq_to_diff), level)

            to_tear_vars_toplevel = filter(var->inv_varlevel[var] >= level, vars)
            to_tear_vars = ascend_dg_all(to_tear_vars_toplevel, invview(var_to_diff), level, maxlevel)

            if old_level_vars !== ()
                # Inherit constraints from previous level.
                # TODO: Is this actually a good idea or do we want full freedom
                # to tear differently on each level? Does it make a difference
                # whether we're using heuristic or optimal tearing?
                removed_eqs = Int[]
                removed_vars = Int[]
                for var in old_level_vars
                    old_assign = ict.graph.matching[var]
                    if !isa(old_assign, Int) || ict.graph.matching[var_to_diff[var]] !== unassigned
                        continue
                    end
                    # Make sure the ict knows about this edge, so it doesn't accidentally introduce
                    # a cycle.
                    ok = try_assign_eq!(ict, var_to_diff[var], eq_to_diff[old_assign])
                    @assert ok
                    var_eq_matching[var_to_diff[var]] = eq_to_diff[old_assign]
                    push!(removed_eqs, eq_to_diff[ict.graph.matching[var]])
                    push!(removed_vars, var_to_diff[var])
                end
                to_tear_eqs = setdiff(to_tear_eqs, removed_eqs)
                to_tear_vars = setdiff(to_tear_vars, removed_vars)
            end
            filter!(var->ict.graph.matching[var] === unassigned, to_tear_vars)
            filter!(eq->invview(ict.graph.matching)[eq] === unassigned, to_tear_eqs)
            tearEquations!(ict, solvable_graph.fadjlist, to_tear_eqs, to_tear_vars)
            for var in to_tear_vars
                var_eq_matching[var] = ict.graph.matching[var]
            end
            old_level_vars = to_tear_vars
            level -= 1
        end
    end
    for var in 1:ndsts(graph)
        if varlevel[var] !== 0 && var_eq_matching[var] === unassigned
            var_eq_matching[var] = SelectedState()
        end
    end
    return var_eq_matching
end

struct SelectedState; end
function partial_state_selection_graph!(structure::SystemStructure, var_eq_matching)
    @unpack eq_to_diff, var_to_diff, graph, solvable_graph = structure
    eq_to_diff = complete(eq_to_diff)

    inv_eqlevel = map(1:nsrcs(graph)) do eq
        level = 0
        while invview(eq_to_diff)[eq] !== nothing
            eq = invview(eq_to_diff)[eq]
            level += 1
        end
        level
    end

    varlevel = map(1:ndsts(graph)) do var
        level = 0
        while var_to_diff[var] !== nothing
            var = var_to_diff[var]
            level += 1
        end
        level
    end

    inv_varlevel = map(1:ndsts(graph)) do var
        level = 0
        while invview(var_to_diff)[var] !== nothing
            var = invview(var_to_diff)[var]
            level += 1
        end
        level
    end

    # TODO: Should pantelides just return this?
    for var in 1:ndsts(graph)
        if var_to_diff[var] !== nothing
            var_eq_matching[var] = unassigned
        end
    end

    var_eq_matching = pss_graph_modia!(structure,
        complete(var_eq_matching), varlevel, inv_varlevel, inv_eqlevel)

    var_eq_matching
end
