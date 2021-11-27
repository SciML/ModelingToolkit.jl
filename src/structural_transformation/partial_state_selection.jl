function partial_state_selection_graph!(sys::ODESystem)
    s = get_structure(sys)
    (s isa SystemStructure) || (sys = initialize_system_structure(sys))
    s = structure(sys)
    find_solvables!(sys)
    @set! s.graph = complete(s.graph)
    @set! sys.structure = s
    (sys, partial_state_selection_graph!(s.graph, s.solvable_graph, s.var_to_diff)...)
end

struct SelectedState; end
function partial_state_selection_graph!(graph, solvable_graph, var_to_diff)
    var_eq_matching, eq_to_diff = pantelides!(graph, solvable_graph, var_to_diff)
    eq_to_diff = complete(eq_to_diff)

    eqlevel = map(1:nsrcs(graph)) do eq
        level = 0
        while eq_to_diff[eq] !== nothing
            eq = eq_to_diff[eq]
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

    all_selected_states = Int[]

    level = 0
    level_vars = [var for var in 1:ndsts(graph) if varlevel[var] == 0 && invview(var_to_diff)[var] !== nothing]

    # TODO: Is this actually useful or should we just compute another maximal matching?
    for var in 1:ndsts(graph)
        if !(var in level_vars)
            var_eq_matching[var] = unassigned
        end
    end

    while level < maximum(eqlevel)
        var_eq_matching = tear_graph_modia(graph, solvable_graph;
            eqfilter  = eq->eqlevel[eq] == level && invview(eq_to_diff)[eq] !== nothing,
            varfilter = var->(var in level_vars && !(var in all_selected_states)))
        for var in level_vars
            if var_eq_matching[var] === unassigned
                selected_state = invview(var_to_diff)[var]
                push!(all_selected_states, selected_state)
                #=
                    # TODO: This is what the Matteson paper says, but it doesn't
                    #       quite seem to work.
                    while selected_state !== nothing
                        push!(all_selected_states, selected_state)
                        selected_state = invview(var_to_diff)[selected_state]
                    end
                =#
            end
        end
        level += 1
        level_vars = [var for var = 1:ndsts(graph) if varlevel[var] == level && invview(var_to_diff)[var] !== nothing]
    end

    var_eq_matching = tear_graph_modia(graph, solvable_graph;
        varfilter = var->!(var in all_selected_states))
    var_eq_matching = Matching{Union{Unassigned, SelectedState}}(var_eq_matching)
    for var in all_selected_states
        var_eq_matching[var] = SelectedState()
    end
    return var_eq_matching, eq_to_diff
end
