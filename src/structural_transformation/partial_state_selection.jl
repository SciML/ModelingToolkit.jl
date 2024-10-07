function partial_state_selection_graph!(state::TransformationState)
    find_solvables!(state; allow_symbolic = true)
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

function pss_graph_modia!(structure::SystemStructure, maximal_top_matching, varlevel,
        inv_varlevel, inv_eqlevel)
    @unpack eq_to_diff, var_to_diff, graph, solvable_graph = structure

    # var_eq_matching is a maximal matching on the top-differentiated variables.
    # Find Strongly connected components. Note that after pantelides, we expect
    # a balanced system, so a maximal matching should be possible.
    var_sccs::Vector{Union{Vector{Int}, Int}} = find_var_sccs(graph, maximal_top_matching)
    var_eq_matching = Matching{Union{Unassigned, SelectedState}}(ndsts(graph))
    for vars in var_sccs
        # TODO: We should have a way to not have the scc code look at unassigned vars.
        if length(vars) == 1 && maximal_top_matching[vars[1]] === unassigned
            continue
        end

        # Now proceed level by level from lowest to highest and tear the graph.
        eqs = [maximal_top_matching[var]
               for var in vars if maximal_top_matching[var] !== unassigned]
        isempty(eqs) && continue
        maxeqlevel = maximum(map(x -> inv_eqlevel[x], eqs))
        maxvarlevel = level = maximum(map(x -> inv_varlevel[x], vars))
        old_level_vars = ()
        ict = IncrementalCycleTracker(
            DiCMOBiGraph{true}(graph,
                complete(Matching(ndsts(graph)), nsrcs(graph))),
            dir = :in)

        while level >= 0
            to_tear_eqs_toplevel = filter(eq -> inv_eqlevel[eq] >= level, eqs)
            to_tear_eqs = ascend_dg(to_tear_eqs_toplevel, invview(eq_to_diff), level)

            to_tear_vars_toplevel = filter(var -> inv_varlevel[var] >= level, vars)
            to_tear_vars = ascend_dg(to_tear_vars_toplevel, invview(var_to_diff), level)

            assigned_eqs = Int[]

            if old_level_vars !== ()
                # Inherit constraints from previous level.
                # TODO: Is this actually a good idea or do we want full freedom
                # to tear differently on each level? Does it make a difference
                # whether we're using heuristic or optimal tearing?
                removed_eqs = Int[]
                removed_vars = Int[]
                for var in old_level_vars
                    old_assign = var_eq_matching[var]
                    if isa(old_assign, SelectedState)
                        push!(removed_vars, var)
                        continue
                    elseif !isa(old_assign, Int) ||
                           ict.graph.matching[var_to_diff[var]] !== unassigned
                        continue
                    end
                    # Make sure the ict knows about this edge, so it doesn't accidentally introduce
                    # a cycle.
                    assgned_eq = eq_to_diff[old_assign]
                    ok = try_assign_eq!(ict, var_to_diff[var], assgned_eq)
                    @assert ok
                    var_eq_matching[var_to_diff[var]] = assgned_eq
                    push!(removed_eqs, eq_to_diff[ict.graph.matching[var]])
                    push!(removed_vars, var_to_diff[var])
                    push!(removed_vars, var)
                end
                to_tear_eqs = setdiff(to_tear_eqs, removed_eqs)
                to_tear_vars = setdiff(to_tear_vars, removed_vars)
            end
            tearEquations!(ict, solvable_graph.fadjlist, to_tear_eqs, BitSet(to_tear_vars),
                nothing)

            for var in to_tear_vars
                @assert var_eq_matching[var] === unassigned
                assgned_eq = ict.graph.matching[var]
                var_eq_matching[var] = assgned_eq
                isa(assgned_eq, Int) && push!(assigned_eqs, assgned_eq)
            end

            if level != 0
                remaining_vars = collect(v
                for v in to_tear_vars
                if var_eq_matching[v] === unassigned)
                if !isempty(remaining_vars)
                    remaining_eqs = setdiff(to_tear_eqs, assigned_eqs)
                    nlsolve_matching = maximal_matching(graph,
                        Base.Fix2(in, remaining_eqs),
                        Base.Fix2(in, remaining_vars))
                    for var in remaining_vars
                        if nlsolve_matching[var] === unassigned &&
                           var_eq_matching[var] === unassigned
                            var_eq_matching[var] = SelectedState()
                        end
                    end
                end
            end

            old_level_vars = to_tear_vars
            level -= 1
        end
    end
    return complete(var_eq_matching, nsrcs(graph))
end

struct SelectedState end
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
        graph_level = level = 0
        while var_to_diff[var] !== nothing
            var = var_to_diff[var]
            level += 1
            if !isempty(𝑑neighbors(graph, var))
                graph_level = level
            end
        end
        graph_level
    end

    inv_varlevel = map(1:ndsts(graph)) do var
        level = 0
        while invview(var_to_diff)[var] !== nothing
            var = invview(var_to_diff)[var]
            level += 1
        end
        level
    end

    var_eq_matching = pss_graph_modia!(structure,
        complete(var_eq_matching), varlevel, inv_varlevel,
        inv_eqlevel)

    var_eq_matching
end

function dummy_derivative_graph!(state::TransformationState, jac = nothing;
        state_priority = nothing, log = Val(false), kwargs...)
    state.structure.solvable_graph === nothing && find_solvables!(state; kwargs...)
    complete!(state.structure)
    var_eq_matching = complete(pantelides!(state; kwargs...))
    dummy_derivative_graph!(state.structure, var_eq_matching, jac, state_priority, log)
end

struct DummyDerivativeSummary
    var_sccs::Vector{Vector{Int}}
    state_priority::Vector{Vector{Float64}}
end

function dummy_derivative_graph!(
        structure::SystemStructure, var_eq_matching, jac = nothing,
        state_priority = nothing, ::Val{log} = Val(false)) where {log}
    @unpack eq_to_diff, var_to_diff, graph = structure
    diff_to_eq = invview(eq_to_diff)
    diff_to_var = invview(var_to_diff)
    invgraph = invview(graph)
    extended_sp = let state_priority = state_priority, var_to_diff = var_to_diff,
        diff_to_var = diff_to_var

        var -> begin
            min_p = max_p = 0.0
            while var_to_diff[var] !== nothing
                var = var_to_diff[var]
            end
            while true
                p = state_priority(var)
                max_p = max(max_p, p)
                min_p = min(min_p, p)
                (var = diff_to_var[var]) === nothing && break
            end
            min_p < 0 ? min_p : max_p
        end
    end

    var_sccs = find_var_sccs(graph, var_eq_matching)
    var_perm = Int[]
    var_dummy_scc = Vector{Int}[]
    var_state_priority = Vector{Float64}[]
    eqcolor = falses(nsrcs(graph))
    dummy_derivatives = Int[]
    col_order = Int[]
    neqs = nsrcs(graph)
    nvars = ndsts(graph)
    eqs = Int[]
    vars = Int[]
    next_eq_idxs = Int[]
    next_var_idxs = Int[]
    new_eqs = Int[]
    new_vars = Int[]
    eqs_set = BitSet()
    for vars′ in var_sccs
        empty!(eqs)
        empty!(vars)
        for var in vars′
            eq = var_eq_matching[var]
            eq isa Int || continue
            diff_to_eq[eq] === nothing || push!(eqs, eq)
            if var_to_diff[var] !== nothing
                error("Invalid SCC")
            end
            (diff_to_var[var] !== nothing && is_present(structure, var)) && push!(vars, var)
        end
        isempty(eqs) && continue

        rank_matching = Matching(max(nvars, neqs))
        isfirst = true
        if jac === nothing
            J = nothing
        else
            _J = jac(eqs, vars)
            # only accept small integers to avoid overflow
            is_all_small_int = all(_J) do x′
                x = unwrap(x′)
                x isa Number || return false
                isinteger(x) && typemin(Int8) <= x <= typemax(Int8)
            end
            J = is_all_small_int ? Int.(unwrap.(_J)) : nothing
        end
        while true
            nrows = length(eqs)
            iszero(nrows) && break

            if state_priority !== nothing && isfirst
                sp = extended_sp.(vars)
                resize!(var_perm, length(sp))
                sortperm!(var_perm, sp)
                permute!(vars, var_perm)
                permute!(sp, var_perm)
                push!(var_dummy_scc, copy(vars))
                push!(var_state_priority, sp)
            end
            # TODO: making the algorithm more robust
            # 1. If the Jacobian is a integer matrix, use Bareiss to check
            # linear independence. (done)
            #
            # 2. If the Jacobian is a single row, generate pivots. (Dynamic
            # state selection.)
            #
            # 3. If the Jacobian is a polynomial matrix, use Gröbner basis (?)
            if J !== nothing
                if !isfirst
                    J = J[next_eq_idxs, next_var_idxs]
                end
                N = ModelingToolkit.nullspace(J; col_order) # modifies col_order
                rank = length(col_order) - size(N, 2)
                for i in 1:rank
                    push!(dummy_derivatives, vars[col_order[i]])
                end
            else
                empty!(eqs_set)
                union!(eqs_set, eqs)
                rank = 0
                for var in vars
                    eqcolor .= false
                    # We need `invgraph` here because we are matching from
                    # variables to equations.
                    pathfound = construct_augmenting_path!(rank_matching, invgraph, var,
                        Base.Fix2(in, eqs_set), eqcolor)
                    pathfound || continue
                    push!(dummy_derivatives, var)
                    rank += 1
                    rank == nrows && break
                end
                fill!(rank_matching, unassigned)
            end
            if rank != nrows
                @warn "The DAE system is singular!"
            end

            # prepare the next iteration
            if J !== nothing
                empty!(next_eq_idxs)
                empty!(next_var_idxs)
            end
            empty!(new_eqs)
            empty!(new_vars)
            for (i, eq) in enumerate(eqs)
                ∫eq = diff_to_eq[eq]
                # descend by one diff level, but the next iteration of equations
                # must still be differentiated
                ∫eq === nothing && continue
                ∫∫eq = diff_to_eq[∫eq]
                ∫∫eq === nothing && continue
                if J !== nothing
                    push!(next_eq_idxs, i)
                end
                push!(new_eqs, ∫eq)
            end
            for (i, var) in enumerate(vars)
                ∫var = diff_to_var[var]
                ∫var === nothing && continue
                ∫∫var = diff_to_var[∫var]
                ∫∫var === nothing && continue
                if J !== nothing
                    push!(next_var_idxs, i)
                end
                push!(new_vars, ∫var)
            end
            eqs, new_eqs = new_eqs, eqs
            vars, new_vars = new_vars, vars
            isfirst = false
        end
    end

    if (n_diff_eqs = count(!isnothing, diff_to_eq)) !=
       (n_dummys = length(dummy_derivatives))
        @warn "The number of dummy derivatives ($n_dummys) does not match the number of differentiated equations ($n_diff_eqs)."
    end

    ret = tearing_with_dummy_derivatives(structure, BitSet(dummy_derivatives))
    if log
        (ret..., DummyDerivativeSummary(var_dummy_scc, var_state_priority))
    else
        ret[1]
    end
end

function is_present(structure, v)::Bool
    @unpack var_to_diff, graph = structure
    while true
        # if a higher derivative is present, then it's present
        isempty(𝑑neighbors(graph, v)) || return true
        v = var_to_diff[v]
        v === nothing && return false
    end
end

# Derivatives that are either in the dummy derivatives set or ended up not
# participating in the system at all are not considered differential
function is_some_diff(structure, dummy_derivatives, v)::Bool
    !(v in dummy_derivatives) && is_present(structure, v)
end

# We don't want tearing to give us `y_t ~ D(y)`, so we skip equations with
# actually differentiated variables.
function isdiffed((structure, dummy_derivatives), v)::Bool
    @unpack var_to_diff, graph = structure
    diff_to_var = invview(var_to_diff)
    diff_to_var[v] !== nothing && is_some_diff(structure, dummy_derivatives, v)
end

function tearing_with_dummy_derivatives(structure, dummy_derivatives)
    @unpack var_to_diff = structure
    # We can eliminate variables that are not selected (differential
    # variables). Selected unknowns are differentiated variables that are not
    # dummy derivatives.
    can_eliminate = falses(length(var_to_diff))
    for (v, dv) in enumerate(var_to_diff)
        dv = var_to_diff[v]
        if dv === nothing || !is_some_diff(structure, dummy_derivatives, dv)
            can_eliminate[v] = true
        end
    end
    var_eq_matching, full_var_eq_matching, var_sccs = tear_graph_modia(structure,
        Base.Fix1(isdiffed, (structure, dummy_derivatives)),
        Union{Unassigned, SelectedState};
        varfilter = Base.Fix1(getindex, can_eliminate))

    for v in 𝑑vertices(structure.graph)
        is_present(structure, v) || continue
        dv = var_to_diff[v]
        (dv === nothing || !is_some_diff(structure, dummy_derivatives, dv)) && continue
        var_eq_matching[v] = SelectedState()
    end

    return var_eq_matching, full_var_eq_matching, var_sccs, can_eliminate
end
