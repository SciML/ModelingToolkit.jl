using SymbolicUtils: Rewriters
using Graphs.Experimental.Traversals

alias_elimination(sys) = alias_elimination!(TearingState(sys))[1]

"""
    $TYPEDSIGNATURES

Efficiently construct an expression equivalent to
```julia
sum(cf * fullvars[v] for (cf, v) in zip(coeffs, vars))
```

`add_buffer` is a temporary buffer to avoid allocations.
"""
function build_expr_from_coeffs_vars!(
        add_buffer::Vector{SymbolicT}, coeffs::Vector, vars::Vector{Int},
        fullvars::Vector{SymbolicT}
    )
    resize!(add_buffer, length(coeffs))
    for (i, (cf, var)) in enumerate(zip(coeffs, vars))
        add_buffer[i] = cf * fullvars[var]
    end
    return SU.add_worker(VartypeT, add_buffer)
end

"""
    $TYPEDSIGNATURES

Convenience wrapper for `find_perfect_aliases!`.
"""
function eliminate_perfect_aliases!(state::TearingState)
    StateSelection.complete!(state.structure)
    eqs_to_rm = Int[]
    vars_to_rm = Int[]
    aliases = find_perfect_aliases!(state, eqs_to_rm, vars_to_rm)
    old_to_new_eq, old_to_new_var = StateSelection.rm_eqs_vars!(
        state, eqs_to_rm, vars_to_rm; eqs_sorted_and_uniqued = true
    )
    return nothing
end

"""
    $TYPEDSIGNATURES

Pick the target variable for an alias group. Irreducible variables must remain as
unknowns, so one of them is chosen as the target when the group contains any. Otherwise
the variable with the highest `state_priority` wins. When priorities are tied, prefer
the variable appearing in the most equations: visualization-only variables appear in a
single alias equation, while physics variables appear in many dynamics equations.
"""
function pick_alias_target(
        fullvars::Vector{SymbolicT}, group_vars::Vector{Int}, state_priorities, irreducibles::AtomicSetT,
        graph = nothing
    )
    irr_idx = findfirst(
        Base.Fix1(contains_possibly_indexed_element, irreducibles) ∘ Base.Fix1(getindex, fullvars),
        group_vars
    )
    irr_idx === nothing || return group_vars[irr_idx]
    max_priority = maximum(Base.Fix1(getindex, state_priorities), group_vars)
    candidates = filter(v -> state_priorities[v] == max_priority, group_vars)
    if graph !== nothing && length(candidates) > 1
        _, target_idx = findmax(v -> length(𝑑neighbors(graph, v)), candidates)
        return candidates[target_idx]
    end
    return candidates[1]
end

"""
    $TYPEDSIGNATURES

Identify variable aliases in `state`. Intended to act before `solvable_graph` is populated.
`eqs_to_rm` and `vars_to_rm` are buffers which will be appended with equations and
variables to remove from `state` and `mm`. Return `aliases::Dict{Int, Int}` indicating
which variables are aliased to some other variables. Keys of `aliases` are present in
`vars_to_rm`.
"""
function find_perfect_aliases!(
        state::TearingState, eqs_to_rm::Vector{Int}, vars_to_rm::Vector{Int}
    )
    (; sys, fullvars, structure) = state
    (; graph, solvable_graph, var_to_diff, state_priorities) = structure

    @assert solvable_graph === nothing
    diff_to_var = invview(var_to_diff)
    aliases = Dict{Int, Int}()
    subs = Dict{SymbolicT, SymbolicT}()
    eqs = collect(equations(state))
    original_eqs = state.original_eqs
    irreducibles = get_irreducibles(sys)

    # Not `IntDisjointSet` because we don't want singleton sets for every single variable
    alias_groups = DisjointSet{Int}()
    # Candidate alias equations `(ieq, v1_idx, v2_idx, edge_sign)`. `edge_sign == +1`
    # encodes `v1 ~ v2` (coefficients sum to zero); `edge_sign == -1` encodes
    # `v1 ~ -v2` (coefficients are equal). Removal is decided below once each group's
    # target is known: equations with a non-target irreducible endpoint must stay so
    # the remaining irreducibles are still constrained to the target.
    candidate_eqs = Tuple{Int, Int, Int, Int8}[]

    for ieq in 1:nsrcs(graph)
        snbors = 𝑠neighbors(graph, ieq)
        length(snbors) == 2 || continue
        diff_to_var[snbors[1]] === nothing || continue
        diff_to_var[snbors[2]] === nothing || continue

        eq = eqs[ieq]
        v1 = fullvars[snbors[1]]
        v2 = fullvars[snbors[2]]
        local edge_sign::Int8
        Moshi.Match.@match eq.rhs begin
            BSImpl.AddMul(; variant, coeff, dict) && if variant == SU.AddMulVariant.ADD end => begin
                SU._iszero(coeff) || continue
                length(dict) == 2 || continue
                haskey(dict, v1) && haskey(dict, v2) || continue
                if SU._iszero(dict[v1] + dict[v2])
                    edge_sign = Int8(1)   # v1 ~ v2
                elseif SU._iszero(dict[v1] - dict[v2])
                    edge_sign = Int8(-1)  # v1 ~ -v2
                else
                    continue
                end
            end
            _ => continue
        end
        push!(candidate_eqs, (ieq, snbors[1], snbors[2], edge_sign))
        push!(alias_groups, snbors[1])
        push!(alias_groups, snbors[2])
        union!(alias_groups, snbors[1], snbors[2])
    end

    alias_sets = Dict{Int, Vector{Int}}()
    sizehint!(alias_sets, DataStructures.num_groups(alias_groups))
    sizehint!(aliases, length(alias_groups) - DataStructures.num_groups(alias_groups))
    for var in alias_groups
        set = get!(() -> Int[], alias_sets, DataStructures.find_root!(alias_groups, var))
        push!(set, var)
    end

    group_target = Dict{Int, Int}()
    for (root, group_vars) in alias_sets
        group_target[root] = pick_alias_target(fullvars, group_vars, state_priorities, irreducibles, graph)
    end

    # Build per-group adjacency from the (signed) candidate edges and run a BFS from
    # each group's target to assign `var_sign[v] = ±1`, indicating whether `v` is
    # aliased to `+target` or `-target`. Sign propagation correctly handles
    # transitive chains, including double negation: `x ~ -y` and `y ~ -z` yield
    # `var_sign[z] = +1` (i.e. `z ~ x`). If a variable is reached via two paths
    # implying opposite signs (e.g. both `x ~ y` and `x ~ -y` are present), the
    # group is inconsistent: every variable in it is implied to be zero. We record
    # such groups in `conflict_groups`; below we substitute their non-irreducibles
    # with the symbolic literal `0` and rewrite one alias equation per irreducible
    # to `irr ~ 0`.
    group_adj = Dict{Int, Vector{Tuple{Int, Int8}}}()
    for (_, v1, v2, s) in candidate_eqs
        push!(get!(() -> Tuple{Int, Int8}[], group_adj, v1), (v2, s))
        push!(get!(() -> Tuple{Int, Int8}[], group_adj, v2), (v1, s))
    end
    var_sign = Dict{Int, Int8}()
    conflict_groups = Set{Int}()
    bfs_queue = Int[]
    for (root, _) in alias_sets
        target = group_target[root]
        var_sign[target] = Int8(1)
        empty!(bfs_queue)
        push!(bfs_queue, target)
        while !isempty(bfs_queue)
            u = popfirst!(bfs_queue)
            us = var_sign[u]
            nbrs = get(group_adj, u, nothing)
            nbrs === nothing && continue
            for (n, s) in nbrs
                implied = Int8(s * us)
                if haskey(var_sign, n)
                    if var_sign[n] != implied
                        push!(conflict_groups, root)
                    end
                else
                    var_sign[n] = implied
                    push!(bfs_queue, n)
                end
            end
        end
    end

    is_sticky = (v) -> contains_possibly_indexed_element(irreducibles, fullvars[v]) || state_priorities[v] > 0
    is_irreducible_v = (v) -> contains_possibly_indexed_element(irreducibles, fullvars[v])

    # Symbolic literal zero used as the substitution target for variables in conflict
    # (sign-inconsistent) groups.
    zero_sym = Symbolics.COMMON_ZERO
    signed_sym = (s::Int8, x::SymbolicT) -> s == 1 ? x : (-1) * x

    # Conflict-group alias equations: in conflict groups every variable is
    # forced to `0`. Each irreducible must stay as an unknown, so claim one
    # alias equation per irreducible and rewrite it eagerly to `v ~ 0` (and
    # strip the bipartite edge to the other endpoint). A conflict group has
    # at least as many alias equations as vertices (it contains a cycle), so
    # there is always one equation per irreducible. Remaining alias equations
    # become `0 ~ 0` after the substitution pass and are queued for removal.
    # Doing the rewrite here -- before the main loop walks
    # `𝑑neighbors(graph, v)` for non-irreducibles -- keeps the pinned eq out
    # of `eqs_to_substitute`, so it survives untouched.
    let claimed = Set{Int}()
        for (ieq, v1, v2, _) in candidate_eqs
            root = DataStructures.find_root!(alias_groups, v1)
            root in conflict_groups || continue
            v_pin = if is_irreducible_v(v1) && !(v1 in claimed)
                v1
            elseif is_irreducible_v(v2) && !(v2 in claimed)
                v2
            else
                nothing
            end
            if v_pin === nothing
                push!(eqs_to_rm, ieq)
            else
                push!(claimed, v_pin)
                for u in copy(𝑠neighbors(graph, ieq))
                    u == v_pin && continue
                    Graphs.rem_edge!(graph, ieq, u)
                end
                eqs[ieq] = fullvars[v_pin] ~ zero_sym
                original_eqs[ieq] = fullvars[v_pin] ~ zero_sym
            end
        end
    end

    # Consistent-group alias equations: queue for removal only if both
    # endpoints collapse out after substitution (the equation becomes `0 ~ 0`).
    # An equation with a non-target irreducible endpoint must stay so the
    # remaining irreducible is still pinned to the target after substitution.
    for (ieq, v1, v2, _) in candidate_eqs
        root = DataStructures.find_root!(alias_groups, v1)
        root in conflict_groups && continue
        target = group_target[root]
        c1 = is_sticky(v1) ? v1 : target
        c2 = is_sticky(v2) ? v2 : target
        c1 == c2 && push!(eqs_to_rm, ieq)
    end

    eqs_to_substitute = Int[]
    for (root, group_vars) in alias_sets
        target = group_target[root]
        is_conflict = root in conflict_groups
        # Only meaningful for non-conflict groups: the chosen target survives as an
        # unknown. For conflict groups every non-irreducible variable (target
        # included) is replaced by `0`, so we don't pin the target with
        # `always_present` here -- irreducibles get marked individually below.
        is_conflict || (state.always_present[target] = true)
        for v in group_vars
            !is_conflict && v == target && continue
            # In consistent groups, sticky vars (irreducible OR priority>0) other
            # than the target stay as unknowns. In conflict groups, ONLY truly
            # irreducible vars survive as unknowns; priority>0 non-irreducibles
            # are eliminated along with everything else. Irreducibles in conflict
            # groups already had one alias equation rewritten to `v ~ 0` above,
            # so we just mark them present here -- adding them to `subs` would
            # clobber that pin.
            if is_conflict ? is_irreducible_v(v) : is_sticky(v)
                state.always_present[v] = true
                continue
            end

            if is_conflict
                push!(vars_to_rm, v)
                subs[fullvars[v]] = zero_sym
                push!(state.additional_observed, fullvars[v] ~ zero_sym)

                dnbors = copy(𝑑neighbors(graph, v))
                for e in dnbors
                    push!(eqs_to_substitute, e)
                    Graphs.rem_edge!(graph, e, v)
                end

                dv = var_to_diff[v]
                while dv isa Int
                    push!(vars_to_rm, dv)
                    subs[fullvars[dv]] = zero_sym

                    dnbors = copy(𝑑neighbors(graph, dv))
                    for e in dnbors
                        push!(eqs_to_substitute, e)
                        Graphs.rem_edge!(graph, e, dv)
                    end

                    dv = var_to_diff[dv]
                end
            else
                s = var_sign[v]
                push!(vars_to_rm, v)
                rhs_sym = signed_sym(s, fullvars[target])
                subs[fullvars[v]] = rhs_sym
                push!(state.additional_observed, fullvars[v] ~ rhs_sym)
                aliases[v] = target

                dnbors = copy(𝑑neighbors(graph, v))
                for e in dnbors
                    push!(eqs_to_substitute, e)
                    Graphs.rem_edge!(graph, e, v)
                    Graphs.add_edge!(graph, e, target)
                end

                dv = var_to_diff[v]
                dtarget = var_to_diff[target]
                while dv isa Int
                    if dtarget === nothing
                        dtarget = StateSelection.var_derivative!(state, target)
                    end
                    push!(vars_to_rm, dv)

                    subs[fullvars[dv]] = signed_sym(s, fullvars[dtarget])
                    aliases[dv] = dtarget

                    dnbors = copy(𝑑neighbors(graph, dv))
                    for e in dnbors
                        push!(eqs_to_substitute, e)
                        Graphs.rem_edge!(graph, e, dv)
                        Graphs.add_edge!(graph, e, dtarget)
                    end

                    dv = var_to_diff[dv]
                    dtarget = var_to_diff[dtarget]
                end
            end
        end
    end

    # We need to handle unscalarized array variables
    for k in keys(subs)
        k, isarr = split_indexed_var(k)
        isarr || continue
        haskey(subs, k) && continue
        subs[k] = Symbolics.SConst(collect(k))
    end

    ir = get_irstructure(sys)
    subber = SU.IRSubstituter{false}(ir, subs)
    for e in eqs_to_substitute
        # Double substitute to handle unscalarized array variables. First one
        # substitutes `x` to `[x[1], x[2]]`. The second substitutes `x[1]` and `x[2]`.
        eqs[e] = subber(subber(eqs[e]))
        original_eqs[e] = subber(subber(original_eqs[e]))
    end

    # After substitution, alias equations that connected a sticky non-target
    # variable to a zero-priority variable become structurally identical to the
    # direct alias between the sticky variable and the target (since the
    # zero-priority variable was redirected to the target in the graph). Remove
    # all but the first copy of each (v_a, v_b) variable pair.
    let seen = Set{Tuple{Int,Int}}()
        eqs_rm_set = Set(eqs_to_rm)
        removed_additional_eqs = false
        for (ieq, _, _, _) in candidate_eqs
            ieq in eqs_rm_set && continue
            snbors = 𝑠neighbors(graph, ieq)
            length(snbors) == 2 || continue
            pair = minmax(snbors[1], snbors[2])
            if pair in seen
                removed_additional_eqs = true
                push!(eqs_to_rm, ieq)
            else
                push!(seen, pair)
            end
        end
        removed_additional_eqs && sort!(eqs_to_rm)
    end

    @set! sys.eqs = eqs
    state.sys = sys

    return aliases
end

function alias_elimination!(state::TearingState; fully_determined = true,
                            print_underconstrained_variables = false, kwargs...)
    StateSelection.complete!(state.structure)
    eqs_to_rm = Int[]
    vars_to_rm = Int[]
    add_buffer = SymbolicT[]
    aliases = Dict{Int, Int}()

    (; fullvars, sys) = state
    (; graph, var_to_diff, solvable_graph) = state.structure
    # Previously, underconstrained variables were zeroed out. This leads to significant
    # unintuitive behavior, especially for initialization systems. The underconstrained variables
    # should constitute an underdetermined error, and hence the behavior is removed. This pass
    # continues to remove redundant equations, since it is essential for adding analysis points
    # to existing connections.
    variable_underconstrained! = IgnoreUnderconstrainedVariable()
    mm = StateSelection.structural_singularity_removal!(state; variable_underconstrained!, kwargs...)

    if print_underconstrained_variables
        underconstrained_vars = state.fullvars[variable_underconstrained!.underconstrained]
        @info "Found underconstrained variables in the system" underconstrained_vars
    end

    eqs = collect(equations(state))
    fullvars_to_idx = Dict{SymbolicT, Int}(Iterators.map(reverse, enumerate(fullvars)))
    original_eqs = state.original_eqs

    for (ieq, eq) in enumerate(mm.nzrows)
        rcol = mm.row_cols[ieq]
        rval = mm.row_vals[ieq]
        if isempty(rcol)
            push!(eqs_to_rm, eq)
            continue
        end

        rhs = build_expr_from_coeffs_vars!(add_buffer, rval, rcol, fullvars)
        orig = eqs[eq]
        eqs[eq] = Symbolics.COMMON_ZERO ~ rhs
        oeq = original_eqs[eq]
        # NOTE: For discrete systems, `original_eqs` isn't shifted forward by 1
        # like everything else, and we need to maintain that invariant here too.
        if StateSelection.is_only_discrete(state.structure)
            iv = get_iv(sys)::SymbolicT
            lhs = MTKBase.simplify_shifts(Shift(iv, 1)(oeq.lhs))
        else
            lhs = oeq.lhs
        end
        if (idx = get(fullvars_to_idx, lhs, nothing); idx isa Int) && (colidx = findfirst(isequal(idx), rcol); colidx isa Int) && rval[colidx] == -1
            if StateSelection.is_only_discrete(state.structure)
                iv = get_iv(sys)::SymbolicT
                original_eqs[eq] = oeq.lhs ~ MTKTearing.backshift_expr(rhs, iv) + oeq.lhs
            else
                original_eqs[eq] = oeq.lhs ~ rhs + oeq.lhs
            end
        else
            if StateSelection.is_only_discrete(state.structure)
                iv = get_iv(sys)::SymbolicT
                original_eqs[eq] = MTKTearing.backshift_expr(eqs[eq], iv)
            else
                original_eqs[eq] = eqs[eq]
            end
        end
    end

    @set! sys.eqs = eqs
    state.sys = sys

    old_to_new_eq, old_to_new_var = StateSelection.rm_eqs_vars!(state, eqs_to_rm, vars_to_rm)
    sys = state.sys
    mm = StateSelection.get_new_mm(aliases, old_to_new_eq, old_to_new_var, mm)
    # This phrasing infers the return type as `Union{Tuple{...}}` instead of
    # `Tuple{Union{...}, ...}`
    if mm isa CLIL.SparseMatrixCLIL{BigInt, Int}
        return invalidate_cache!(sys), mm
    else
        return invalidate_cache!(sys), mm
    end
end

struct IgnoreUnderconstrainedVariable
    underconstrained::Vector{Int}
end

IgnoreUnderconstrainedVariable() = IgnoreUnderconstrainedVariable(Int[])

function (iuv::IgnoreUnderconstrainedVariable)(structure::SystemStructure, ils::CLIL.SparseMatrixCLIL, v::Int)
    push!(iuv.underconstrained, v)
    return ils
end

function exactdiv(a::Integer, b)
    d, r = divrem(a, b)
    @assert r == 0
    return d
end

swap!(v, i, j) = v[i], v[j] = v[j], v[i]
