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

Weighted union-find augmented with sign tracking. Merges the components of `v1`
and `v2` under the relation `v1 ~ edge_sign · v2` (`edge_sign ∈ ±1`)
maintains the invariant that `parent[v]` is always the current root and `parity[v]`
is the sign of `v` relative to that root so finds are O(1).

Merges are smaller-into-larger using `members` (root → member list)
giving O(E log V) across all unions.
A contradiction (e.g. `x ~ y, x ~ -y`) zeros out every member's `parity` in the
affected component; downstream code treats `parity[v] == 0` as "v is forced to 0".
"""
function union_with_sign!(
        parent::Dict{Int, Int}, parity::Dict{Int, Int8},
        members::Dict{Int, Vector{Int}},
        v1::Int, v2::Int, edge_sign::Int8
    )
    for v in (v1, v2)
        if !haskey(parent, v)
            parent[v] = v
            parity[v] = Int8(1)
            members[v] = [v]
        end
    end
    r1 = parent[v1]; s1 = parity[v1]
    r2 = parent[v2]; s2 = parity[v2]
    if r1 == r2
        if s1 != Int8(edge_sign * s2)
            for m in members[r1]
                parity[m] = Int8(0)
            end
        end
        return
    end
    if length(members[r1]) < length(members[r2])
        r1, r2 = r2, r1
        s1, s2 = s2, s1
    end
    # Edge `v1 = edge_sign · v2` ⇒ sign of r2 relative to r1 = s1 · edge_sign · s2.
    # If either side is already zero (conflict-tainted), `r2_to_r1` is 0 and we
    # propagate the zero through both components.
    r2_to_r1 = Int8(s1 * edge_sign * s2)
    r2_members = members[r2]
    for m in r2_members
        parent[m] = r1
        parity[m] = Int8(parity[m] * r2_to_r1)
    end
    if r2_to_r1 == 0
        for m in members[r1]
            parity[m] = Int8(0)
        end
    end
    append!(members[r1], r2_members)
    delete!(members, r2)
    return
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

    # Weighted union-find tracking each variable's sign relative to its
    # component's current root. `parent[v]` is always the root.
    # `parity[v] ∈ ±1` is the sign of `v` relative to that root, or `0` if the
    # component contains a sign contradiction (forcing every member to 0).
    # `members` maps root → member list (used for smaller-to-larger merges).
    parent = Dict{Int, Int}()
    parity = Dict{Int, Int8}()
    members = Dict{Int, Vector{Int}}()
    # Candidate alias equations `(ieq, v1_idx, v2_idx, edge_sign)`. 
    # `edge_sign == ±1` encodes `v1 ~ edge_sign*v2`.
    # Removal is decided below once each group's target is known:
    # equations with a non-target irreducible endpoint must stay so
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
        union_with_sign!(parent, parity, members, snbors[1], snbors[2], edge_sign)
    end
    sizehint!(aliases, length(parent) - length(members))

    # After the scan above, every variable in any alias group has
    # `parent[v] == root` and `parity[v] ∈ {-1, 0, +1}`
    # (sign relative to root, or 0 if the component is conflict-tainted).
    # The main rewriting loop below picks a target per non-conflict group
    # and uses `s = parity[v] * parity[target]` (±1) as `v`'s sign relative to the target.
    # For Conflict groups: every non-irreducible gets substituted with the
    # symbolic `0`, and one alias eq per irreducible is rewritten to `irr ~ 0`.
    group_target = Dict{Int, Int}()
    sizehint!(group_target, length(members))

    is_sticky = (v) -> contains_possibly_indexed_element(irreducibles, fullvars[v]) || state_priorities[v] > 0
    is_irreducible_v = (v) -> contains_possibly_indexed_element(irreducibles, fullvars[v])

    # Symbolic zero used as substitution target for variables in conflict groups.
    zero_sym = Symbolics.COMMON_ZERO
    signed_sym = (s::Int8, x::SymbolicT) -> s == 1 ? x : (-1) * x

    eqs_to_substitute = Int[]
    irrs_by_root = Dict{Int, Vector{Int}}()
    for (root, group_vars) in members
        if parity[root] == 0 # is conflict group
            # collect the irreducibles that need an equation forcing them to `0`.
            irrs_by_root[root] = filter(is_irreducible_v, group_vars)
            for v in group_vars
                # In conflict groups, ONLY irreducible vars survive as unknowns.
                # The per-equation cleanup below rewrites one alias eq per
                # irreducible to `v ~ 0`; here we just mark them present.
                if is_irreducible_v(v)
                    state.always_present[v] = true
                    continue
                end
                push!(vars_to_rm, v)
                subs[fullvars[v]] = zero_sym
                push!(state.additional_observed, fullvars[v] ~ zero_sym)

                append!(eqs_to_substitute, 𝑑neighbors(graph, v))
                BipartiteGraphs.set_neighbors!(BipartiteGraphs.invview(graph), v, ())

                dv = var_to_diff[v]
                while dv isa Int
                    push!(vars_to_rm, dv)
                    subs[fullvars[dv]] = zero_sym

                    append!(eqs_to_substitute, 𝑑neighbors(graph, dv))
                    BipartiteGraphs.set_neighbors!(BipartiteGraphs.invview(graph), dv, ())

                    dv = var_to_diff[dv]
                end
            end
            continue
        end
        # For consistent groups pick a target (survives as unknown) and rebase
        # parities relative to it via `target_p`.
        target = pick_alias_target(fullvars, group_vars, state_priorities, irreducibles, graph)
        group_target[root] = target
        target_p = parity[target]
        for v in group_vars
            # In consistent groups sticky vars (irreducible OR priority>0) other
            # than the target stay as unknowns.
            if is_sticky(v) || v == target
                state.always_present[v] = true
                continue
            end

            s = Int8(parity[v] * target_p)
            push!(vars_to_rm, v)
            rhs_sym = signed_sym(s, fullvars[target])
            subs[fullvars[v]] = rhs_sym
            push!(state.additional_observed, fullvars[v] ~ rhs_sym)
            aliases[v] = target

            for e in copy(𝑑neighbors(graph, v))
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

                for e in copy(𝑑neighbors(graph, dv))
                    push!(eqs_to_substitute, e)
                    Graphs.rem_edge!(graph, e, dv)
                    Graphs.add_edge!(graph, e, dtarget)
                end

                dv = var_to_diff[dv]
                dtarget = var_to_diff[dtarget]
            end
        end
    end

    # Per-equation cleanup. Two cases:
    #
    # Conflict groups: every variable is forced to `0`.
    # Take the first `length(irrs)` eq in the group to force the irrs to 0
    # (overwriting the eq and replacing its graph row with the edge to `irr`). 
    # Remaining conflict eqs become `0 ~ 0` and are queued for removal.
    # Pinned eqs end up in `eqs_to_substitute` from the var-elim pass above
    # but `subber` leaves `irr ~ 0` untouched (irreducibles aren't in `subs`).
    #
    # Consistent groups: queue for removal iff both endpoints collapse. An eq
    # with a non-target irreducible endpoint must stay so the remaining
    # irreducible is still pinned to the target after substitution.
    for (ieq, v1, v2, _) in candidate_eqs
        if parity[v1] == 0 # if conflict group
            irrs = irrs_by_root[parent[v1]]
            if isempty(irrs)
                push!(eqs_to_rm, ieq)
            else
                v_pin = pop!(irrs) # irreducible var to pin to 0
                BipartiteGraphs.set_neighbors!(graph, ieq, [v_pin])
                eqs[ieq] = fullvars[v_pin] ~ zero_sym
                original_eqs[ieq] = fullvars[v_pin] ~ zero_sym
            end
        else
            target = group_target[parent[v1]]
            c1 = is_sticky(v1) ? v1 : target
            c2 = is_sticky(v2) ? v2 : target
            c1 == c2 && push!(eqs_to_rm, ieq)
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
    subber = SU.IRSubstituter{true}(ir, subs)
    new_vars = SU.IRStructureSearchBuffer(ir, Set{SymbolicT}())
    for e in eqs_to_substitute
        # Double substitute to handle unscalarized array variables. First one
        # substitutes `x` to `[x[1], x[2]]`. The second substitutes `x[1]` and `x[2]`.
        eqs[e] = subber(subber(eqs[e]))
        original_eqs[e] = subber(subber(original_eqs[e]))
        # Substitution can drop vars beyond the one substituted: zero
        # substitution annihilates multiplicative cofactors (`v*w` with `v→0`
        # also removes `w`); alias substitution can cancel the target
        # (`v - target` with `v→target` leaves neither). Substitution never
        # introduces new vars, so the new incidence is a subset of the old —
        # prune the row to entries actually present in the simplified RHS.
        empty!(new_vars)
        SU.search_variables!(new_vars, eqs[e])
        new_row = filter(
            v_idx -> fullvars[v_idx] in new_vars || split_indexed_var(fullvars[v_idx])[1] in new_vars,
            𝑠neighbors(graph, e)
        )
        BipartiteGraphs.set_neighbors!(graph, e, new_row)
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
