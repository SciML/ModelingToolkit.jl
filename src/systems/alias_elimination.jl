using SymbolicUtils: Rewriters
using Graphs.Experimental.Traversals

"""
    $TYPEDSIGNATURES

Return a system with perfect aliases eliminated.
"""
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
    return old_to_new_eq, old_to_new_var, aliases
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
the variable with the highest `state_priority` wins; all other group members are
eliminated in favour of the target. When positive priorities are tied, emit a warning and prefer the
variable appearing in the most equations (visualization-only variables appear in a single alias
equation, while physics variables appear in many dynamics equations); if that is also
tied, one of the tied variables is chosen arbitrarily.
"""
function pick_alias_target(
        state::TearingState, group_vars::Vector{Int}, irreducibles::AtomicSetT
    )
    (; fullvars, structure) = state
    (; graph, canonical_ranks, state_priorities) = structure
    irr_idx = findfirst(
        Base.Fix1(contains_possibly_indexed_element, irreducibles) ∘ Base.Fix1(getindex, fullvars),
        group_vars
    )
    irr_idx === nothing || return group_vars[irr_idx]
    max_priority = maximum(Base.Fix1(getindex, state_priorities), group_vars)
    candidates = filter(v -> state_priorities[v] == max_priority, group_vars)
    if length(candidates) > 1 && max_priority > 0
        if max_priority >= 100
            tied_names = getindex.(Ref(fullvars), candidates)
            @warn "Multiple variables in an alias group share the highest state_priority \
            ($max_priority); choosing alias target by equation count. Tied variables: $tied_names"
        end
        max_degree = maximum(v -> length(𝑑neighbors(graph, v)), candidates)
        filter!(v -> length(𝑑neighbors(graph, v)) == max_degree, candidates)
        sort!(candidates; by = Base.Fix1(getindex, canonical_ranks))
    end
    return candidates[1]
end

"""
    $TYPEDSIGNATURES

If the substitution rules `subrules` contain a rule for an indexed array variable
`v[i]`, add the rule `v => SConst(collect(v))`.
"""
function __add_unscalarized_array_subs!(subrules::Dict{SymbolicT, SymbolicT})
    for k in collect(keys(subrules))
        k, isarr = split_indexed_var(k)
        isarr || continue
        haskey(subrules, k) && continue
        # We could do `SConst(collect(k))` but `collect` is unstable and slow. We can build the
        # `array_literal` directly.
        args = Symbolics.SArgsT()
        sizehint!(args, length(k)::Int + 1)
        push!(args, Symbolics.SConst(size(k)))
        for i in SU.stable_eachindex(k)
            push!(args, k[i])
        end
        scal_k = Symbolics.STerm(SU.array_literal, args; type = SU.symtype(k), shape = SU.shape(k))
        subrules[k] = scal_k
    end

    return nothing
end

"""
    $TYPEDSIGNATURES

Perform the substitutions `subrules` on the subset of equations `eqs` and original equations
`original_eqs` indicated by `eqs_to_substitute` (an iterable of indices). `eqs` and `original_eqs`
should correspond to `𝑠vertices(state.structure.graph)`. Also update the incidence of
`state.structure.graph` and `state.structure.solvable_graph`.

NOTE: This assumes that the substitution in `state` does not introduce new incidence,
and that `𝑠neighbors(state.structure.graph, ieq)` after substitution is a subset of its value
before substitution.
"""
function __substitute_and_update_incidence!(
        eqs::Vector{Equation}, original_eqs::Vector{Equation}, state::TearingState,
        eqs_to_substitute, subrules::AbstractDict{SymbolicT, SymbolicT}
    )
    (; fullvars, sys, structure) = state
    (; graph, solvable_graph) = structure

    ir = get_irstructure(sys)
    subber = SU.IRSubstituter{true}(ir, subrules)
    new_vars = SU.IRStructureSearchBuffer(ir, Set{SymbolicT}())
    for e in eqs_to_substitute
        olde = eqs[e]
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
        if solvable_graph !== nothing
            filter!(v -> Graphs.has_edge(solvable_graph, BipartiteEdge(e, v)), new_row)
            BipartiteGraphs.set_neighbors!(solvable_graph, e, new_row)
        end
    end

    return nothing
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
                solvable_graph === nothing || BipartiteGraphs.set_neighbors!(BipartiteGraphs.invview(solvable_graph), v, ())

                dv = var_to_diff[v]
                while dv isa Int
                    push!(vars_to_rm, dv)
                    subs[fullvars[dv]] = zero_sym

                    append!(eqs_to_substitute, 𝑑neighbors(graph, dv))
                    BipartiteGraphs.set_neighbors!(BipartiteGraphs.invview(graph), dv, ())
                    solvable_graph === nothing || BipartiteGraphs.set_neighbors!(BipartiteGraphs.invview(solvable_graph), dv, ())

                    dv = var_to_diff[dv]
                end
            end
            continue
        end
        # For consistent groups pick a target (survives as unknown) and rebase
        # parities relative to it via `target_p`.
        target = pick_alias_target(state, group_vars, irreducibles)
        group_target[root] = target
        target_p = parity[target]
        for v in group_vars
            # In consistent groups, irreducible vars other than the target stay as unknowns.
            # All other non-target vars (including those with positive state_priority) are
            # eliminated in favour of the highest-priority target.
            if is_irreducible_v(v) || v == target
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
                if solvable_graph !== nothing && Graphs.has_edge(solvable_graph, BipartiteEdge(e, v))
                    Graphs.rem_edge!(solvable_graph, e, v)
                    Graphs.add_edge!(solvable_graph, e, target)
                end
            end

            dv = var_to_diff[v]
            # `prev_dtarget` must always be one level of differentiation below `dtarget`.
            prev_dtarget = target
            dtarget = var_to_diff[target]
            while dv isa Int
                if dtarget === nothing
                    dtarget = StateSelection.var_derivative!(state, prev_dtarget)
                end
                push!(vars_to_rm, dv)

                subs[fullvars[dv]] = signed_sym(s, fullvars[dtarget])
                aliases[dv] = dtarget

                for e in copy(𝑑neighbors(graph, dv))
                    push!(eqs_to_substitute, e)
                    Graphs.rem_edge!(graph, e, dv)
                    Graphs.add_edge!(graph, e, dtarget)
                    if solvable_graph !== nothing && Graphs.has_edge(solvable_graph, BipartiteEdge(e, dv))
                        Graphs.rem_edge!(solvable_graph, e, dv)
                        Graphs.add_edge!(solvable_graph, e, dtarget)
                    end
                end

                dv = var_to_diff[dv]
                prev_dtarget = dtarget
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
                if solvable_graph !== nothing
                    BipartiteGraphs.set_neighbors!(solvable_graph, ieq, [v_pin])
                end
                eqs[ieq] = fullvars[v_pin] ~ zero_sym
                original_eqs[ieq] = fullvars[v_pin] ~ zero_sym
            end
        else
            target = group_target[parent[v1]]
            c1 = is_irreducible_v(v1) ? v1 : target
            c2 = is_irreducible_v(v2) ? v2 : target
            c1 == c2 && push!(eqs_to_rm, ieq)
        end
    end

    # We need to handle unscalarized array variables
    __add_unscalarized_array_subs!(subs)

    __substitute_and_update_incidence!(eqs, original_eqs, state, eqs_to_substitute, subs)

    # After substitution, alias equations that connected a sticky non-target
    # variable to a zero-priority variable become structurally identical to the
    # direct alias between the sticky variable and the target (since the
    # zero-priority variable was redirected to the target in the graph). Remove
    # all but the first copy of each (v_a, v_b) variable pair.
    let seen = Set{Tuple{Int, Int}}()
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

function __filter_mm_eqs!(fn::F, mm::CLIL.SparseMatrixCLIL) where {F}
    ptr = firstindex(mm.nzrows)
    nnz = 0
    for i in eachindex(mm.nzrows)
        eqidx = mm.nzrows[i]
        fn(eqidx) && continue
        nnz += 1
        mm.nzrows[ptr] = mm.nzrows[i]
        mm.row_cols[ptr] = mm.row_cols[i]
        mm.row_vals[ptr] = mm.row_vals[i]
        ptr = nextind(mm.nzrows, ptr)
    end
    resize!(mm.nzrows, nnz)
    resize!(mm.row_cols, nnz)
    resize!(mm.row_vals, nnz)

    return mm
end

"""
    $TYPEDSIGNATURES

Runs `eliminate_zero_variables!` for a maximum of `maxiters` iterations. Returns the
updated `mm`.
"""
function eliminate_zero_variables_fixpoint!(state::TearingState, mm::CLIL.SparseMatrixCLIL; maxiters = 4, kwargs...)
    for i in 1:maxiters
        mm, modified = eliminate_zero_variables!(state, mm; kwargs...)
        modified || break
    end
    return mm
end

"""
    $TYPEDSIGNATURES

Create a parameter representing the value of `var` at `t = 0`. Requires that `var` is
of the form `x(t)` or `x(t)[i, j, ...]`.
"""
function __create_t0_parameter_for(var::SymbolicT)
    # FIXME: the parameter created here is of the form `x#0(t)`, but really it should be
    # `x#0`. The only reason it is declared with `(t)` is in case it is updated in a callback.
    # Callback affects are discrete systems, which require their unknowns to be `(t)`.
    # Realistically, callbacks use none of the discrete infrastructure and actually are
    # very similar to initialization systems. If we reuse that infrastructure, the
    # `(t)` here can be dropped.
    T = SU.symtype(var)
    sh = SU.shape(var)
    arrvar, isidx = split_indexed_var(var)
    result = Moshi.Match.@match arrvar begin
        BSImpl.Term(; f) && if f isa SymbolicT end => Moshi.Match.@match f begin
            BSImpl.Sym(; name) => begin
                io = IOBuffer()
                print(io, name)
                if isidx
                    idxs = @view(arguments(var)[2:end])
                    print(io, '[')
                    for idx in idxs
                        print(io, unwrap_const(idx)::Int, ", ")
                    end
                    seek(io, position(io) - 2)
                    truncate(io, position(io))
                    print(io, ']')
                end
                print(io, "#0")
                newname = Symbol(take!(io))
                toparam(
                    Symbolics.SSym(
                        newname; type = SU.FnType{Tuple, T, Nothing}, shape = sh
                    )(arguments(arrvar)[1])
                )
            end
        end
    end
    return result
end

"""
    $TYPEDSIGNATURES

If a row in `mm` only contains a single variable, that variable is identically zero.
Eliminate such any such zero variables which can be identified from `mm`. Return the
new `mm` and a boolean indicating if any variables were thus eliminated.
"""
function eliminate_zero_variables!(state::TearingState, mm::CLIL.SparseMatrixCLIL; allow_symbolic = false, allow_parameter = true, kws...)
    StateSelection.complete!(state.structure)
    (; sys, fullvars, structure, original_eqs) = state
    (; var_to_diff, graph) = structure
    diff_to_var = invview(var_to_diff)
    zero_vars, zero_eqs, eqs_to_substitute, protected_vars, protected_eqs = __eliminate_zero_variables!(state, mm)
    isempty(zero_vars) && return mm, false

    # Zero variables that are actually eliminated (not protected by irreducibility).
    unprotected_zero_vars = setdiff(zero_vars, protected_vars)
    # Equations we substitute the zeros into. Protected `var ~ 0` equations are retained
    # verbatim, so we must not substitute into them (that would turn them into `0 ~ 0`).
    setdiff!(eqs_to_substitute, protected_eqs)
    # If every zero variable is protected and none of them appears in any other equation,
    # there is nothing to eliminate or substitute, so the system is unchanged.
    isempty(unprotected_zero_vars) && isempty(eqs_to_substitute) && return mm, false

    eqs = collect(equations(state))
    aliases = Dict{Int, Int}()
    vars_to_rm = Int[]
    eqs_to_rm = Int[]
    # Substitute the zeros into all other equations
    subrules = Dict{SymbolicT, SymbolicT}()
    sizehint!(subrules, length(zero_vars))
    # Equations in `mm` to remove because they contain the integrated form of
    # a variable eliminated as zero.
    mm_eqs_to_rm = Set{Int}()
    # List all analytically integrated variables
    analytic_integ = Dict{SymbolicT, SymbolicT}()
    # Create new parameters representing the `t = 0` values of analytically integrated
    # variables. We do not reuse `Initial` parameters since they have their own distinct
    # semantics. `Initial` parameters take the value of the variable at `tspan[1]`. The
    # parameters added here represent the values at `t = 0`. For example, if `D(D(x)) ~ 0`
    # then `x ~ dx_0 * t + x0`. Here, `x0` (`dx_0`) is the value of `x` (`D(x)`) at
    # `t = 0` and not `tspan[1]`.
    new_ps = SymbolicT[]
    _guesses = copy(get_guesses(sys))
    for v in zero_vars
        sym = fullvars[v]
        ttsym = StateSelection.is_only_discrete(structure) ? sym : default_toterm(sym)
        subrules[sym] = Symbolics.COMMON_ZERO
        subrules[ttsym] = Symbolics.COMMON_ZERO
        # Protected variables (whose differential chain contains an irreducible variable)
        # are retained as unknowns along with their `var ~ 0` equation. We still substitute
        # their zero value into other equations (via `subrules`), but do not eliminate them,
        # mark them observed, or analytically integrate their antiderivatives.
        v in protected_vars && continue
        push!(state.additional_observed, ttsym ~ Symbolics.COMMON_ZERO)
        # Also need to handle corresponding integrated forms
        ∫var = diff_to_var[v]
        # If this is already the lowest order derivative, do nothing
        ∫var isa Int || continue
        iv = get_iv(sys)::SymbolicT
        D = if StateSelection.is_only_discrete(structure)
            Shift(iv, 1)
        else
            Differential(iv)
        end
        # Only integrated forms of the lowest order zero are polynomials in `t`.
        # E.g. if `D(x)` and `D(D(x))` are in `zero_vars` and we process `D(D(x))`
        # first, this check prevents us from writing `D(x) ~ Initial(D(x))`.
        if ∫var in zero_vars
            ∫sym = fullvars[∫var]
            tt∫sym = default_toterm(∫sym)
            # This must still be specified to correctly populate `schedule.dummy_sub`.
            state.analytical_derivatives[D(tt∫sym)] = Symbolics.COMMON_ZERO
            continue
        end
        rhs = Symbolics.COMMON_ZERO
        root_var = sym

        while ∫var isa Int
            # This variable needs to be removed too, though it's not zero
            push!(vars_to_rm, ∫var)
            sym = fullvars[∫var]
            ttsym = D isa Shift ? sym : default_toterm(sym)

            # This will populate `get_schedule(sys)` during reassembly. This allows us to
            # track derivative information even after eliminating it, and will add
            # `ttsym ~ rhs` as an initialization equation.
            state.analytical_derivatives[D(ttsym)] = rhs
            # Reuse the parameter derivative mechanism for specifying derivatives of these
            # variables. This prevents us from having to substitute every equation where this
            # is present. It also allows the initial condition of the variable to be solved for
            # from initial conditions of other observed variables.
            state.param_derivative_map[D(sym)] = rhs
            state.param_derivative_map[D(ttsym)] = rhs

            # This pattern generates the Horner form of the polynomial. We know all
            # `Initial` parameters will be present because `complete` adds `Initial`s
            # for every observable.
            t0_param = __create_t0_parameter_for(ttsym)
            analytic_integ[sym] = t0_param
            rhs = rhs * iv + t0_param
            push!(new_ps, t0_param)
            push!(state.additional_observed, ttsym ~ rhs)
            _guesses[t0_param] = sym

            nbors = 𝑑neighbors(graph, ∫var)
            # If any of them are in `mm`, they shouldn't be
            union!(mm_eqs_to_rm, nbors)

            v = ∫var
            ∫var = diff_to_var[∫var]
        end
    end
    # All `t0` parameters are solvable
    new_binds = copy(parent(get_bindings(sys)))
    for p in new_ps
        new_binds[p] = COMMON_MISSING
    end

    __filter_mm_eqs!(!in(mm_eqs_to_rm), mm)
    append!(vars_to_rm, unprotected_zero_vars)
    append!(eqs_to_rm, setdiff(zero_eqs, protected_eqs))
    __add_unscalarized_array_subs!(subrules)
    __substitute_and_update_incidence!(eqs, original_eqs, state, eqs_to_substitute, subrules)
    # It's possible that after substitution some higher order derivatives are not present
    # anymore. For example, `D(x)` might only be present in an equation as `y * D(x)`, and
    # we just substituted `y => 0`. Without the following, this will result in `D(x)` not being
    # incident on any equation. This causes dummy derivatives to complain, since it sees
    # `D(x)` as an unmatched highest differentiated variable, and thinks Pantelides made a
    # mistake.
    for v in 𝑑vertices(graph)
        StateSelection.is_present(structure, v) || push!(vars_to_rm, v)
    end

    sys = ConstructionBase.setproperties(
        sys;
        analytically_integrated = merge(analytically_integrated(sys), analytic_integ),
        ps = [get_ps(sys); new_ps], bindings = ROSymmapT(new_binds), eqs = eqs,
        guesses = _guesses
    )
    state.sys = sys

    old_to_new_eq, old_to_new_var = StateSelection.rm_eqs_vars!(state, eqs_to_rm, vars_to_rm)
    sys = state.sys
    mm = StateSelection.get_new_mm(aliases, old_to_new_eq, old_to_new_var, mm)

    eqs = equations(sys)

    # After substitution, some equations may now be integer-coefficient linear combinations,
    # and thus should be added to `mm`.
    mm_rows = BitSet(mm.nzrows)
    # Reuse buffers
    empty!(vars_to_rm)
    to_rm = vars_to_rm
    empty!(eqs_to_rm)
    coeffs = eqs_to_rm
    for e in eqs_to_substitute
        # Refer to equations by their new indices now
        e = old_to_new_eq[e]
        iszero(e) && continue
        # Ignore ones already in `mm`
        e in mm_rows && continue
        eq = eqs[e]
        all_int_vars, resid = StateSelection.find_eq_solvables!(state, e, to_rm, coeffs; allow_symbolic, allow_parameter)
        if all_int_vars && SU._iszero(resid)
            push!(mm.nzrows, e)
            push!(mm.row_cols, copy(𝑠neighbors(state.structure.solvable_graph, e)))
            push!(mm.row_vals, copy(coeffs))
        end
    end

    return mm, true
end

"""
    $TYPEDSIGNATURES

Helper for `eliminate_zero_variables!`. Returns:
- `zero_vars`: A `Set{Int}` of variables identified to be identically zero.
- `zero_eqs`: A `Set{Int}` of equations identifying such variables to be zero.
- `eqs_to_substitute`: A `Set{Int}` of equations in `state` which contain variables
  in `zero_vars`. These equations should be substituted to remove the zero variables.
- `protected_vars`: A `Set{Int}`, subset of `zero_vars`, whose differential chain (the
  variable together with all its derivatives and antiderivatives) contains an irreducible
  variable. These zero variables should still be substituted into other equations, but
  neither they nor their `var ~ 0` equations should be eliminated from the system.
- `protected_eqs`: A `Set{Int}`, subset of `zero_eqs`, of the `var ~ 0` equations
  identifying variables in `protected_vars` to be zero. These must be retained.

Also updates `mm` in-place to remove non-protected zero variables/equations.
"""
function __eliminate_zero_variables!(state::TearingState, mm::CLIL.SparseMatrixCLIL)
    (; graph, var_to_diff) = state.structure
    fullvars = state.fullvars
    diff_to_var = invview(var_to_diff)
    irreducibles = get_irreducibles(state.sys)

    # Inverse mapping of `mm.nzrows`
    nzrow_to_idx = Dict{Int, Int}()
    sizehint!(nzrow_to_idx, length(mm.nzrows))
    for (i, eqidx) in enumerate(mm.nzrows)
        nzrow_to_idx[eqidx] = i
    end
    # Variables we can identify are zero via `mm`
    zero_vars = Set{Int}()
    # Equations in `mm` of the form `var ~ 0`
    zero_eqs = Set{Int}()
    # Maps each zero variable to the `var ~ 0` equation (in `zero_eqs`) identifying it.
    zero_var_to_eq = Dict{Int, Int}()
    # Equations that `zero_vars` are incident on and need to be substituted
    eqs_to_substitute = Set{Int}()
    # Queue of indices in `mm.nzrows` containing rows to check for being
    # `var ~ 0`
    queue = Queue{Int}()
    # Initially check all rows with just one element
    for i in eachindex(mm.nzrows)
        isone(length(mm.row_cols[i])) && push!(queue, i)
    end

    process_neighbors! = let graph = graph, nzrow_to_idx = nzrow_to_idx,
            eqs_to_substitute = eqs_to_substitute, queue = queue
        function __process_neighbors!(zvar::Int)
            nbors = 𝑑neighbors(graph, zvar)
            union!(eqs_to_substitute, nbors)
            for nbor in nbors
                nzrow_idx = get(nzrow_to_idx, nbor, 0)
                # If `nbor` is not present in `mm`, it must be substituted later
                iszero(nzrow_idx) || push!(queue, nzrow_idx)
            end
            return
        end
    end

    while !isempty(queue)
        row_i = popfirst!(queue)
        eqidx = mm.nzrows[row_i]
        # Skip rows we already processed in case they show up twice
        eqidx in zero_eqs && continue

        # If a row only contains one non-zero element, that element is also zero.
        # We filter `rcol` and `rval` to remove any zero elements. This could be
        # done in the `for nbor in nbors` loop below. However, it's possible we process
        # two zero variables both present in the same row before processing that row.
        # This approach avoids the row being processed by each of the two zero variables
        # individually in their `nbors` loop. We delay removing all zero elements to
        # the latest possible time to aggregate updates.
        rcol = mm.row_cols[row_i]
        rval = mm.row_vals[row_i]
        ptr = firstindex(rcol)
        nnz = 0
        for i in eachindex(rcol)
            rcol[i] in zero_vars && continue
            rcol[ptr] = rcol[i]
            rval[ptr] = rval[i]
            ptr = nextind(rcol, ptr)
            nnz += 1
        end
        resize!(rcol, nnz)
        resize!(rval, nnz)
        isone(nnz) || continue

        zvar = rcol[1]

        push!(zero_vars, zvar)
        push!(zero_eqs, eqidx)
        zero_var_to_eq[zvar] = eqidx
        process_neighbors!(zvar)
        # All derivatives of this variable are also zero
        dzvar = var_to_diff[zvar]
        while dzvar isa Int
            push!(zero_vars, dzvar)
            process_neighbors!(dzvar)
            dzvar = var_to_diff[dzvar]
        end
    end


    # A zero variable is protected (retained as an unknown instead of being analytically
    # integrated) if its lowest order derivative is marked as irreducible, or if it (or any
    # variable in its differential chain) is modified by a discrete/continuous event. A
    # variable an event impulsively modifies is not a smooth function of time, so it must keep
    # being numerically integrated rather than being replaced by a polynomial in the iv.
    event_modified = MTKBase.variables_modified_by_events(state.sys)
    protected_vars = Set{Int}()
    if !isempty(irreducibles) || !isempty(event_modified)
        for v in zero_vars
            # Walk down to the lowest order antiderivative.
            root::Int = v
            while (∫root = diff_to_var[root]) isa Int
                root = ∫root
            end
            if contains_possibly_indexed_element(irreducibles, fullvars[root]) ||
                    contains_possibly_indexed_element(event_modified, fullvars[root]) ||
                    contains_possibly_indexed_element(event_modified, fullvars[v])
                push!(protected_vars, v)
            end
        end
    end
    # `var ~ 0` equations of protected variables must be retained. Only the lowest-order
    # zero of a chain has such an equation; higher derivatives are zeroed by inference.
    protected_eqs = Set{Int}()
    for v in protected_vars
        eq = get(zero_var_to_eq, v, 0)
        iszero(eq) || push!(protected_eqs, eq)
    end

    # Remove non-protected zero equations from `mm`. Protected `var ~ 0` equations are
    # retained since they still define the (retained) variable.
    eqs_to_remove = setdiff(zero_eqs, protected_eqs)
    isempty(eqs_to_remove) || __filter_mm_eqs!(!in(eqs_to_remove), mm)

    return zero_vars, zero_eqs, eqs_to_substitute, protected_vars, protected_eqs
end

function alias_elimination!(
        state::TearingState;
        print_underconstrained_variables = false, kwargs...
    )
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
