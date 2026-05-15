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

    diff_to_var = invview(var_to_diff)
    aliases = Dict{Int, Int}()
    subs = Dict{SymbolicT, SymbolicT}()
    eqs = collect(equations(state))
    original_eqs = state.original_eqs
    irreducibles = get_irreducibles(sys)

    # Not `IntDisjointSet` because we don't want singleton sets for every single variable
    alias_groups = DisjointSet{Int}()
    # Candidate alias equations `(ieq, v1_idx, v2_idx)`. Removal is decided below once
    # each group's target is known: equations with a non-target irreducible endpoint
    # must stay so the remaining irreducibles are still constrained to the target.
    candidate_eqs = Tuple{Int, Int, Int}[]

    for ieq in 1:nsrcs(graph)
        snbors = 𝑠neighbors(graph, ieq)
        length(snbors) == 2 || continue
        diff_to_var[snbors[1]] === nothing || continue
        diff_to_var[snbors[2]] === nothing || continue

        eq = eqs[ieq]
        v1 = fullvars[snbors[1]]
        v2 = fullvars[snbors[2]]
        Moshi.Match.@match eq.rhs begin
            BSImpl.AddMul(; variant, coeff, dict) && if variant == SU.AddMulVariant.ADD end => begin
                SU._iszero(coeff) || continue
                length(dict) == 2 || continue
                haskey(dict, v1) && haskey(dict, v2) && SU._iszero(dict[v1] + dict[v2]) || continue
            end
            _ => continue
        end
        push!(candidate_eqs, (ieq, snbors[1], snbors[2]))
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

    # Queue an alias equation for removal only if both endpoints collapse onto the
    # target after non-irreducibles are substituted -- i.e. the equation becomes
    # `T ~ T`. Any equation with a non-target irreducible endpoint is kept; when the
    # other endpoint is a non-irreducible, the existing substitution machinery below
    # rewrites the kept equation into `I ~ T` form automatically.
    for (ieq, v1, v2) in candidate_eqs
        target = group_target[DataStructures.find_root!(alias_groups, v1)]
        c1 = contains_possibly_indexed_element(irreducibles, fullvars[v1]) || state_priorities[v1] > 0 ? v1 : target
        c2 = contains_possibly_indexed_element(irreducibles, fullvars[v2]) || state_priorities[v2] > 0 ? v2 : target
        c1 == c2 && push!(eqs_to_rm, ieq)
    end

    eqs_to_substitute = Int[]
    for (root, group_vars) in alias_sets
        target = group_target[root]
        state.always_present[target] = true
        for v in group_vars
            v == target && continue
            # Irreducibles other than the target stay as unknowns; only non-irreducibles
            # are eliminated in favor of the target.
            if contains_possibly_indexed_element(irreducibles, fullvars[v]) || state_priorities[v] > 0
                state.always_present[v] = true
                continue
            end
            
            push!(vars_to_rm, v)
            subs[fullvars[v]] = fullvars[target]
            push!(state.additional_observed, fullvars[v] ~ fullvars[target])
            aliases[v] = target

            dnbors = copy(𝑑neighbors(graph, v))
            for e in dnbors
                snbors = 𝑠neighbors(graph, e)
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

                subs[fullvars[dv]] = fullvars[dtarget]
                aliases[dv] = dtarget

                dnbors = copy(𝑑neighbors(graph, dv))
                for e in dnbors
                    snbors = 𝑠neighbors(graph, e)
                    push!(eqs_to_substitute, e)
                    Graphs.rem_edge!(graph, e, dv)
                    Graphs.add_edge!(graph, e, dtarget)
                end

                dv = var_to_diff[dv]
                dtarget = var_to_diff[dtarget]
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
        for (ieq, _, _) in candidate_eqs
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
    # unintuitive behavior, especially for intialization systems. The underconstrained variables
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
