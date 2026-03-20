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
    empty!(add_buffer)
    for (cf, var) in zip(coeffs, vars)
        push!(add_buffer, cf * fullvars[var])
    end
    return SU.add_worker(VartypeT, add_buffer)
end

"""
    $TYPEDSIGNATURES

Identify variable aliases in `state`. Intended to act before `solvable_graph` is populated. `add_buffer`
is a temporary mutable buffer to avoid allocations. `eqs_to_rm` and `vars_to_rm` are buffers
which will be appended with equations and variables to remove from `state` and `mm`.
Return `aliases::Dict{Int, Int}` indicating which variables are aliased to some other
variables. Keys of `aliases` are present in `vars_to_rm`.
"""
function find_perfect_aliases!(
        state::TearingState, add_buffer::Vector{SymbolicT},
        eqs_to_rm::Vector{Int}, vars_to_rm::Vector{Int}
    )
    (; sys, fullvars, structure) = state
    (; graph, solvable_graph, var_to_diff) = structure

    @assert solvable_graph === nothing
    diff_to_var = invview(var_to_diff)
    aliases = Dict{Int, Int}()
    subs = Dict{SymbolicT, SymbolicT}()
    eqs = collect(equations(state))
    original_eqs = state.original_eqs

    # Not `IntDisjointSet` because we don't want singleton sets for every single variable
    alias_groups = DisjointSet{Int}()

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

    for aset in values(alias_sets)
        _, target_idx = findmax(Base.Fix1(getindex, state.structure.state_priorities), aset)
        target = aset[target_idx]
        state.always_present[target] = true
        for v in aset
            v == target && continue
            push!(vars_to_rm, v)
            subs[fullvars[v]] = fullvars[target]
            push!(state.additional_observed, fullvars[v] ~ fullvars[target])
            aliases[v] = target

            dnbors = copy(𝑑neighbors(graph, v))
            alias_eq = 0
            for e in dnbors
                snbors = 𝑠neighbors(graph, e)
                if alias_eq == 0 && length(snbors) == 2 && (snbors[1] == target || snbors[2] == target)
                    push!(eqs_to_rm, e)
                    alias_eq = e
                    continue
                end
                eqs[e] = substitute(eqs[e], subs)
                original_eqs[e] = substitute(original_eqs[e], subs)
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
                    if length(snbors) == 2 && (snbors[1] == dtarget || snbors[2] == dtarget)
                        push!(eqs_to_rm, e)
                        alias_eq = e
                        continue
                    end
                    eqs[e] = substitute(eqs[e], subs)
                    original_eqs[e] = substitute(original_eqs[e], subs)
                    Graphs.rem_edge!(graph, e, dv)
                    Graphs.add_edge!(graph, e, dtarget)
                end

                dv = var_to_diff[dv]
                dtarget = var_to_diff[dtarget]
            end
        end
    end

    @set! sys.eqs = eqs
    state.sys = sys

    return aliases
end

function alias_elimination!(state::TearingState; fully_determined = true,
                            print_underconstrained_variables = false, kwargs...)
    sys = state.sys
    StateSelection.complete!(state.structure)
    eqs_to_rm = Int[]
    vars_to_rm = Int[]
    add_buffer = SymbolicT[]
    aliases = find_perfect_aliases!(state, add_buffer, eqs_to_rm, vars_to_rm)

    old_to_new_eq, old_to_new_var = StateSelection.rm_eqs_vars!(state, eqs_to_rm, vars_to_rm)
    StateSelection.complete!(state.structure)
    mm = StateSelection.linear_subsys_adjmat!(state; kwargs...)
    if isempty(mm.nzrows)
        return sys, mm
    end

    empty!(eqs_to_rm)
    empty!(vars_to_rm)
    empty!(aliases)

    (; fullvars) = state
    (; graph, var_to_diff, solvable_graph) = state.structure
    # Previously, underconstrained variables were zeroed out. This leads to significant
    # unintuitive behavior, especially for intialization systems. The underconstrained variables
    # should constitute an underdetermined error, and hence the behavior is removed. This pass
    # continues to remove redundant equations, since it is essential for adding analysis points
    # to existing connections.
    variable_underconstrained! = IgnoreUnderconstrainedVariable()
    mm, pivotinfo = StateSelection.structural_singularity_removal!(state, mm, Val{true}(); variable_underconstrained!)
    for (ei, e) in enumerate(mm.nzrows)
        set_neighbors!(graph, e, mm.row_cols[ei])
    end
    if solvable_graph isa BipartiteGraph{Int, Nothing}
        for (ei, e) in enumerate(mm.nzrows)
            set_neighbors!(solvable_graph, e, mm.row_cols[ei])
        end
    end

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
