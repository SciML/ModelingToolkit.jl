using SymbolicUtils: Rewriters
using Graphs.Experimental.Traversals

alias_elimination(sys) = alias_elimination!(TearingState(sys))[1]
function alias_elimination!(state::TearingState; fully_determined = true, kwargs...)
    sys = state.sys
    StateSelection.complete!(state.structure)
    variable_underconstrained! = ZeroVariablesIfFullyDetermined!(fully_determined === true)
    mm = StateSelection.structural_singularity_removal!(state; variable_underconstrained!, kwargs...)

    fullvars = state.fullvars
    @unpack var_to_diff, graph, solvable_graph = state.structure

    subs = Dict{SymbolicT, SymbolicT}()
    # If we encounter y = -D(x), then we need to expand the derivative when
    # D(y) appears in the equation, so that D(-D(x)) becomes -D(D(x)).
    to_expand = Int[]
    diff_to_var = invview(var_to_diff)

    dels = Int[]
    eqs = collect(equations(state))
    resize!(eqs, nsrcs(graph))

    __trivial_eq_rhs = let fullvars = fullvars
        function trivial_eq_rhs(pair)
            var, coeff = pair
            iszero(coeff) && return Symbolics.COMMON_ZERO
            return coeff * fullvars[var]
        end
    end
    for (ei, e) in enumerate(mm.nzrows)
        vs = 𝑠neighbors(graph, e)
        if isempty(vs)
            # remove empty equations
            push!(dels, e)
        else
            rhs = mapfoldl(__trivial_eq_rhs, +, pairs(CLIL.nonzerosmap(@view mm[ei, :])))
            eqs[e] = Symbolics.COMMON_ZERO ~ rhs
        end
    end
    deleteat!(eqs, sort!(dels))
    old_to_new_eq = Vector{Int}(undef, nsrcs(graph))
    idx = 0
    cursor = 1
    ndels = length(dels)
    for i in eachindex(old_to_new_eq)
        if cursor <= ndels && i == dels[cursor]
            cursor += 1
            old_to_new_eq[i] = -1
            continue
        end
        idx += 1
        old_to_new_eq[i] = idx
    end

    n_new_eqs = idx

    eqs_to_update = BitSet()
    for ieq in eqs_to_update
        eq = eqs[ieq]
        eqs[ieq] = substitute(eq, subs)
    end
    new_nparentrows = nsrcs(graph)
    new_row_cols = eltype(mm.row_cols)[]
    new_row_vals = eltype(mm.row_vals)[]
    new_nzrows = Int[]
    for (i, eq) in enumerate(mm.nzrows)
        old_to_new_eq[eq] > 0 || continue
        push!(new_row_cols, mm.row_cols[i])
        push!(new_row_vals, mm.row_vals[i])
        push!(new_nzrows, old_to_new_eq[eq])
    end
    mm = typeof(mm)(new_nparentrows, mm.ncols, new_nzrows, new_row_cols, new_row_vals)

    for old_ieq in to_expand
        ieq = old_to_new_eq[old_ieq]
        eqs[ieq] = expand_derivatives(eqs[ieq])
    end

    diff_to_var = invview(var_to_diff)
    new_graph = BipartiteGraph(n_new_eqs, ndsts(graph))
    new_solvable_graph = BipartiteGraph(n_new_eqs, ndsts(graph))
    new_eq_to_diff = StateSelection.DiffGraph(n_new_eqs)
    eq_to_diff = state.structure.eq_to_diff
    for (i, ieq) in enumerate(old_to_new_eq)
        ieq > 0 || continue
        set_neighbors!(new_graph, ieq, 𝑠neighbors(graph, i))
        set_neighbors!(new_solvable_graph, ieq, 𝑠neighbors(solvable_graph, i))
        new_eq_to_diff[ieq] = eq_to_diff[i]
    end

    # update DiffGraph
    new_var_to_diff = StateSelection.DiffGraph(length(var_to_diff))
    for v in 1:length(var_to_diff)
        new_var_to_diff[v] = var_to_diff[v]
    end
    state.structure.graph = new_graph
    state.structure.solvable_graph = new_solvable_graph
    state.structure.eq_to_diff = new_eq_to_diff
    state.structure.var_to_diff = new_var_to_diff

    sys = state.sys
    @set! sys.eqs = eqs
    state.sys = sys
    # This phrasing infers the return type as `Union{Tuple{...}}` instead of
    # `Tuple{Union{...}, ...}`
    if mm isa CLIL.SparseMatrixCLIL{BigInt, Int}
        return invalidate_cache!(sys), mm
    else
        return invalidate_cache!(sys), mm
    end
end

struct ZeroVariablesIfFullyDetermined!
    fully_determined::Bool
end

function (zvifd::ZeroVariablesIfFullyDetermined!)(structure::SystemStructure, ils::CLIL.SparseMatrixCLIL, v::Int)
    return zvifd.fully_determined ? StateSelection.force_var_to_zero!(structure, ils, v) : ils
end

function exactdiv(a::Integer, b)
    d, r = divrem(a, b)
    @assert r == 0
    return d
end

swap!(v, i, j) = v[i], v[j] = v[j], v[i]
