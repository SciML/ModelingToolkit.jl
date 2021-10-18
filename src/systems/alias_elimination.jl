using SymbolicUtils: Rewriters

const KEEP = typemin(Int)

include("compat/bareiss.jl")

function alias_elimination(sys)
    sys = initialize_system_structure(sys; quick_cancel=true)
    s = structure(sys)
    is_linear_equations, eadj, cadj = find_linear_equations(sys)

    v_eliminated, v_types, n_null_vars, degenerate_equations, linear_equations = alias_eliminate_graph(
        s, is_linear_equations, eadj, cadj
    )

    s = structure(sys)
    @unpack fullvars, graph = s

    n_reduced_states = length(v_eliminated) - n_null_vars
    subs = OrderedDict()
    if n_reduced_states > 0
        for (i, v) in enumerate(@view v_eliminated[n_null_vars+1:end])
            subs[fullvars[v]] = iszeroterm(v_types, v) ? 0.0 :
                                isalias(v_types, v) ? fullvars[alias(v_types, v)] :
                                -fullvars[negalias(v_types, v)]
        end
    end

    dels = Set{Int}()
    eqs = copy(equations(sys))
    for (ei, e) in enumerate(linear_equations)
        vs = 𝑠neighbors(graph, e)
        if isempty(vs)
            push!(dels, e)
        else
            rhs = 0
            for vj in eachindex(vs)
                var = fullvars[vs[vj]]
                rhs += cadj[ei][vj] * var
            end
            eqs[e] = 0 ~ rhs
        end
    end
    dels = sort(collect(dels))
    deleteat!(eqs, dels)

    dict = Dict(subs)
    for (ieq, eq) in enumerate(eqs)
        eqs[ieq] = eq.lhs ~ fixpoint_sub(eq.rhs, dict)
    end

    newstates = []
    sts = states(sys)
    for j in eachindex(fullvars)
        if isirreducible(v_types, j)
            isdervar(s, j) || push!(newstates, fullvars[j])
        end
    end

    @set! sys.eqs = eqs
    @set! sys.states = newstates
    @set! sys.observed = [observed(sys); [lhs ~ rhs for (lhs, rhs) in pairs(subs)]]
    @set! sys.structure = nothing
    return sys
end

"""
    SparseMatrixCLIL{T, Ti}

The SparseMatrixCLIL represents a sparse matrix in two distinct ways:

1. As a sparse (in both row and column) n x m matrix
2. As a row-dense, column-sparse k x m matrix

The data structure keeps a permutation between the row order of the two representations.
Swapping the rows in one does not affect the other.

On construction, the second representation is equivalent to the first with fully-sparse
rows removed, though this may cease being true as row permutations are being applied
to the matrix.

The default structure of the `SparseMatrixCLIL` type is the second structure, while
the first is available via the thin `AsSubMatrix` wrapper.
"""
struct SparseMatrixCLIL{T, Ti<:Integer} <: AbstractSparseMatrix{T, Ti}
    nparentrows::Int
    ncols::Int
    nzrows::Vector{Ti}
    row_cols::Vector{Vector{Ti}}
    row_vals::Vector{Vector{T}}
end
Base.size(S::SparseMatrixCLIL) = (length(S.nzrows), S.ncols)
Base.copy(S::SparseMatrixCLIL{T, Ti}) where {T, Ti} =
    SparseMatrixCLIL(S.nparentrows, S.ncols, copy(S.nzrows), copy(S.row_cols), copy(S.row_vals))
function swaprows!(S::SparseMatrixCLIL, i, j)
    swap!(S.nzrows, i, j)
    swap!(S.row_cols, i, j)
    swap!(S.row_vals, i, j)
end

function bareiss_update_virtual_colswap_mtk!(zero!, M::SparseMatrixCLIL, k, swapto, pivot, last_pivot; pivot_equal_optimization=true)
    # for ei in nzrows(>= k)
    eadj = M.row_cols
    old_cadj = M.row_vals
    vpivot = swapto[2]

    ## N.B.: Micro-optimization
    #
    # For rows that do not have an entry in the eliminated column, all this
    # update does is multiply the row in question by `pivot/last_pivot` (the
    # result of which is guaranteed to be integer by general properties of the
    # bareiss algorithm, even if `pivot/last_pivot` is not).
    #
    # Thus, when `pivot == last pivot`, we can skip the update for any rows that
    # do not have an entry in the eliminated column (because we'd simply be
    # multiplying by 1).
    #
    # As an additional MTK-specific enhancement, we further allow the case
    # when the absolute values are equal, i.e. effectively multiplying the row
    # by `-1`. To ensure this is legal, we need to show two things.
    # 1. The multiplication does not change the answer and
    # 2. The multiplication does not affect the fraction-freeness of the Bareiss
    #    algorithm.
    #
    # For point 1, remember that we're working on a system of linear equations,
    # so it is always legal for us to multiply any row by a sclar without changing
    # the underlying system of equations.
    #
    # For point 2, note that the factorization we're now computing is the same
    # as if we had multiplied the corresponding row (accounting for row swaps)
    # in the original matrix by `last_pivot/pivot`, ensuring that the matrix
    # itself has integral entries when `last_pivot/pivot` is integral (here we
    # have -1, which counts). We could use the more general integrality
    # condition, but that would in turn disturb the growth bounds on the
    # factorization matrix entries that the bareiss algorithm guarantees. To be
    # conservative, we leave it at this, as this captures the most important
    # case for MTK (where most pivots are `1` or `-1`).
    pivot_equal = pivot_equal_optimization && abs(pivot) == abs(last_pivot)

    for ei in k+1:size(M, 1)
        # elimate `v`
        coeff = 0
        ivars = eadj[ei]
        vj = findfirst(isequal(vpivot), ivars)
        if vj !== nothing
            coeff = old_cadj[ei][vj]
            deleteat!(old_cadj[ei], vj)
            deleteat!(eadj[ei], vj)
        elseif pivot_equal
            continue
        end

        # the pivot row
        kvars = eadj[k]
        kcoeffs = old_cadj[k]
        # the elimination target
        ivars = eadj[ei]
        icoeffs = old_cadj[ei]

        tmp_incidence = similar(eadj[ei], 0)
        tmp_coeffs = similar(old_cadj[ei], 0)
        vars = union(ivars, kvars)

        for v in vars
            v == vpivot && continue
            ck = getcoeff(kvars, kcoeffs, v)
            ci = getcoeff(ivars, icoeffs, v)
            ci = (pivot*ci - coeff*ck) ÷ last_pivot
            if ci !== 0
                push!(tmp_incidence, v)
                push!(tmp_coeffs, ci)
            end
        end

        eadj[ei] = tmp_incidence
        old_cadj[ei] = tmp_coeffs
    end

    # Swap pivots to the front of the coefficient list
    # TODO: This prevents the coefficient list from being sorted, making
    # the rest of the algorithm much more expensive
    pivot_idx = findfirst(==(vpivot), eadj[k])
    deleteat!(eadj[k], pivot_idx)
    deleteat!(old_cadj[k], pivot_idx)
    pushfirst!(eadj[k], vpivot)
    pushfirst!(old_cadj[k], pivot)
end

function bareiss_update_virtual_colswap_mtk!(zero!, M::AbstractMatrix, k, swapto, pivot, last_pivot; pivot_equal_optimization=true)
    if pivot_equal_optimization
        error("MTK pivot micro-optimization not implemented for `$(typeof(M))`.
            Turn off the optimization for debugging or use a different matrix type.")
    end
    bareiss_update_virtual_colswap!(zero!, M, k, swapto, pivot, last_pivot)
end

struct AsSubMatrix{T, Ti<:Integer} <: AbstractSparseMatrix{T, Ti}
    M::SparseMatrixCLIL{T, Ti}
end
Base.size(S::AsSubMatrix) = (S.M.nparentrows, S.M.ncols)

function Base.getindex(S::SparseMatrixCLIL{T}, i1, i2) where {T}
    checkbounds(S, i1, i2)

    nncol = findfirst(==(i2), S.row_cols[i1])
    isnothing(nncol) && return zero(T)

    return S.row_vals[i1][nncol]
end

function Base.getindex(S::AsSubMatrix{T}, i1, i2) where {T}
    checkbounds(S, i1, i2)
    S = S.M

    nnrow = findfirst(==(i1), S.nzrows)
    isnothing(nnrow) && return zero(T)

    nncol = findfirst(==(i2), S.row_cols[nnrow])
    isnothing(nncol) && return zero(T)

    return S.row_vals[nnrow][nncol]
end

"""
$(SIGNATURES)

Find the first linear variable such that `𝑠neighbors(adj, i)[j]` is true given
the `constraint`.
"""
@inline function find_first_linear_variable(
        M::SparseMatrixCLIL,
        range,
        mask,
        constraint,
    )
    eadj = M.row_cols
    for i in range
        vertices = eadj[i]
        if constraint(length(vertices))
            for (j, v) in enumerate(vertices)
                (mask === nothing || mask[v]) && return (CartesianIndex(i, v), M.row_vals[i][j])
            end
        end
    end
    return nothing
end

@inline function find_first_linear_variable(
        M::AbstractMatrix,
        range,
        mask,
        constraint,
    )
    for i in range
        row = @view M[i, :]
        if constraint(count(!iszero, row))
            for (v, val) in enumerate(row)
                iszero(val) && continue
                if mask === nothing || mask[v]
                    return CartesianIndex(i, v), val
                end
            end
        end
    end
    return nothing
end

function find_masked_pivot(variables, M, k)
    r = find_first_linear_variable(M, k:size(M,1), variables, isequal(1))
    r !== nothing && return r
    r = find_first_linear_variable(M, k:size(M,1), variables, isequal(2))
    r !== nothing && return r
    r = find_first_linear_variable(M, k:size(M,1), variables, _->true)
    return r
end

function alias_eliminate_graph(s::SystemStructure, is_linear_equations, eadj, cadj)
    @unpack graph, varassoc = s
    invvarassoc = inverse_mapping(varassoc)

    old_cadj = map(copy, cadj)

    is_not_potential_state = iszero.(varassoc)
    is_linear_variables = copy(is_not_potential_state)
    for i in 𝑠vertices(graph); is_linear_equations[i] && continue
        for j in 𝑠neighbors(graph, i)
            is_linear_variables[j] = false
        end
    end
    solvable_variables = findall(is_linear_variables)

    linear_equations = findall(is_linear_equations)

    mm = SparseMatrixCLIL(nsrcs(graph),
                          ndsts(graph),
                          linear_equations, eadj, old_cadj)

    function do_bareiss!(M, cadj=nothing)
        rank1 = rank2 = nothing
        function find_pivot(M, k)
            if rank1 === nothing
                r = find_masked_pivot(is_linear_variables, M, k)
                r !== nothing && return r
                rank1 = k - 1
            end
            if rank2 === nothing
                r = find_masked_pivot(is_not_potential_state, M, k)
                r !== nothing && return r
                rank2 = k - 1
            end
            return find_masked_pivot(nothing, M, k)
        end
        function myswaprows!(M, i, j)
            cadj !== nothing && swap!(cadj, i, j)
            swaprows!(M, i, j)
        end
        bareiss_ops = ((M,i,j)->nothing, myswaprows!, bareiss_update_virtual_colswap_mtk!, bareiss_zero!)
        rank3 = bareiss!(M, bareiss_ops; find_pivot)
        rank1 = something(rank1, rank3)
        rank2 = something(rank2, rank3)
        (rank1, rank2, rank3)
    end

    # mm2 = Array(copy(mm))
    # @show do_bareiss!(mm2)
    # display(mm2)

    (rank1, rank2, rank3) = do_bareiss!(mm, cadj)

    v_solved = [eadj[i][1] for i in 1:rank1]
    v_eliminated = setdiff(solvable_variables, v_solved)
    n_null_vars = length(v_eliminated)

    v_types = fill(KEEP, ndsts(graph))
    for v in v_eliminated
        v_types[v] = 0
    end

    # kind of like the backward substitution
    for ei in reverse(1:rank2)
        locally_structure_simplify!(
                                    (eadj[ei], old_cadj[ei]),
                                    invvarassoc, v_eliminated, v_types
                                   )
    end

    reduced = false
    for ei in 1:rank2
        if length(cadj[ei]) >= length(old_cadj[ei])
            cadj[ei] = old_cadj[ei]
        else
            # MEMORY ALIAS of a vector
            eadj[ei] = 𝑠neighbors(graph, linear_equations[ei])
            reduced |= locally_structure_simplify!(
                                                   (eadj[ei], cadj[ei]),
                                                   invvarassoc, v_eliminated, v_types
                                                  )
        end
    end

    while reduced
        for ei in 1:rank2
            if !isempty(eadj[ei])
                reduced |= locally_structure_simplify!(
                                                       (eadj[ei], cadj[ei]),
                                                       invvarassoc, v_eliminated, v_types
                                                      )
                reduced && break # go back to the begining of equations
            end
        end
    end

    for ei in rank2+1:length(linear_equations)
        cadj[ei] = old_cadj[ei]
    end

    for (ei, e) in enumerate(linear_equations)
        graph.fadjlist[e] = eadj[ei]
    end

    degenerate_equations = rank3 < length(linear_equations) ? linear_equations[rank3+1:end] : Int[]
    return v_eliminated, v_types, n_null_vars, degenerate_equations, linear_equations
end

iszeroterm(v_types, v) = v_types[v] == 0
isirreducible(v_types, v) = v_types[v] == KEEP
isalias(v_types, v) = v_types[v] > 0 && !isirreducible(v_types, v)
alias(v_types, v) = v_types[v]
negalias(v_types, v) = -v_types[v]

function locally_structure_simplify!(
        (vars, coeffs),
        invvarassoc, v_eliminated, v_types
       )
    while length(vars) > 1 && any(!isequal(KEEP), (v_types[v] for v in @view vars[2:end]))
        for vj in 2:length(vars)
            v = vars[vj]
            if isirreducible(v_types, v)
                continue
            elseif iszeroterm(v_types, v)
                deleteat!(vars, vj)
                deleteat!(coeffs, vj)
                break
            else
                coeff = coeffs[vj]
                if isalias(v_types, v)
                    v = alias(v_types, v)
                else
                    v = negalias(v_types, v)
                    coeff = -coeff
                end

                has_v = false
                for vi in 2:length(vars)
                    (vi !== vj && vars[vi] == v) || continue
                    has_v = true
                    c = (coeffs[vi] += coeff)
                    if c == 0
                        if vi < vj
                            deleteat!(vars, [vi, vj])
                            deleteat!(coeffs, [vi, vj])
                        else
                            deleteat!(vars, [vj, vi])
                            deleteat!(coeffs, [vj, vi])
                        end
                    end
                    break
                end # for vi

                if has_v
                    break
                else
                    vars[vj] = v
                    coeffs[vj] = coeff
                end # if
            end # else
        end # for
    end # while

    v = first(vars)

    # Do not attempt to eliminate derivatives
    invvarassoc[v] != 0 && return false

    if length(vars) == 1
        push!(v_eliminated, v)
        v_types[v] = 0
        empty!(vars); empty!(coeffs)
        return true
    elseif length(vars) == 2 && abs(coeffs[1]) == abs(coeffs[2])
        if (coeffs[1] > 0 && coeffs[2] < 0) || (coeffs[1] < 0 && coeffs[2] > 0)
            # positive alias
            push!(v_eliminated, v)
            v_types[v] = vars[2]
        else
            # negative alias
            push!(v_eliminated, v)
            v_types[v] = -vars[2]
        end
        empty!(vars); empty!(coeffs)
        return true
    end
    return false
end

swap!(v, i, j) = v[i], v[j] = v[j], v[i]

function getcoeff(vars, coeffs, var)
    for (vj, v) in enumerate(vars)
        v == var && return coeffs[vj]
    end
    return 0
end

function inverse_mapping(assoc)
    invassoc = zeros(Int, length(assoc))
    for (i, v) in enumerate(assoc)
        v <= 0 && continue
        invassoc[v] = i
    end
    return invassoc
end

"""
$(SIGNATURES)

Use Kahn's algorithm to topologically sort observed equations.

Example:
```julia
julia> @variables t x(t) y(t) z(t) k(t)
(t, x(t), y(t), z(t), k(t))

julia> eqs = [
           x ~ y + z
           z ~ 2
           y ~ 2z + k
       ];

julia> ModelingToolkit.topsort_equations(eqs, [x, y, z, k])
3-element Vector{Equation}:
 Equation(z(t), 2)
 Equation(y(t), k(t) + 2z(t))
 Equation(x(t), y(t) + z(t))
```
"""
function topsort_equations(eqs, states; check=true)
    graph, assigns = observed2graph(eqs, states)
    neqs = length(eqs)
    degrees = zeros(Int, neqs)

    for 𝑠eq in 1:length(eqs); var = assigns[𝑠eq]
        for 𝑑eq in 𝑑neighbors(graph, var)
            # 𝑠eq => 𝑑eq
            degrees[𝑑eq] += 1
        end
    end

    q = Queue{Int}(neqs)
    for (i, d) in enumerate(degrees)
        d == 0 && enqueue!(q, i)
    end

    idx = 0
    ordered_eqs = similar(eqs, 0); sizehint!(ordered_eqs, neqs)
    while !isempty(q)
        𝑠eq = dequeue!(q)
        idx+=1
        push!(ordered_eqs, eqs[𝑠eq])
        var = assigns[𝑠eq]
        for 𝑑eq in 𝑑neighbors(graph, var)
            degree = degrees[𝑑eq] = degrees[𝑑eq] - 1
            degree == 0 && enqueue!(q, 𝑑eq)
        end
    end

    (check && idx != neqs) && throw(ArgumentError("The equations have at least one cycle."))

    return ordered_eqs
end

function observed2graph(eqs, states)
    graph = BipartiteGraph(length(eqs), length(states))
    v2j = Dict(states .=> 1:length(states))

    # `assigns: eq -> var`, `eq` defines `var`
    assigns = similar(eqs, Int)

    for (i, eq) in enumerate(eqs)
        lhs_j = get(v2j, eq.lhs, nothing)
        lhs_j === nothing && throw(ArgumentError("The lhs $(eq.lhs) of $eq, doesn't appear in states."))
        assigns[i] = lhs_j
        vs = vars(eq.rhs)
        for v in vs
            j = get(v2j, v, nothing)
            j !== nothing && add_edge!(graph, i, j)
        end
    end

    return graph, assigns
end

function fixpoint_sub(x, dict)
    y = substitute(x, dict)
    while !isequal(x, y)
        y = x
        x = substitute(y, dict)
    end

    return x
end

function substitute_aliases(eqs, dict)
    sub = Base.Fix2(fixpoint_sub, dict)
    map(eq->eq.lhs ~ sub(eq.rhs), eqs)
end
