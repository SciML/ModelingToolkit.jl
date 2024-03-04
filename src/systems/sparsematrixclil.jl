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
struct SparseMatrixCLIL{T, Ti <: Integer} <: AbstractSparseMatrix{T, Ti}
    nparentrows::Int
    ncols::Int
    nzrows::Vector{Ti}
    row_cols::Vector{Vector{Ti}} # issorted
    row_vals::Vector{Vector{T}}
end
Base.size(S::SparseMatrixCLIL) = (length(S.nzrows), S.ncols)
function Base.copy(S::SparseMatrixCLIL{T, Ti}) where {T, Ti}
    SparseMatrixCLIL(S.nparentrows, S.ncols, copy(S.nzrows), map(copy, S.row_cols),
        map(copy, S.row_vals))
end
function swaprows!(S::SparseMatrixCLIL, i, j)
    i == j && return
    swap!(S.nzrows, i, j)
    swap!(S.row_cols, i, j)
    swap!(S.row_vals, i, j)
end

function Base.convert(::Type{SparseMatrixCLIL{T, Ti}}, S::SparseMatrixCLIL) where {T, Ti}
    return SparseMatrixCLIL(S.nparentrows,
        S.ncols,
        copy.(S.nzrows),
        copy.(S.row_cols),
        [T.(row) for row in S.row_vals])
end

function SparseMatrixCLIL(mm::AbstractMatrix)
    nrows, ncols = size(mm)
    row_cols = [findall(!iszero, row) for row in eachrow(mm)]
    row_vals = [row[cols] for (row, cols) in zip(eachrow(mm), row_cols)]
    SparseMatrixCLIL(nrows, ncols, Int[1:length(row_cols);], row_cols, row_vals)
end

struct CLILVector{T, Ti} <: AbstractSparseVector{T, Ti}
    vec::SparseVector{T, Ti}
end
Base.hash(v::CLILVector, s::UInt) = hash(v.vec, s) ⊻ 0xc71be0e9ccb75fbd
Base.size(v::CLILVector) = Base.size(v.vec)
Base.getindex(v::CLILVector, idx::Integer...) = Base.getindex(v.vec, idx...)
Base.setindex!(vec::CLILVector, v, idx::Integer...) = Base.setindex!(vec.vec, v, idx...)
function Base.view(a::SparseMatrixCLIL, i::Integer, ::Colon)
    CLILVector(SparseVector(a.ncols, a.row_cols[i], a.row_vals[i]))
end
SparseArrays.nonzeroinds(a::CLILVector) = SparseArrays.nonzeroinds(a.vec)
SparseArrays.nonzeros(a::CLILVector) = SparseArrays.nonzeros(a.vec)
SparseArrays.nnz(a::CLILVector) = nnz(a.vec)

function Base.setindex!(S::SparseMatrixCLIL, v::CLILVector, i::Integer, c::Colon)
    if v.vec.n != S.ncols
        throw(BoundsError(v, 1:(S.ncols)))
    end
    any(iszero, v.vec.nzval) && error("setindex failed")
    S.row_cols[i] = copy(v.vec.nzind)
    S.row_vals[i] = copy(v.vec.nzval)
    return v
end

zero!(a::AbstractArray{T}) where {T} = a[:] .= zero(T)
zero!(a::SparseVector) = (empty!(a.nzind); empty!(a.nzval))
zero!(a::CLILVector) = zero!(a.vec)
SparseArrays.dropzeros!(a::CLILVector) = SparseArrays.dropzeros!(a.vec)

struct NonZeros{T <: AbstractArray}
    v::T
end
Base.pairs(nz::NonZeros{<:CLILVector}) = NonZerosPairs(nz.v)

struct NonZerosPairs{T <: AbstractArray}
    v::T
end

Base.IteratorSize(::Type{<:NonZerosPairs}) = Base.SizeUnknown()
# N.B.: Because of how we're using this, this must be robust to modification of
# the underlying vector. As such, we treat this as an iteration over indices
# that happens to short cut using the sparse structure and sortedness of the
# array.
function Base.iterate(nzp::NonZerosPairs{<:CLILVector}, (idx, col))
    v = nzp.v.vec
    nzind = v.nzind
    nzval = v.nzval
    if idx > length(nzind)
        idx = length(col)
    end
    oldcol = nzind[idx]
    if col != oldcol
        # The vector was changed since the last iteration. Find our
        # place in the vector again.
        tail = col > oldcol ? (@view nzind[(idx + 1):end]) : (@view nzind[1:idx])
        tail_i = searchsortedfirst(tail, col + 1)
        # No remaining indices.
        tail_i > length(tail) && return nothing
        new_idx = col > oldcol ? idx + tail_i : tail_i
        new_col = nzind[new_idx]
        return (new_col => nzval[new_idx], (new_idx, new_col))
    end
    idx == length(nzind) && return nothing
    new_col = nzind[idx + 1]
    return (new_col => nzval[idx + 1], (idx + 1, new_col))
end

function Base.iterate(nzp::NonZerosPairs{<:CLILVector})
    v = nzp.v.vec
    nzind = v.nzind
    nzval = v.nzval
    isempty(nzind) && return nothing
    return nzind[1] => nzval[1], (1, nzind[1])
end

# Arguably this is how nonzeros should behave in the first place, but let's
# build something that works for us here and worry about it later.
nonzerosmap(a::CLILVector) = NonZeros(a)

using FindFirstFunctions: findfirstequal

function bareiss_update_virtual_colswap_mtk!(zero!, M::SparseMatrixCLIL, k, swapto, pivot,
        last_pivot; pivot_equal_optimization = true)
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
    # so it is always legal for us to multiply any row by a scalar without changing
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
    @inbounds for ei in (k + 1):size(M, 1)
        # eliminate `v`
        coeff = 0
        ivars = eadj[ei]
        vj = findfirstequal(vpivot, ivars)
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

        numkvars = length(kvars)
        numivars = length(ivars)
        tmp_incidence = similar(eadj[ei], numkvars + numivars)
        tmp_coeffs = similar(old_cadj[ei], numkvars + numivars)
        tmp_len = 0
        kvind = ivind = 0
        if _debug_mode
            # in debug mode, we at least check to confirm we're iterating over
            # `v`s in the correct order
            vars = sort(union(ivars, kvars))
            vi = 0
        end
        if numivars > 0 && numkvars > 0
            kvv = kvars[kvind += 1]
            ivv = ivars[ivind += 1]
            dobreak = false
            while true
                if kvv == ivv
                    v = kvv
                    ck = kcoeffs[kvind]
                    ci = icoeffs[ivind]
                    kvind += 1
                    ivind += 1
                    if kvind > numkvars
                        dobreak = true
                    else
                        kvv = kvars[kvind]
                    end
                    if ivind > numivars
                        dobreak = true
                    else
                        ivv = ivars[ivind]
                    end
                    p1 = Base.Checked.checked_mul(pivot, ci)
                    p2 = Base.Checked.checked_mul(coeff, ck)
                    ci = exactdiv(Base.Checked.checked_sub(p1, p2), last_pivot)
                elseif kvv < ivv
                    v = kvv
                    ck = kcoeffs[kvind]
                    kvind += 1
                    if kvind > numkvars
                        dobreak = true
                    else
                        kvv = kvars[kvind]
                    end
                    p2 = Base.Checked.checked_mul(coeff, ck)
                    ci = exactdiv(Base.Checked.checked_neg(p2), last_pivot)
                else # kvv > ivv
                    v = ivv
                    ci = icoeffs[ivind]
                    ivind += 1
                    if ivind > numivars
                        dobreak = true
                    else
                        ivv = ivars[ivind]
                    end
                    ci = exactdiv(Base.Checked.checked_mul(pivot, ci), last_pivot)
                end
                if _debug_mode
                    @assert v == vars[vi += 1]
                end
                if v != vpivot && !iszero(ci)
                    tmp_incidence[tmp_len += 1] = v
                    tmp_coeffs[tmp_len] = ci
                end
                dobreak && break
            end
        elseif numkvars > 0
            ivind = 1
            kvv = kvars[kvind += 1]
        elseif numivars > 0
            kvind = 1
            ivv = ivars[ivind += 1]
        end
        if kvind <= numkvars
            v = kvv
            while true
                if _debug_mode
                    @assert v == vars[vi += 1]
                end
                if v != vpivot
                    ck = kcoeffs[kvind]
                    p2 = Base.Checked.checked_mul(coeff, ck)
                    ci = exactdiv(Base.Checked.checked_neg(p2), last_pivot)
                    if !iszero(ci)
                        tmp_incidence[tmp_len += 1] = v
                        tmp_coeffs[tmp_len] = ci
                    end
                end
                (kvind == numkvars) && break
                v = kvars[kvind += 1]
            end
        elseif ivind <= numivars
            v = ivv
            while true
                if _debug_mode
                    @assert v == vars[vi += 1]
                end
                if v != vpivot
                    p1 = Base.Checked.checked_mul(pivot, icoeffs[ivind])
                    ci = exactdiv(p1, last_pivot)
                    if !iszero(ci)
                        tmp_incidence[tmp_len += 1] = v
                        tmp_coeffs[tmp_len] = ci
                    end
                end
                (ivind == numivars) && break
                v = ivars[ivind += 1]
            end
        end
        resize!(tmp_incidence, tmp_len)
        resize!(tmp_coeffs, tmp_len)
        eadj[ei] = tmp_incidence
        old_cadj[ei] = tmp_coeffs
    end
end

function bareiss_update_virtual_colswap_mtk!(zero!, M::AbstractMatrix, k, swapto, pivot,
        last_pivot; pivot_equal_optimization = true)
    if pivot_equal_optimization
        error("MTK pivot micro-optimization not implemented for `$(typeof(M))`.
            Turn off the optimization for debugging or use a different matrix type.")
    end
    bareiss_update_virtual_colswap!(zero!, M, k, swapto, pivot, last_pivot)
end

struct AsSubMatrix{T, Ti <: Integer} <: AbstractSparseMatrix{T, Ti}
    M::SparseMatrixCLIL{T, Ti}
end
Base.size(S::AsSubMatrix) = (S.M.nparentrows, S.M.ncols)

function Base.getindex(S::SparseMatrixCLIL{T}, i1::Integer, i2::Integer) where {T}
    checkbounds(S, i1, i2)

    col = S.row_cols[i1]
    nncol = searchsortedfirst(col, i2)
    (nncol > length(col) || col[nncol] != i2) && return zero(T)

    return S.row_vals[i1][nncol]
end

function Base.getindex(S::AsSubMatrix{T}, i1::Integer, i2::Integer) where {T}
    checkbounds(S, i1, i2)
    S = S.M

    nnrow = findfirst(==(i1), S.nzrows)
    isnothing(nnrow) && return zero(T)

    col = S.row_cols[nnrow]
    nncol = searchsortedfirst(col, i2)
    (nncol > length(col) || col[nncol] != i2) && return zero(T)

    return S.row_vals[nnrow][nncol]
end
