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
    row_cols::Vector{Vector{Ti}} # issorted
    row_vals::Vector{Vector{T}}
end
Base.size(S::SparseMatrixCLIL) = (length(S.nzrows), S.ncols)
Base.copy(S::SparseMatrixCLIL{T, Ti}) where {T, Ti} =
    SparseMatrixCLIL(S.nparentrows, S.ncols, copy(S.nzrows), map(copy, S.row_cols), map(copy, S.row_vals))
function swaprows!(S::SparseMatrixCLIL, i, j)
    swap!(S.nzrows, i, j)
    swap!(S.row_cols, i, j)
    swap!(S.row_vals, i, j)
end

struct CLILVector{T, Ti} <: AbstractSparseVector{T, Ti}
    vec::SparseVector{T, Ti}
end
Base.size(v::CLILVector) = Base.size(v.vec)
Base.getindex(v::CLILVector, idx::Integer...) = Base.getindex(v.vec, idx...)
Base.setindex!(vec::CLILVector, v, idx::Integer...) = Base.setindex!(vec.vec, v, idx...)
Base.view(a::SparseMatrixCLIL, i::Integer, ::Colon) =
    CLILVector(SparseVector(a.ncols, a.row_cols[i], a.row_vals[i]))
SparseArrays.nonzeroinds(a::CLILVector) = SparseArrays.nonzeroinds(a.vec)
SparseArrays.nonzeros(a::CLILVector) = SparseArrays.nonzeros(a.vec)

function Base.setindex!(S::SparseMatrixCLIL, v::CLILVector, i::Integer, c::Colon)
    if v.vec.n != S.ncols
        throw(BoundsError(v, 1:S.ncols))
    end
    S.row_cols[i] = copy(v.vec.nzind)
    S.row_vals[i] = copy(v.vec.nzval)
    return v
end

zero!(a::AbstractArray{T}) where {T} = a[:] .= zero(T)
zero!(a::SparseVector) = (empty!(a.nzind); empty!(a.nzval))
zero!(a::CLILVector) = zero!(a.vec)

struct NonZeros{T <: AbstractArray}
    v::T
end
Base.pairs(nz::NonZeros{<:CLILVector}) = NonZerosPairs(nz.v)

struct NonZerosPairs{T <: AbstractArray}
    v::T
end

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
    if col !== oldcol
        # The vector was changed since the last iteration. Find our
        # place in the vector again.
        tail = col > oldcol ? (@view nzind[idx+1:end]) : (@view nzind[1:idx])
        tail_i = searchsortedfirst(tail, col + 1)
        # No remaining indices.
        tail_i > length(tail) && return nothing
        new_idx = col > oldcol ? idx + tail_i : tail_i
        new_col = nzind[new_idx]
        return (new_col=>nzval[new_idx], (new_idx, new_col))
    end
    idx == length(nzind) && return nothing
    new_col = nzind[idx+1]
    return (new_col=>nzval[idx+1], (idx+1, new_col))
end

function Base.iterate(nzp::NonZerosPairs{<:CLILVector})
    v = nzp.v.vec
    nzind = v.nzind
    nzval = v.nzval
    isempty(nzind) && return nothing
    return nzind[1]=>nzval[1], (1, nzind[1])
end

# Arguably this is how nonzeros should behave in the first place, but let's
# build something that works for us here and worry about it later.
nonzerosmap(a::CLILVector) = NonZeros(a)

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

    col = S.row_cols[i1]
    nncol = searchsortedfirst(col, i2)
    (nncol > length(col) || col[nncol] != i2) && return zero(T)

    return S.row_vals[i1][nncol]
end

function Base.getindex(S::AsSubMatrix{T}, i1, i2) where {T}
    checkbounds(S, i1, i2)
    S = S.M

    nnrow = findfirst(==(i1), S.nzrows)
    isnothing(nnrow) && return zero(T)

    col = S.row_cols[nnrow]
    nncol = searchsortedfirst(col, i2)
    (nncol > length(col) || col[nncol] != i2) && return zero(T)

    return S.row_vals[nnrow][nncol]
end
