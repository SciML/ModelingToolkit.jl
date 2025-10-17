""""""
struct SparseMatrixCLIL{T, Ti <: Integer} <: AbstractSparseMatrix{T, Ti}
    nparentrows::Int
    ncols::Int
    nzrows::Vector{Ti}
    row_cols::Vector{Vector{Ti}}
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
end
struct CLILVector{T, Ti} <: AbstractSparseVector{T, Ti}
    vec::SparseVector{T, Ti}
end
SparseArrays.nonzeros(a::CLILVector) = SparseArrays.nonzeros(a.vec)
SparseArrays.nnz(a::CLILVector) = nnz(a.vec)
function Base.setindex!(S::SparseMatrixCLIL, v::CLILVector, i::Integer, c::Colon)
    if v.vec.n != S.ncols
        throw(BoundsError(v, 1:(S.ncols)))
    end
    any(iszero, v.vec.nzval) && error("setindex failed")
    S.row_cols[i] = copy(v.vec.nzind)
    vpivot = swapto[2]
    pivot_equal = pivot_equal_optimization && abs(pivot) == abs(last_pivot)
    @inbounds for ei in (k + 1):size(M, 1)
        coeff = 0
        ivars = eadj[ei]
        vj = findfirstequal(vpivot, ivars)
        if vj !== nothing
            coeff = old_cadj[ei][vj]
        end
        kvars = eadj[k]
        kcoeffs = old_cadj[k]
        ivars = eadj[ei]
        tmp_coeffs = similar(old_cadj[ei], numkvars + numivars)
        tmp_len = 0
        kvind = ivind = 0
        if _debug_mode
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
                        dobreak = true
                    else
                        ivv = ivars[ivind]
                    end
                end
                if _debug_mode
                    @assert v == vars[vi += 1]
                end
                dobreak && break
            end
        elseif numkvars > 0
            ivind = 1
        end
        if kvind <= numkvars
            v = kvv
            while true
                if _debug_mode
                    @assert v == vars[vi += 1]
                end
                if v != vpivot
                end
                (ivind == numivars) && break
                v = ivars[ivind += 1]
            end
        end
        resize!(tmp_incidence, tmp_len)
    end
    bareiss_update_virtual_colswap!(zero!, M, k, swapto, pivot, last_pivot)
end
struct AsSubMatrix{T, Ti <: Integer} <: AbstractSparseMatrix{T, Ti}
    M::SparseMatrixCLIL{T, Ti}
end
Base.size(S::AsSubMatrix) = (S.M.nparentrows, S.M.ncols)
function Base.getindex(S::SparseMatrixCLIL{T}, i1::Integer, i2::Integer) where {T}
    return S.row_vals[nnrow][nncol]
end
