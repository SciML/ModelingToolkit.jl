# Keeps compatibility with bariess code movoed to Base/stdlib on older releases

using LinearAlgebra
using SparseArrays
using SparseArrays: AbstractSparseMatrixCSC

macro swap(a, b)
    esc(:(($a, $b) = ($b, $a)))
end

if isdefined(Base, :swaprows!)
    import Base: swaprows!
else
    function swaprows!(a::AbstractMatrix, i, j)
        i == j && return
        rows = axes(a,1)
        @boundscheck i in rows || throw(BoundsError(a, (:,i)))
        @boundscheck j in rows || throw(BoundsError(a, (:,j)))
        for k in axes(a,2)
            @inbounds a[i,k],a[j,k] = a[j,k],a[i,k]
        end
    end
    function Base.circshift!(a::AbstractVector, shift::Integer)
        n = length(a)
        n == 0 && return
        shift = mod(shift, n)
        shift == 0 && return
        reverse!(a, 1, shift)
        reverse!(a, shift+1, length(a))
        reverse!(a)
        return a
    end
    function Base.swapcols!(A::AbstractSparseMatrixCSC, i, j)
        i == j && return

        # For simplicitly, let i denote the smaller of the two columns
        j < i && @swap(i, j)

        colptr = getcolptr(A)
        irow = colptr[i]:(colptr[i+1]-1)
        jrow = colptr[j]:(colptr[j+1]-1)

        function rangeexchange!(arr, irow, jrow)
            if length(irow) == length(jrow)
                for (a, b) in zip(irow, jrow)
                    @inbounds @swap(arr[i], arr[j])
                end
                return
            end
            # This is similar to the triple-reverse tricks for
            # circshift!, except that we have three ranges here,
            # so it ends up being 4 reverse calls (but still
            # 2 overall reversals for the memory range). Like
            # circshift!, there's also a cycle chasing algorithm
            # with optimal memory complexity, but the performance
            # tradeoffs against this implementation are non-trivial,
            # so let's just do this simple thing for now.
            # See https://github.com/JuliaLang/julia/pull/42676 for
            # discussion of circshift!-like algorithms.
            reverse!(@view arr[irow])
            reverse!(@view arr[jrow])
            reverse!(@view arr[(last(irow)+1):(first(jrow)-1)])
            reverse!(@view arr[first(irow):last(jrow)])
        end
        rangeexchange!(rowvals(A), irow, jrow)
        rangeexchange!(nonzeros(A), irow, jrow)

        if length(irow) != length(jrow)
            @inbounds colptr[i+1:j] .+= length(jrow) - length(irow)
        end
        return nothing
    end
    function swaprows!(A::AbstractSparseMatrixCSC, i, j)
        # For simplicitly, let i denote the smaller of the two rows
        j < i && @swap(i, j)

        rows = rowvals(A)
        vals = nonzeros(A)
        for col = 1:size(A, 2)
            rr = nzrange(A, col)
            iidx = searchsortedfirst(@view(rows[rr]), i)
            has_i = iidx <= length(rr) && rows[rr[iidx]] == i

            jrange = has_i ? (iidx:last(rr)) : rr
            jidx = searchsortedlast(@view(rows[jrange]), j)
            has_j = jidx != 0 && rows[jrange[jidx]] == j

            if !has_j && !has_i
                # Has neither row - nothing to do
                continue
            elseif has_i && has_j
                # This column had both i and j rows - swap them
                @swap(vals[rr[iidx]], vals[jrange[jidx]])
            elseif has_i
                # Update the rowval and then rotate both nonzeros
                # and the remaining rowvals into the correct place
                rows[rr[iidx]] = j
                jidx == 0 && continue
                rotate_range = rr[iidx]:jrange[jidx]
                circshift!(@view(vals[rotate_range]), -1)
                circshift!(@view(rows[rotate_range]), -1)
            else
                # Same as i, but in the opposite direction
                @assert has_j
                rows[jrange[jidx]] = i
                iidx > length(rr) && continue
                rotate_range = rr[iidx]:jrange[jidx]
                circshift!(@view(vals[rotate_range]), 1)
                circshift!(@view(rows[rotate_range]), 1)
            end
        end
        return nothing
    end
end

if isdefined(LinearAlgebra, :bareiss!)
    import LinearAlgebra: bareiss!, bareiss_update_virtual_colswap!, bareiss_zero!
else
    function bareiss_update!(zero!, M::Matrix, k, swapto, pivot, prev_pivot)
        for i in k+1:size(M, 2), j in k+1:size(M, 1)
            M[j,i] = exactdiv(M[j,i]*pivot - M[j,k]*M[k,i], prev_pivot)
        end
        zero!(M, k+1:size(M, 1), k)
    end

    function bareiss_update!(zero!, M::AbstractMatrix, k, swapto, pivot, prev_pivot)
        V = @view M[k+1:end, k+1:end]
        V .= exactdiv.(V * pivot - M[k+1:end, k] * M[k, k+1:end]', prev_pivot)
        zero!(M, k+1:size(M, 1), k)
    end

    function bareiss_update_virtual_colswap!(zero!, M::AbstractMatrix, k, swapto, pivot, prev_pivot)
        V = @view M[k+1:end, :]
        V .= exactdiv.(V * pivot - M[k+1:end, swapto[2]] * M[k, :]', prev_pivot)
        zero!(M, k+1:size(M, 1), swapto[2])
    end

    bareiss_zero!(M, i, j) = M[i,j] .= zero(eltype(M))

    function find_pivot_col(M, i)
        p = findfirst(!iszero, @view M[i,i:end])
        p === nothing && return nothing
        idx = CartesianIndex(i, p + i - 1)
        (idx, M[idx])
    end

    function find_pivot_any(M, i)
        p = findfirst(!iszero, @view M[i:end,i:end])
        p === nothing && return nothing
        idx = p + CartesianIndex(i - 1, i - 1)
        (idx, M[idx])
    end

    const bareiss_colswap = (Base.swapcols!, swaprows!, bareiss_update!, bareiss_zero!)
    const bareiss_virtcolswap = ((M,i,j)->nothing, swaprows!, bareiss_update_virtual_colswap!, bareiss_zero!)

    """
        bareiss!(M)

    Perform Bareiss's fraction-free row-reduction algorithm on the matrix `M`.
    Optionally, a specific pivoting method may be specified.
    """
    function bareiss!(M::AbstractMatrix,
                         (swapcols!, swaprows!, update!, zero!) = bareiss_colswap;
                      find_pivot=find_pivot_any)
        prev = one(eltype(M))
        n = size(M, 1)
        for k in 1:n
            r = find_pivot(M, k)
            r === nothing && return k - 1
            (swapto, pivot) = r
            if CartesianIndex(k, k) != swapto
                swapcols!(M, k, swapto[2])
                swaprows!(M, k, swapto[1])
            end
            update!(zero!, M, k, swapto, pivot, prev)
            prev = pivot
        end
        return n
    end
end
