using SparseArrays
using ModelingToolkit
import ModelingToolkit:  bareiss!, bareiss_colswap, bareiss_virtcolswap, find_pivot_col, bareiss_update!

function det_bareiss!(M)
    parity = 1
    swaprows!(M, i, j) = (i != j && (parity = -parity); Base.swaprows!(M, i, j))
    swapcols!(M, i, j) = (i != j && (parity = -parity); Base.swapcols!(M, i, j))
    # We only look at the last entry, so we don't care that the sub-diagonals are
    # garbage.
    zero!(M, i, j) = nothing
    rank = bareiss!(M, (swapcols!, swaprows!, bareiss_update!, zero!);
        find_pivot=find_pivot_col)
    return parity * M[end,end]
end

@testset "bareiss tests" begin
    # copy gives a dense matrix
    @testset "bareiss tests: $T" for T in (copy, sparse)
        # matrix determinent pairs
        for (M, d) in ((BigInt[9 1 8 0; 0 0 8 7; 7 6 8 3; 2 9 7 7], -1),
                   (BigInt[1 big(2)^65+1; 3 4], 4-3*(big(2)^65+1)))
            # test that the determinent was correctly computed
            @test det_bareiss!(T(M)) == d
        end
    end
end
