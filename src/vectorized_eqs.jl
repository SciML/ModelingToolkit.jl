struct HybridEquationView <: AbstractVector{Equation}
    vectorized::Vector{Any}
    scalarized::Vector{Equation}
    lens::Vector{Int}
end

function HybridEquationView(vectorized, scalarized)
    n = length(vectorized)
    lens = Vector{Int}(undef, n + 2)
    lens[1] = 1
    for i in 1:n
        lens[i + 1] = lens[i] + length(vectorized[i])
    end
    lens[n + 2] = lens[n + 1] + length(scalarized)
    HybridEquationView(vectorized, scalarized, lens)
end

function Base.size(h::HybridEquationView)
    (h.lens[end] - 1,)
end

function Base.getindex(h::HybridEquationView, i::Int)
    @boundscheck 1 <= i <= length(h) || throw(BoundsError(h, i))
    lens = h.lens
    m = length(lens)
    j = searchsortedfirst(lens, i)
    lj = lens[j] == i ? j : j - 1
    k = lens[lj]
    l = i - k + 1
    return lj + 1 == m ? h.scalarized[l] : scalarize(h.vectorized[k][l])
end
