struct BiGraph{T}
    data::Vector{Vector{T}}
end

function sys2bigraph(sys)
    ss = states(sys)
    data = Operation[]
    for eq in sys.eqs
        es = []
        lhs = eq.lhs
        lhs.op isa Differential && push!(eq, lhs)
        push!(data, es)
    end
end
