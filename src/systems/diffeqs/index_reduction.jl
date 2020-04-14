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

function get_vnodes(sys)
    diffnodes = []
    diffedges = Tuple{Int, Int}[]
    for (i, eq) in enumerate(sys.eqs)
        if !(eq.lhs isa Constant)
            # Make sure that the LHS is a first order derivative of a var.
            @assert eq.lhs.op isa Differential
            @assert !(eq.lhs.args[1] isa Differential) # first order

            push!(diffnodes, eq.lhs)
            # For efficiency we note down the diff edges here
            push!(diffedges, (i, length(diffnodes)))
        end
    end

    diffvars = (first ∘ var_from_nested_derivative).(diffnodes)
    algvars  = setdiff(states(sys), diffvars)
    return diffnodes, diffedges, algvars
end

function sys2bigraph2(sys)
    diffvars, edges, algvars = get_vnodes(sys)
    varnumber_offset = length(diffvars)

    for (i, eq) in enumerate(sys.eqs)
        # T or D(x):
        # We assume no derivatives appear on the RHS at this point
        vs = vars(eq.rhs)
        for v in vs
            for (j, target_v) in enumerate(algvars)
                if v == target_v
                    push!(edges, (i, j+varnumber_offset))
                end
            end
        end
    end
    vcat(diffvars, algvars), edges
end

print_bigraph(sys, vars, edges) = print_bigraph(stdout, sys, vars, edges)
function print_bigraph(io::IO, sys, vars, edges)
    println(io, "Equations:")
    foreach(x->println(io, x), [i => sys.eqs[i] for i in 1:length(sys.eqs)])
    println(io)
    for (i, j) in edges
        println(io, "Eq $i has $(vars[j])")
    end
end
