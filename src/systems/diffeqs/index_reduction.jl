struct BiGraph{T}
    data::Vector{Vector{T}}
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

function sys2bigraph(sys)
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


function matching_equation!(edges, i, assignments=Dict{Int, Int}(), colored=Set{Int}())
    push!(colored, i)
    # select a v
    vars = unique(last.(filter(isequal(i)∘first, edges)))
    for v in vars
        if !haskey(assignments, v)# && !(v in colored)
            # v found
            assignments[v] = i
            return true
        end
    end
    # Else
    remaining = setdiff(vars, colored)
    for v in remaining
        push!(colored, v)
        if match_equation!(edges, assignments[v], colored, assignments)
            assignments[v] = i
            return true
        end
    end
    return false
end

function matching(sys, vars, edges)
    assignments=Dict{Int, Int}()
    colored=Set{Int}()
    for i in 1:length(sys.eqs)
        @show matching_equation!(edges, i, assignments, colored)
    end
    assignments
end
