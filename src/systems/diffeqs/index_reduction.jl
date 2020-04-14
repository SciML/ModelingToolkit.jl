function get_vnodes(sys)
    diffnodes = Operation[]
    edges = map(_->Int[], 1:length(sys.eqs))
    for (i, eq) in enumerate(sys.eqs)
        if !(eq.lhs isa Constant)
            # Make sure that the LHS is a first order derivative of a var.
            @assert eq.lhs.op isa Differential
            @assert !(eq.lhs.args[1] isa Differential) # first order

            push!(diffnodes, eq.lhs)
            # For efficiency we note down the diff edges here
            push!(edges[i], length(diffnodes))
        end
    end

    diffvars = (first ∘ var_from_nested_derivative).(diffnodes)
    algvars  = setdiff(states(sys), diffvars)
    return diffnodes, edges, algvars
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
                    push!(edges[i], j+varnumber_offset)
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


function matching_equation!(edges, i, assignments, active, vcolor=falses(length(active)), ecolor=falses(length(edges)))
    # `edge[active]` are active edges
    # i: variables
    # j: equations
    # assignments: assignments[j] == i means (i-j) is assigned
    #
    # color the equation
    ecolor[i] = true
    # if a V-node j exists s.t. edge (i-j) exists and assignments[j] == 0
    for j in edges[i]
        if active[j] && assignments[j] == 0
            assignments[j] = i
            return true
        end
    end
    # for every j such that edge (i-j) exists and j is uncolored
    for j in edges[i]
        (active[j] && !vcolor[j]) || continue
        # color the variable
        vcolor[j] = true
        if match_equation!(edges, assignments[j], assignments, active, vcolor, ecolor)
            assignments[v] = i
            return true
        end
    end
    return false
end

function matching(edges, nvars, active=trues(nvars))
    assignments = zeros(Int, nvars)
    for i in 1:length(edges)
        matching_equation!(edges, i, assignments, active)
    end
    return assignments
end
