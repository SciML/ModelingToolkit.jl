# V-nodes `[x_1, x_2, x_3, ..., dx_1, dx_2, ..., y_1, y_2, ...]` where `x`s are
# differential variables and `y`s are algebraic variables.
function get_vnodes(sys)
    dxvars = Operation[]
    edges = map(_->Int[], 1:length(sys.eqs))
    for (i, eq) in enumerate(sys.eqs)
        if !(eq.lhs isa Constant)
            # Make sure that the LHS is a first order derivative of a var.
            @assert eq.lhs.op isa Differential
            @assert !(eq.lhs.args[1] isa Differential) # first order

            push!(dxvars, eq.lhs)
            # For efficiency we note down the diff edges here
            push!(edges[i], length(dxvars))
        end
    end

    xvars = (first ∘ var_from_nested_derivative).(dxvars)
    algvars  = setdiff(states(sys), xvars)
    return xvars, dxvars, edges, algvars
end

function sys2bigraph(sys)
    xvars, dxvars, edges, algvars = get_vnodes(sys)
    xvar_offset = length(xvars)
    algvar_offset = 2xvar_offset
    for edge in edges
        isempty(edge) || (edge .+= xvar_offset)
    end

    for (i, eq) in enumerate(sys.eqs)
        # T or D(x):
        # We assume no derivatives appear on the RHS at this point
        vs = vars(eq.rhs)
        for v in vs
            for (j, target_v) in enumerate(xvars)
                if v == target_v
                    push!(edges[i], j)
                end
            end
            for (j, target_v) in enumerate(algvars)
                if v == target_v
                    push!(edges[i], j+algvar_offset)
                end
            end
        end
    end

    fullvars = [xvars; dxvars; algvars] # full list of variables
    vars_asso = [(1:xvar_offset) .+ xvar_offset; zeros(Int, length(fullvars) - xvar_offset)] # variable association list
    return edges, fullvars, vars_asso
end

print_bigraph(sys, vars, edges) = print_bigraph(stdout, sys, vars, edges)
function print_bigraph(io::IO, sys, vars, edges)
    println(io, "Equations:")
    foreach(x->println(io, x), [i => sys.eqs[i] for i in 1:length(sys.eqs)])
    for (i, edge) in enumerate(edges)
        println(io, "\nEq $i has:")
        print(io, '[')
        for e in edge
            print(io, "$(vars[e]), ")
        end
        print(io, ']')
    end
    return nothing
end

function match_equation!(edges, i, assign, active, vcolor=falses(length(active)), ecolor=falses(length(edges)))
    # `edge[active]` are active edges
    # i: equations
    # j: variables
    # assign: assign[j] == i means (i-j) is assigned
    #
    # color the equation
    ecolor[i] = true
    # if a V-node j exists s.t. edge (i-j) exists and assign[j] == 0
    for j in edges[i]
        if active[j] && assign[j] == 0
            assign[j] = i
            return true
        end
    end
    # for every j such that edge (i-j) exists and j is uncolored
    for j in edges[i]
        (active[j] && !vcolor[j]) || continue
        # color the variable
        vcolor[j] = true
        if match_equation!(edges, assign[j], assign, active, vcolor, ecolor)
            assign[j] = i
            return true
        end
    end
    return false
end

function matching(edges, nvars, active=trues(nvars))
    assign = zeros(Int, nvars)
    for i in 1:length(edges)
        match_equation!(edges, i, assign, active)
    end
    return assign
end

function pantelides(sys::ODESystem)
    edges, fullvars, vars_asso = sys2bigraph(sys)
    return pantelides!(edges, fullvars, vars_asso)
end

function pantelides!(edges, vars, vars_asso)
    neqs = length(edges)
    nvars = length(vars)
    assign = zeros(Int, nvars)
    eqs_asso = fill(0, neqs)
    neqs′ = neqs
    for k in 1:neqs′
        i = k
        pathfound = false
        iter = 0
        while !pathfound && iter < 3
            iter += 1
            # run matching on (dx, y) variables
            active = vars_asso .== 0
            vcolor = falses(nvars)
            ecolor = falses(neqs)
            pathfound = match_equation!(edges, i, assign, active, vcolor, ecolor)
            if !pathfound
                # for every colored V-node j
                for j in eachindex(vcolor); vcolor[j] || continue
                    # introduce a new variable
                    nvars += 1
                    push!(vars_asso, 0)
                    vars_asso[j] = nvars
                    push!(assign, 0)
                end

                # for every colored E-node l
                for l in eachindex(ecolor); ecolor[l] || continue
                    neqs += 1
                    # create new E-node
                    push!(edges, copy(edges[l]))
                    # create edges from E-node `neqs` to all V-nodes `j` and
                    # `vars_asso[j]` s.t. edge `(l-j)` exists
                    for j in edges[l]
                        if !(vars_asso[j] in edges[neqs])
                            push!(edges[neqs], vars_asso[j])
                        end
                    end
                    push!(eqs_asso, 0)
                    eqs_asso[l] = neqs
                end

                # for every colored V-node j
                for j in eachindex(vcolor); vcolor[j] || continue
                    assign[vars_asso[j]] = eqs_asso[assign[j]]
                end
                i = eqs_asso[i]
            end # if !pathfound
        end # while !pathfound
    end # for k in 1:neqs′
    return edges, assign, vars_asso, eqs_asso
end
