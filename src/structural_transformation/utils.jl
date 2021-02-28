###
### Bipartite graph utilities
###

"""
    find_augmenting_path(g::BipartiteGraph, eq, assign, varwhitelist, vcolor=falses(length(varwhitelist)), ecolor=falses(nsrcs(g))) -> path_found::Bool

Try to find augmenting paths.
"""
function find_augmenting_path(g, eq, assign, varwhitelist, vcolor=falses(length(varwhitelist)), ecolor=falses(nsrcs(g)))
    ecolor[eq] = true

    # if a `var` is unassigned and the edge `eq <=> var` exists
    for var in 𝑠neighbors(g, eq)
        if varwhitelist[var] && assign[var] == UNASSIGNED
            assign[var] = eq
            return true
        end
    end

    # for every `var` such that edge `eq <=> var` exists and `var` is uncolored
    for var in 𝑠neighbors(g, eq)
        (varwhitelist[var] && !vcolor[var]) || continue
        vcolor[var] = true
        if find_augmenting_path(g, assign[var], assign, varwhitelist, vcolor, ecolor)
            assign[var] = eq
            return true
        end
    end
    return false
end

"""
    matching(s::Union{SystemStructure,BipartiteGraph}, varwhitelist=trues(nvars), eqwhitelist=nothing) -> assign

Find equation-variable bipartite matching. `s.graph` is a bipartite graph.
`nvars` is the number of variables. Matched variables and equations are stored
in like `assign[var] => eq`.
"""
matching(s::SystemStructure, varwhitelist=trues(ndsts(s.graph)), eqwhitelist=nothing) = matching(s.graph, varwhitelist, eqwhitelist)
function matching(g::BipartiteGraph, varwhitelist=trues(ndsts(g)), eqwhitelist=nothing)
    assign = fill(UNASSIGNED, ndsts(g))
    for eq in 𝑠vertices(g)
        if eqwhitelist !== nothing
            eqwhitelist[eq] || continue
        end
        find_augmenting_path(g, eq, assign, varwhitelist)
    end
    return assign
end

###
### BLT ordering
###

"""
    find_scc(g::BipartiteGraph, assign=nothing)

Find strongly connected components of the equations defined by `g`. `assign`
gives the undirected bipartite graph a direction. When `assign === nothing`, we
assume that the ``i``-th variable is assigned to the ``i``-th equation.
"""
function find_scc(g::BipartiteGraph, assign=nothing)
    id = 0
    stack = Int[]
    components = Vector{Int}[]
    n = nsrcs(g)
    onstack = falses(n)
    lowlink = zeros(Int, n)
    ids = fill(UNVISITED, n)

    for eq in 𝑠vertices(g)
        if ids[eq] == UNVISITED
            id = strongly_connected!(stack, onstack, components, lowlink, ids, g, assign, eq, id)
        end
    end
    return components
end

"""
    strongly_connected!(stack, onstack, components, lowlink, ids, g, assign, eq, id)

Use Tarjan's algorithm to find strongly connected components.
"""
function strongly_connected!(stack, onstack, components, lowlink, ids, g, assign, eq, id)
    id += 1
    lowlink[eq] = ids[eq] = id

    # add `eq` to the stack
    push!(stack, eq)
    onstack[eq] = true

    # for `adjeq` in the adjacency list of `eq`
    for var in 𝑠neighbors(g, eq)
        if assign === nothing
            adjeq = var
        else
            # assign[var] => the equation that's assigned to var
            adjeq = assign[var]
            # skip equations that are not assigned
            adjeq == UNASSIGNED && continue
        end

        # if `adjeq` is not yet idsed
        if ids[adjeq] == UNVISITED # visit unvisited nodes
            id = strongly_connected!(stack, onstack, components, lowlink, ids, g, assign, adjeq, id)
        end
        # at the callback of the DFS
        if onstack[adjeq]
            lowlink[eq] = min(lowlink[eq], lowlink[adjeq])
        end
    end

    # if we are at a start of a strongly connected component
    if lowlink[eq] == ids[eq]
        component = Int[]
        repeat = true
        # pop until we are at the start of the strongly connected component
        while repeat
            w = pop!(stack)
            onstack[w] = false
            lowlink[w] = ids[eq]
            # put `w` in current component
            push!(component, w)
            repeat = w != eq
        end
        push!(components, sort!(component))
    end
    return id
end

function sorted_incidence_matrix(sys, val=true; only_algeqs=false, only_algvars=false)
    sys = algebraic_equations_scc(sys)
    s = structure(sys)
    @unpack assign, inv_assign, fullvars, scc, graph = s
    g = graph
    varsmap = zeros(Int, ndsts(graph))
    eqsmap = zeros(Int, nsrcs(graph))
    varidx = 0
    eqidx = 0
    for c in scc, eq in c
        var = inv_assign[eq]
        if var != 0
            eqsmap[eq] = (eqidx += 1)
            varsmap[var] = (varidx += 1)
        end
    end
    for i in diffvars_range(s)
        varsmap[i] = (varidx += 1)
    end
    for i in dervars_range(s)
        varsmap[i] = (varidx += 1)
    end
    for i in 1:nsrcs(graph)
        if eqsmap[i] == 0
            eqsmap[i] = (eqidx += 1)
        end
    end

    I = Int[]
    J = Int[]
    for eq in 𝑠vertices(g)
        only_algeqs && (isalgeq(s, eq) || continue)
        for var in 𝑠neighbors(g, eq)
            only_algvars && (isalgvar(s, var) || continue)
            push!(I, eqsmap[eq])
            push!(J, varsmap[var])
        end
    end
    #sparse(I, J, val, nsrcs(g), ndsts(g))
    sparse(I, J, val)
end

###
### Structural and symbolic utilities
###

function find_solvables!(sys)
    s = structure(sys)
    @unpack fullvars, graph, solvable_graph = s
    eqs = equations(sys)
    empty!(solvable_graph)
    for (i, eq) in enumerate(eqs)
        isdiffeq(eq) && continue
        term = value(eq.rhs - eq.lhs)
        for j in 𝑠neighbors(graph, i)
            isalgvar(s, j) || continue
            D = Differential(fullvars[j])
            c = expand_derivatives(D(term), false)
            if !(c isa Symbolic) && c isa Number && c != 0
                add_edge!(solvable_graph, i, j)
            end
        end
    end
    s
end

###
### Miscellaneous
###

function inverse_mapping(assign)
    invassign = zeros(Int, length(assign))
    for (i, eq) in enumerate(assign)
        eq <= 0 && continue
        invassign[eq] = i
    end
    return invassign
end

# debugging use
function reordered_matrix(sys, partitions=structure(sys).partitions)
    s = structure(sys)
    @unpack graph = s
    eqs = equations(sys)
    nvars = ndsts(graph)
    I, J = Int[], Int[]
    ii = 0
    M = Int[]
    for (e_solved, v_solved, e_residue, v_tear) in partitions
        append!(M, v_solved)
        append!(M, v_tear)
    end
    M = inverse_mapping(vcat(M, setdiff(1:nvars, M)))
    for (e_solved, v_solved, e_residue, v_tear) in partitions
        for es in e_solved
            isdiffeq(eqs[es]) && continue
            ii += 1
            js = [M[x] for x in 𝑠neighbors(graph, es) if isalgvar(s, x)]
            append!(I, fill(ii, length(js)))
            append!(J, js)
        end

        for er in e_residue
            isdiffeq(eqs[er]) && continue
            ii += 1
            js = [M[x] for x in 𝑠neighbors(graph, er) if isalgvar(s, x)]
            append!(I, fill(ii, length(js)))
            append!(J, js)
        end
    end
    # only plot algebraic variables and equations
    sparse(I, J, true)
end

###
### Nonlinear equation(s) solving
###

@noinline nlsolve_failure(rc) = error("The nonlinear solver failed with the return code $rc.")

function numerical_nlsolve(f, u0, p)
    prob = NonlinearProblem{false}(f, u0, p)
    sol = solve(prob, NewtonRaphson())
    rc = sol.retcode
    rc === :DEFAULT || nlsolve_failure(rc)
    # TODO: robust initial guess, better debugging info, and residual check
    sol.u
end
