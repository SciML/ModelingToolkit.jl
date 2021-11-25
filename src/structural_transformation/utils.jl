###
### Bipartite graph utilities
###

"""
    maximal_matching(s::SystemStructure, eqfilter=eq->true, varfilter=v->true) -> Matching

Find equation-variable maximal bipartite matching. `s.graph` is a bipartite graph.
"""
BipartiteGraphs.maximal_matching(s::SystemStructure, eqfilter=eq->true, varfilter=v->true) =
    maximal_matching(s.graph, eqfilter, varfilter)

function error_reporting(sys, bad_idxs, n_highest_vars, iseqs)
    io = IOBuffer()
    if iseqs
        error_title = "More equations than variables, here are the potential extra equation(s):\n"
        out_arr = equations(sys)[bad_idxs]
    else
        error_title = "More variables than equations, here are the potential extra variable(s):\n"
        out_arr = structure(sys).fullvars[bad_idxs]
    end

    Base.print_array(io, out_arr)
    msg = String(take!(io))
    neqs = length(equations(sys))
    if iseqs
        throw(ExtraEquationsSystemException(
            "The system is unbalanced. "
            * "There are $n_highest_vars highest order derivative variables "
            * "and $neqs equations.\n"
            * error_title
            * msg
        ))
    else
        throw(ExtraVariablesSystemException(
            "The system is unbalanced. "
            * "There are $n_highest_vars highest order derivative variables "
            * "and $neqs equations.\n"
            * error_title
            * msg
        ))
    end
end

###
### Structural check
###
function check_consistency(sys::AbstractSystem)
    s = structure(sys)
    @unpack graph, var_to_diff, fullvars = s
    n_highest_vars = count(v->length(outneighbors(s.var_to_diff, v)) == 0, vertices(s.var_to_diff))
    neqs = nsrcs(graph)
    is_balanced = n_highest_vars == neqs

    if neqs > 0 && !is_balanced
        varwhitelist = var_to_diff .== nothing
        var_eq_matching = maximal_matching(graph, eq->true, v->varwhitelist[v]) # not assigned
        # Just use `error_reporting` to do conditional
        iseqs = n_highest_vars < neqs
        if iseqs
            eq_var_matching = invview(complete(var_eq_matching)) # extra equations
            bad_idxs = findall(isnothing, @view eq_var_matching[1:nsrcs(graph)])
        else
            bad_idxs = findall(isequal(unassigned), var_eq_matching)
        end
        error_reporting(sys, bad_idxs, n_highest_vars, iseqs)
    end

    # This is defined to check if Pantelides algorithm terminates. For more
    # details, check the equation (15) of the original paper.
    extended_graph = (@set graph.fadjlist = Vector{Int}[graph.fadjlist; map(collect, edges(var_to_diff))])
    extended_var_eq_matching = maximal_matching(extended_graph)

    unassigned_var = []
    for (vj, eq) in enumerate(extended_var_eq_matching)
        if eq === unassigned
            push!(unassigned_var, fullvars[vj])
        end
    end

    if !isempty(unassigned_var) || !is_balanced
        io = IOBuffer()
        Base.print_array(io, unassigned_var)
        unassigned_var_str = String(take!(io))
        errmsg = "The system is structurally singular! " *
                 "Here are the problematic variables: \n" *
                 unassigned_var_str
        throw(InvalidSystemException(errmsg))
    end

    return nothing
end

###
### BLT ordering
###

"""
    find_var_sccs(g::BipartiteGraph, assign=nothing)

Find strongly connected components of the variables defined by `g`. `assign`
gives the undirected bipartite graph a direction. When `assign === nothing`, we
assume that the ``i``-th variable is assigned to the ``i``-th equation.
"""
function find_var_sccs(g::BipartiteGraph, assign=nothing)
    cmog = DiCMOBiGraph{true}(g, Matching(assign === nothing ? Base.OneTo(nsrcs(g)) : assign))
    sccs = Graphs.strongly_connected_components(cmog)
    foreach(sort!, sccs)
    return sccs
end

function sorted_incidence_matrix(sys, val=true; only_algeqs=false, only_algvars=false)
    var_eq_matching, var_scc = algebraic_variables_scc(sys)
    s = structure(sys)
    @unpack fullvars, graph = s
    g = graph
    varsmap = zeros(Int, ndsts(graph))
    eqsmap = zeros(Int, nsrcs(graph))
    varidx = 0
    eqidx = 0
    for vs in scc, v in vs
        eq = var_eq_matching[v]
        if eq !== unassigned
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
            i = eqsmap[eq]
            j = varsmap[var]
            (iszero(i) || iszero(j)) && continue
            push!(I, i)
            push!(J, j)
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
            var = fullvars[j]
            isinput(var) && continue
            a, b, islinear = linear_expansion(term, var)
            a = unwrap(a)
            if islinear && (!(a isa Symbolic) && a isa Number && a != 0)
                add_edge!(solvable_graph, i, j)
            end
        end
    end
    s
end

# debugging use
function reordered_matrix(sys, torn_matching)
    s = structure(sys)
    @unpack graph = s
    eqs = equations(sys)
    nvars = ndsts(graph)
    max_matching = complete(maximal_matching(s))
    torn_matching = complete(torn_matching)
    sccs = find_var_sccs(s.graph, max_matching)
    I, J = Int[], Int[]
    ii = 0
    M = Int[]
    solved = BitSet(findall(torn_matching .!== unassigned))
    for vars in sccs
        append!(M, filter(in(solved), vars))
        append!(M, filter(!in(solved), vars))
    end
    M = invperm(vcat(M, setdiff(1:nvars, M)))
    for vars in sccs
        e_solved = [torn_matching[v] for v in vars if torn_matching[v] !== unassigned]
        for es in e_solved
            isdiffeq(eqs[es]) && continue
            ii += 1
            js = [M[x] for x in 𝑠neighbors(graph, es) if isalgvar(s, x)]
            append!(I, fill(ii, length(js)))
            append!(J, js)
        end

        e_residual = setdiff([max_matching[v] for v in vars if max_matching[v] !== unassigned], e_solved)
        for er in e_residual
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
