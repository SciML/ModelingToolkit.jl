###
### BLT ordering
###

"""
    $(TYPEDSIGNATURES)

Obtain the incidence matrix of the system sorted by the SCCs. Requires that the system is
simplified and has a `schedule`.
"""
function sorted_incidence_matrix(sys::AbstractSystem)
    if !iscomplete(sys) || get_tearing_state(sys) === nothing ||
            get_schedule(sys) === nothing
        error("A simplified `System` is required. Call `mtkcompile` on the system before creating an `SCCNonlinearProblem`.")
    end
    sched = get_schedule(sys)
    var_sccs = sched.var_sccs

    ts = get_tearing_state(sys)
    imat = Graphs.incidence_matrix(ts.structure.graph)
    buffer = similar(imat)
    permute!(buffer, imat, 1:size(imat, 2), reduce(vcat, var_sccs))
    return buffer
end

###
### Structural and symbolic utilities
###
function highest_order_variable_mask(ts)
    return let v2d = ts.structure.var_to_diff
        v -> isempty(outneighbors(v2d, v))
    end
end

function lowest_order_variable_mask(ts)
    return let v2d = ts.structure.var_to_diff
        v -> isempty(inneighbors(v2d, v))
    end
end

function but_ordered_incidence(ts::TearingState, varmask = highest_order_variable_mask(ts))
    graph = complete(ts.structure.graph)
    var_eq_matching = complete(maximal_matching(graph; srcfilter = _ -> true, dstfilter = varmask))
    scc = StateSelection.find_var_sccs(graph, var_eq_matching)
    vordering = Vector{Int}(undef, 0)
    bb = Int[1]
    sizehint!(vordering, ndsts(graph))
    sizehint!(bb, ndsts(graph))
    l = 1
    for c in scc
        isemptyc = true
        for v in c
            if varmask(v)
                push!(vordering, v)
                l += 1
                isemptyc = false
            end
        end
        isemptyc || push!(bb, l)
    end
    mm = incidence_matrix(graph)
    reverse!(vordering)
    return mm[[var_eq_matching[v] for v in vordering if var_eq_matching[v] isa Int], vordering], bb
end

"""
    $(TYPEDSIGNATURES)

Given a system `sys` and torn variable-equation matching `torn_matching`, return the sparse
incidence matrix of the system with SCCs grouped together, and each SCC sorted to contain
the analytically solved equations/variables before the unsolved ones.
"""
function reordered_matrix(sys::System, torn_matching)
    s = TearingState(sys)
    StateSelection.complete!(s.structure)
    @unpack graph = s.structure
    eqs = equations(sys)
    nvars = ndsts(graph)
    max_matching = complete(maximal_matching(graph))
    torn_matching = complete(torn_matching)
    sccs = StateSelection.find_var_sccs(graph, max_matching)
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
            js = [M[x] for x in 𝑠neighbors(graph, es) if StateSelection.isalgvar(s.structure, x)]
            append!(I, fill(ii, length(js)))
            append!(J, js)
        end

        e_residual = setdiff(
            [
                max_matching[v]
                    for v in vars if max_matching[v] !== unassigned
            ], e_solved
        )
        for er in e_residual
            isdiffeq(eqs[er]) && continue
            ii += 1
            js = [M[x] for x in 𝑠neighbors(graph, er) if StateSelection.isalgvar(s.structure, x)]
            append!(I, fill(ii, length(js)))
            append!(J, js)
        end
    end
    # only plot algebraic variables and equations
    return sparse(I, J, true)
end

"""
    uneven_invmap(n::Int, list)

returns an uneven inv map with length `n`.
"""
function uneven_invmap(n::Int, list)
    rename = zeros(Int, n)
    for (i, v) in enumerate(list)
        rename[v] = i
    end
    return rename
end
