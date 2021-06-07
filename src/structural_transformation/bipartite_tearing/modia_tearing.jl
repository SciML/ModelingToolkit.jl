# This code is from the Modia project and is licensed as follows:
# https://github.com/ModiaSim/Modia.jl/blob/b61daad643ef7edd0c1ccce6bf462c6acfb4ad1a/LICENSE

################################################
#
# Functions to tear systems of equations
#
# Author: Martin Otter, DLR-SR (first version: Jan. 14, 2017)
#
# Details are described in the paper:
#   Otter, Elmqvist (2017): Transformation of Differential Algebraic Array Equations to
#                           Index One Form. Modelica'2017 Conference.
#
# The following utility algorithm is used below to incrementally added edges to a
# DAG (Directed Acyclic Graph). This algorithm leads to an O(n*m) worst time complexity of the
# tearing (instead of O(m*m)) where n is the number of equations and m is the number of
# variable incidences. Note, the traversals of the DAG are not performed with recursive function
# calls but with while loops and an explicit stack, in order to avoid function stack overflow
# for large algebraic loops.
#
#   Bender, Fineman, Gilbert, Tarjan:
#      A New Approach to Incremental Cycle Detection and Related Problems.
#      ACM Transactions on Algorithms, Volume 12, Issue 2, Feb. 2016
#      http://dl.acm.org/citation.cfm?id=2756553
#
#   Text excerpt from this paper (the advantage of Algorithm N is that it
#   is simple to implement and needs no sophisticated data structures)
#
#       3. A ONE-WAY-SEARCH ALGORITHM FOR DENSE GRAPHS
#           [..]
#       To develop the algorithm, we begin with a simple algorithm and then modify it to
#       improve its running time. We call the simple algorithm Algorithm N (for “naive”). The
#       algorithm initializes k(v) = 1 for each vertex v. The initial vertex levels are a weak
#       topological numbering since the initial graph contains no arcs. To add an arc (v,w),
#       if k(v) < k(w), merely add the arc. If, on the other hand, k(v) ≥ k(w), add the arc and
#       then do a selective forward search that begins by traversing (v,w) and continues until
#       the search traverses an arc into v (there is a cycle) or there are no more candidate
#       arcs to traverse. To traverse an arc (x, y), if y = v, stop (there is a cycle); otherwise, if
#       k(x) ≥ k(y), increase k(y) to k(x)+1 and add all arcs (y, z) to the set of candidate arcs to
#       traverse.
#       It is easy to show that (1) after each arc addition that does not form a cycle
#       the vertex levels are a weak topological numbering, (2) Algorithm N is correct, and
#       (3) 1 ≤ k(v) ≤ size(v) ≤ n for all v. Since an arc (x, y) notEqual (v,w) is only traversed as a
#       result of k(x) increasing, each arc is traversed at most n times, resulting in a total time
#       bound of O(nm) for all arc additions.
#
################################################

const Undefined = typemin(Int)


"""
    td = TraverseDAG(G,nv)

Generate an object td to traverse a set of equations that are
represented as a DAG (Directed Acyclic Graph).
G is the bipartite graph of all relevant equations
and nv is the largest variable number used in G (or larger).
"""
mutable struct TraverseDAG
    minlevel::Int
    curlevel::Int
    level::Vector{Int}
    lastlevel::Vector{Int}
    levelStack::Vector{Int}
    visitedStack::Vector{Int}
    vActive::Vector{Bool}
    visited::Vector{Bool}
    check::Vector{Bool}
    stack::Vector{Int}
    eSolved::Vector{Int}
    vSolved::Vector{Int}
    G::Vector{ Vector{Int} }
    assign::Vector{Int}
    es::Vector{Int}
    vs::Vector{Int}

    function TraverseDAG(G, nv::Int)
        visited      = fill(false, length(G))
        check        = fill(false, length(G))
        vActive      = fill(false, nv)
        level        = fill(Undefined, nv)
        lastlevel    = fill(Undefined, nv)
        levelStack   = fill(0, 0)
        stack        = fill(0, 0)
        visitedStack = fill(0, 0)
        eSolved      = fill(0, 0)
        vSolved      = fill(0, 0)
        assign       = fill(0, nv)

        new(0, Undefined, level, lastlevel, levelStack, visitedStack, vActive,
          visited, check, stack, eSolved, vSolved, G, assign)
    end
end


"""
    initAlgebraicSystem(td::TraverseDAG,es,vs)

Define the set of equations and the set variables for which the equations shall be solved for
(equations es shall be solved for variables vs) and re-initialize td.

eSolvedFixed/vSolvedFixed must be a DAG starting at eSolvedFixed/SolvedFixed[1]
"""
function initAlgebraicSystem(td::TraverseDAG, es::Vector{Int}, vs::Vector{Int},
                             eSolvedFixed::Vector{Int}, vSolvedFixed::Vector{Int}, vTearFixed::Vector{Int})
    # check arguments
    for i in eachindex(es)
        if es[i] <= 0
            error("\n\n... Internal error in Tearing.jl: es[", i, "] = ", es[i], ".\n")
        end
    end

    for i in eachindex(vs)
        if vs[i] <= 0
            error("\n\n... Internal error in Tearing.jl: vs[", i, "] = ", vs[i], ".\n")
        end
    end

    # check that all elements of eSolvedFixed are in es, vSolvedFixed in vs,
    # vTearFixed in vs and that vSolvedFixed and vTearFixed have no variables in common
    @assert( length(eSolvedFixed) == length(vSolvedFixed) )
    ediff = setdiff(eSolvedFixed, es)
    @assert(length(ediff) == 0)
    vdiff = setdiff(vSolvedFixed, vs)
    @assert(length(vdiff) == 0)
    vdiff2 = setdiff(vTearFixed, vs)
    @assert(length(vdiff2) == 0)
    vdiff3 = intersect(vSolvedFixed, vTearFixed)
    @assert(length(vdiff3) == 0)

    # Re-initialize td
    td.minlevel = 0
    td.curlevel = Undefined
    for i in eachindex(td.visited)
        td.visited[i] = false
        td.check[i]   = false
    end

    for i in eachindex(td.vActive)
        td.vActive[i]   = false
        td.assign[i]    = 0
        td.level[i]     = Undefined
        td.lastlevel[i] = Undefined
    end

    for i in eachindex(vs)
        td.vActive[ vs[i] ] = true
    end

    empty!(td.levelStack)
    empty!(td.stack)
    empty!(td.visitedStack)
    empty!(td.eSolved)
    empty!(td.vSolved)

    # Define initial DAG
    vs2 = Int[]
    for i in eachindex(vSolvedFixed)
        vFixed = vSolvedFixed[i]
        td.assign[vFixed] = eSolvedFixed[i]
        td.level[ vFixed] = i
        push!(vs2, vFixed)
    end

    for i in eachindex(vTearFixed)
        td.vActive[ vTearFixed[i] ] = false   # vTearFixed shall not be assigned
    end

    # Store es, vs in td
    td.es = es
    td.vs = vs

    return vs2
end


in_vs(td, v) = td.vActive[v]

function setlevel(td::TraverseDAG, v::Int, parentLevel::Int)
    td.lastlevel[v] = td.level[v]
    td.level[v]     = parentLevel + 1
    push!(td.visitedStack, v)
end


"""
    success = visit!(td::TraverseDAG, v)

Traverse potential DAG starting from new variable node v.
If no cycle is detected return true, otherwise return false.
"""
function visit!(td::TraverseDAG, vcheck::Int)
    empty!(td.stack)
    empty!(td.levelStack)
    empty!(td.visitedStack)
    td.curlevel = td.level[vcheck]
    push!(td.levelStack, td.curlevel)
    push!(td.stack, vcheck)
    first = true

    while length(td.stack) > 0
        parentLevel = pop!(td.levelStack)
        veq = pop!(td.stack)
        eq  = td.assign[veq]
        if first
            first = false
        else
            if td.level[veq] == td.curlevel
            # cycle detected
                return false
            elseif td.level[veq] == Undefined || td.level[veq] <= parentLevel
                setlevel(td, veq, parentLevel)
            end
        end

        if eq > 0
            # Push all child nodes on stack
            parentLevel = td.level[veq]
            for v in td.G[eq]
                if in_vs(td, v) && v != veq   # v is an element of td.vs and is not the variable to solve for
                    eq2 = td.assign[v]
                    if eq2 == 0 || td.level[v] <= parentLevel
                        push!(td.levelStack, parentLevel)
                        push!(td.stack, v)
                    end
                end
            end
        end
    end

    return true
end


"""
    visit2!(td::TraverseDAG,v)

Traverse DAG starting from variable v and store visited equations and variables in stacks
eSolved, vSolved. If a cycle is deteced, raise an error (signals a programming error).
"""
function visit2!(td::TraverseDAG, vVisit::Int)
    push!(td.stack, vVisit)
    while length(td.stack) > 0
        veq = td.stack[end]
        eq  = td.assign[veq]
        if !td.visited[eq]
            td.visited[eq] = true
            td.check[eq]   = true
            for v in td.G[eq]
                if in_vs(td, v) && v != veq  # v is an element of td.vs and is not the variable to solve for
                    eq2 = td.assign[v]
                    if eq2 != 0
                        if !td.visited[eq2]   # visit eq2 if not yet visited
                            push!(td.stack, v)
                        elseif td.check[eq2]  # cycle detected
                            error("... error in Tearing.jl code: \n",
                           "    cycle detected (should not occur): eq = ", eq, ", veq = ", veq, ", eq2 = ", eq2, ", v = ", v)
                        end
                    end
                end
            end
        else
            td.check[eq] = false
            push!(td.eSolved, eq)
            push!(td.vSolved, veq)
            pop!(td.stack)
        end
    end
    nothing
end


"""
    (eSolved, vSolved) = sortDAG!(td::TraverseDAG, vs)

Sort the equations that are assigned by variables vs using object td of type TraverseDAG
and return the sorted equations eSolved and assigned variables vSolved.
"""
function sortDAG!(td::TraverseDAG, vs::Vector{Int})
    # initialize data structure
    empty!(td.stack)
    empty!(td.eSolved)
    empty!(td.vSolved)

    for i in eachindex(td.visited)
        td.visited[i] = false
        td.check[i]   = false
    end

    # visit all assigned variables and equations
    for veq in vs
        if !td.visited[ td.assign[veq] ]
            visit2!(td, veq)
        end
    end

    return (td.eSolved, td.vSolved)
end


"""
    (eSolved, vSolved, eResidue, vTear) = tearEquations!(td, Gsolvable, es, vs; eSolvedFixed=Int[], vSolvedFixed=Int[], vTearFixed=Int[])

Equations es shall be solved with respect to variables vs. The function returns
the teared equation so that if vTear is given, vSolved can be computed from eSolved
in a forward sequence (so solving eSolved[1] for vSolved[1], eSolved[2] for vSolved[2],
and so on). vTear must be selected, so that the equations eResidues are fulfilled.
Equations es are the union of eSolved and eResidue.
Variables vs are the union of vSolved and vTear.

Input argument td is an object of type TraverseDAG. Gsolvable defines the variables
that can be explicitly solved in every equation without influencing the solution space
(= rank preserving operation).

eSolvedFixed/vSolvedFixed must be a DAG starting at eSolvedFixed/SolvedFixed[1]
"""
function tearEquations!(td::TraverseDAG, Gsolvable, es::Vector{Int}, vs::Vector{Int};
                        eSolvedFixed::Vector{Int}=Int[], vSolvedFixed::Vector{Int}=Int[], vTearFixed::Vector{Int}=Int[])
    vs2 = initAlgebraicSystem(td, es, vs, eSolvedFixed, vSolvedFixed, vTearFixed)
    # eResidue = fill(0,0)
    residue  = true
    esReduced = setdiff(es, eSolvedFixed)
    # println("    es = ", es, ", eSolvedFixed = ", eSolvedFixed, ", esReduced = ", esReduced)
    # println("    vs = ", vs, ", vSolvedFixed = ", vSolvedFixed)
    for eq in esReduced  # iterate only over equations that are not in eSolvedFixed
        residue = true
        for vj in Gsolvable[eq]
            if td.assign[vj] == 0 && in_vs(td, vj)
                # vj is an element of vs that is not yet assigned
                # Add equation to graph
                td.assign[vj] = eq

                # Check for cycles
                if td.level[vj] == Undefined
                    # (eq,vj) cannot introduce a cycle
                    # Introduce a new level (the smallest level that exists yet)
                    td.minlevel += -1
                    td.level[vj] = td.minlevel

                    # Inspect all childs and use level+1, if child has no level yet
                    for v in td.G[eq]
                        if in_vs(td, v) && v != vj &&
                     (td.level[v] == Undefined || td.level[v] <= td.level[vj]) # v is an element of td.vs and is not the variable to solve for and no level yet defined
                            setlevel(td, v, td.level[vj])
                        end
                    end

                    push!(vs2, vj)
                    residue = false
                    break # continue with next equation

                else # Traverse DAG starting from eq
                    if visit!(td, vj)
                        # accept vj
                        push!(vs2, vj)
                        residue = false
                        break   # continue with next equation
                    else
                        # cycle; remove vj from DAG and undo its changes
                        for vv in td.visitedStack
                            td.level[vv] = td.lastlevel[vv]
                        end
                        td.assign[vj] = 0
                        # continue with next variable in equation eq
                    end
                end
            end
        end
      #if residue
      #   push!(eResidue, eq)
      #end
    end

    # Determine solved equations and variables
    (eSolved, vSolved) = sortDAG!(td, vs2)
    vTear = setdiff(vs, vSolved)
    eResidue = setdiff(es, eSolved)
    return (eSolved, vSolved, eResidue, vTear)
end

# New Tearing

# License for this file: MIT (expat)
# Copyright 2017-2021, DLR Institute of System Dynamics and Control
# Author: Martin Otter, DLR-SR
#
# Functions to tear systems of equations.

const Undefined = typemin(Int)


#=
Bender, Fineman, Gilbert, Tarjan (2016):
A New Approach to Incremental Cycle Detection and Related Problems.

ACM Transactions on Algorithms, Volume 12, Issue 2, Feb. 2016
http://dl.acm.org/citation.cfm?id=2756553

Text excerpt from this paper (the advantage of Algorithm N is that it
is simple to implement and needs no sophisticated data structures)

3. A ONE-WAY-SEARCH ALGORITHM FOR DENSE GRAPHS
[..]
To develop the algorithm, we begin with a simple algorithm and then modify it to
improve its running time. We call the simple algorithm Algorithm N (for “naive”). The
algorithm initializes k(v) = 1 for each vertex v. The initial vertex levels are a weak
topological numbering since the initial graph contains no arcs. To add an arc (v,w),
if k(v) < k(w), merely add the arc. If, on the other hand, k(v) ≥ k(w), add the arc and
then do a selective forward search that begins by traversing (v,w) and continues until
the search traverses an arc into v (there is a cycle) or there are no more candidate
arcs to traverse. To traverse an arc (x, y), if y = v, stop (there is a cycle); otherwise, if
k(x) ≥ k(y), increase k(y) to k(x)+1 and add all arcs (y, z) to the set of candidate arcs to
traverse.
It is easy to show that (1) after each arc addition that does not form a cycle
the vertex levels are a weak topological numbering, (2) Algorithm N is correct, and
(3) 1 ≤ k(v) ≤ size(v) ≤ n for all v. Since an arc (x, y) notEqual (v,w) is only traversed as a
result of k(x) increasing, each arc is traversed at most n times, resulting in a total time
bound of O(nm) for all arc additions.
=#


"""
    ts = TearingSetup(G [,nv])

Generate a setup object `ts` from a bi-partite graph `G`
to reduce one or more systems of equations that are present in `G`.
The returned object `ts` holds internal auxiliary storage that is
used in analysis with function [`tearEquations!`](@ref).

`G` is expected to be a vector of integer vectors
(for example of type `Vector{ Vector{Int} }`).
The optional argument `nv` is the largest integer number occuring in `G`
(= number of variables). If `nv` is not provided, it is deduced from `G`.

`TearingSetup` is useful to reduce the amount of memory to be allocated,
if function [`tearEquations!`](@ref) shall be called several times on
equation systems from `G`.
"""
mutable struct TearingSetup
    minlevel::Int
    curlevel::Int
    level::Vector{Int}
    lastlevel::Vector{Int}
    levelStack::Vector{Int}
    visitedStack::Vector{Int}
    vActive::Vector{Bool}
    vAssignable::Vector{Bool}
    visited::Vector{Bool}
    check::Vector{Bool}
    stack::Vector{Int}
    eSolved::Vector{Int}
    vSolved::Vector{Int}
    G::Vector{Vector{Int}}
    assign::Vector{Int}

    function TearingSetup(G::AbstractVector, nv::Int)
        visited      = fill(false, length(G))
        check        = fill(false, length(G))
        vActive      = fill(false, nv)
        vAssignable  = fill(false, nv)
        level        = fill(Undefined, nv)
        lastlevel    = fill(Undefined, nv)
        levelStack   = fill(0, 0)
        stack        = fill(0, 0)
        visitedStack = fill(0, 0)
        eSolved      = fill(0, 0)
        vSolved      = fill(0, 0)
        assign       = fill(0, nv)

        new(0, Undefined, level, lastlevel, levelStack, visitedStack, vActive,
            vAssignable, visited, check, stack, eSolved, vSolved, G, assign)
    end

    function TearingSetup(G::AbstractVector)
        nv = 1
        for g in G
            nv = max(nv, maximum(g))
        end
        TearingSetup(G,nv)
    end
end


"""
    reInitializeTearingSetup!(td::TearingSetup, es, vs, eSolvedFixed, vSolvedFixed)

Define the set of equations and the set of variables for which the equations shall be solved for
(equations es shall be solved for variables vs) and re-initialize td.

eSolvedFixed/vSolvedFixed must be a DAG starting at eSolvedFixed/vSolvedFixed[1]
"""
function reInitializeTearingSetup!(td::TearingSetup, es::Vector{Int}, vs::Vector{ Vector{Int} },
        eSolvedFixed::Vector{Int}, vSolvedFixed::Vector{Int}, check::Bool)
    # check arguments
    if check
        # Check that eSolvedFixed, vSolvedFixed are a DAG
        @assert( length(eSolvedFixed) == length(vSolvedFixed) )
        if length(eSolvedFixed) > 0
            checkDAG!(td, eSolvedFixed, vSolvedFixed)
        end

        # Check es
        neq = length(td.G)
        for esi in es
            @assert(esi > 0 && esi <= neq)
        end

        # Check that eSolvedFixed and es have no elements in common
        ediff = intersect(eSolvedFixed, es)
        if length(ediff) != 0
            error("... Error in Tearing.jl: Function tearEquations!(...) was not called correctly:\n",
                  "eSolvedFixed and es are not allowed to have elements in common.")
        end

        # Check vs and that vSolvedFixed and vs have no elements in common
        nv = length(td.vActive)
        for vse in vs
            for vsi in vse
                @assert(vsi > 0 && vsi <= nv)
            end
            vdiff = intersect(vSolvedFixed, vse)
            if length(vdiff) != 0
                error("... Error in Tearing.jl: Function tearEquations!(...) was not called correctly:\n",
                      "vSolvedFixed and vs are not allowed to have elements in common.")
            end
        end
    end

    # Re-initialize td
    td.minlevel = 0
    td.curlevel = Undefined
    for i in eachindex(td.visited)
        td.visited[i] = false
        td.check[i]   = false
    end

    for i in eachindex(td.vActive)
        td.vActive[i]     = false
        td.vAssignable[i] = false
        td.assign[i]      = 0
        td.level[i]       = Undefined
        td.lastlevel[i]   = Undefined
    end

    empty!(td.levelStack)
    empty!(td.stack)
    empty!(td.visitedStack)
    empty!(td.eSolved)
    empty!(td.vSolved)

    # Define initial DAG
    vs2 = Int[]
    for i in eachindex(vSolvedFixed)
        vFixed = vSolvedFixed[i]
        td.assign[vFixed] = eSolvedFixed[i]
        td.level[ vFixed] = i
        push!(vs2, vFixed)
    end

    # Define that vs and vSolvedFixed shall be treated as the unknown variables during traversal
    for vse in vs
        for vsi in vse
            td.vActive[vsi] = true
        end
    end
    for vsi in vSolvedFixed
        td.vActive[vsi] = true
    end

    return vs2
end


in_vActive(td, v)     = td.vActive[v]
in_vAssignable(td, v) = td.vAssignable[v]


"""
    result = isDAG!(td::TearingSetup, vStart::Int)

Traverse potential DAG starting from variable node vStart.
If no cycle is detected return true, otherwise return false.
"""
function isDAG!(td::TearingSetup, vStart::Int)::Bool
    # Initialize stacks and flags
    empty!(td.stack)
    for i in eachindex(td.visited)
        td.visited[i] = false
        td.check[i]   = false
    end

    # Traverse DAG
    push!(td.stack, vStart)
    while length(td.stack) > 0
        veq = td.stack[end]
        eq  = td.assign[veq]
        if !td.visited[eq]
            td.visited[eq] = true
            td.check[eq]   = true
            for v in td.G[eq]
                if in_vActive(td, v) && v != veq  # v is an element of the unknown variables and is not the variable to solve for
                    eq2 = td.assign[v]
                    if eq2 != 0
                        if !td.visited[eq2]
                            # eq2 not yet visited
                            push!(td.stack, v)
                        elseif td.check[eq2]
                            # cycle detected
                            return false
                        end
                        #else
                        #    error("... error in Tearing.jl code in function isDAG! (should not occur): assign[$v] = 0")
                    end
                end
            end
        else
            td.check[eq] = false
            pop!(td.stack)
        end
    end
    return true
end



"""
    visitDAG!(td::TearingSetup,v)

Traverse DAG starting from variable v and store visited equations and variables in stacks
eSolved, vSolved. If a cycle is deteced, raise an error (signals a programming error).
"""
function visitDAG!(td::TearingSetup, vVisit::Int)
    push!(td.stack, vVisit)
    while length(td.stack) > 0
        veq = td.stack[end]
        eq  = td.assign[veq]
        if !td.visited[eq]
            td.visited[eq] = true
            td.check[eq]   = true
            for v in td.G[eq]
                if in_vActive(td, v) && v != veq  # v is an element of the unknown variables and is not the variable to solve for
                    eq2 = td.assign[v]
                    if eq2 != 0
                        if !td.visited[eq2]   # visit eq2 if not yet visited
                            push!(td.stack, v)
                        elseif td.check[eq2]  # cycle detected
                            error("... error in Tearing.jl code: \n",
                                  "    cycle detected (should not occur): eq = ", eq, ", veq = ", veq, ", eq2 = ", eq2, ", v = ", v)
                        end
                    end
                end
            end
        else
            td.check[eq] = false
            push!(td.eSolved, eq)
            push!(td.vSolved, veq)
            pop!(td.stack)
        end
    end
    nothing
end



"""
    (eSolved, vSolved) = sortDAG!(td::TearingSetup, vs)

Sort the equations that are assigned by variables vs using object td of type TearingSetup
and return the sorted equations eSolved and assigned variables vSolved.
"""
function sortDAG!(td::TearingSetup, vs::Vector{Int})
    # initialize data structure
    empty!(td.stack)
    empty!(td.eSolved)
    empty!(td.vSolved)

    for i in eachindex(td.visited)
        td.visited[i] = false
        td.check[i]   = false
    end

    # visit all assigned variables and equations
    for veq in vs
        if !td.visited[ td.assign[veq] ]
            visitDAG!(td, veq)
        end
    end

    return (td.eSolved, td.vSolved)
end


"""
    checkDAG!(td::TearingSetup, es::Vector{Int}, vs::Vector{Int})

A DAG is defined by equations es and variables vs, such that es[1] is solved for vs[1],
afterwards es[2] for vs[2] etc.

This function checks whether es/vs defines a DAG.
If a cycle is detected, an error is raised (signals a programming error).
"""
function checkDAG!(td::TearingSetup, es::Vector{Int}, vs::Vector{Int})::Nothing
    if length(vs) == 0
        return nothing
    end

    # Initialize stacks and flags
    empty!(td.stack)
    for i in eachindex(td.visited)
        td.visited[i]  = false
        td.check[i]    = false
        td.assign[i] = 0
    end
    for i in eachindex(td.vActive)
        td.vActive[i] = false
    end
    for vsi in vs
        td.vActive[vsi] = true
    end
    for i in eachindex(es)
        td.assign[vs[i]] = es[i]
    end

    # Traverse DAG
    push!(td.stack, vs[1])
    while length(td.stack) > 0
        veq = td.stack[end]
        eq  = td.assign[veq]
        if !td.visited[eq]
            td.visited[eq] = true
            td.check[eq]   = true
            for v in td.G[eq]
                if in_vActive(td, v) && v != veq  # v is an element of the unknown variables and is not the variable to solve for
                    eq2 = td.assign[v]
                    if eq2 != 0
                        if !td.visited[eq2]   # visit eq2 if not yet visited
                            push!(td.stack, v)
                        elseif td.check[eq2]  # cycle detected
                            error("... error in Tearing.jl code (should not occur): \n",
                                  "    Cycle detected: eq = ", eq, ", veq = ", veq, ", eq2 = ", eq2, ", v = ", v)
                        end
                    else
                        error("... error in Tearing.jl code (should not occur): assign[$v] = 0")
                    end
                end
            end
        else
            td.check[eq] = false
            pop!(td.stack)
        end
    end
    return nothing
end




function tearEquationsCore!(vs2, td::TearingSetup, isSolvable::Function, es::Vector{Int}, vs::Vector{Int}, log::Bool)::Nothing
    G = td.G
    vAssignable = false
    for eq in es  # iterate only over equations that are not in eSolvedFixed
        for vj in G[eq]
            vAssignable = in_vAssignable(td, vj)
            vAssigned   = td.assign[vj] > 0
            if log
                if !vAssigned
                    if vAssignable
                        if isSolvable(eq,vj)
                            println("        Equation $eq can be solved for variable $vj")
                        else
                            println("        Equation $eq cannot be solved for variable $vj")
                        end
                    end
                end
            end
            if vAssignable && !vAssigned && isSolvable(eq,vj)
                # vj is an element of the variables that can be assigned in this stage, but is not yet assigned
                # Add equation to graph
                td.assign[vj] = eq

                # Check for cycles (traverse DAG starting from vj)
                if isDAG!(td, vj)
                    # accept vj
                    push!(vs2, vj)
                    if log
                        println("            -> solved for variable $vj without cycles")
                    end
                    break   # continue with next equation
                else
                    # cycle; remove vj from DAG and undo its changes
                    td.assign[vj] = 0
                    if log
                        println("            -> solving for variable $vj gives a cycle (so not done)")
                    end
                    # continue with next variable in equation eq
                end
            end
        end
    end

    return nothing
end



"""
    (eSolved, vSolved, eResidue, vTear) = tearEquations!(GorTs, isSolvable, es, vs;
        eSolvedFixed=Int[], vSolvedFixed=Int[], check=true)

This function tears a system of equations consisting of equations `es` and `eSolvedFixed`
that are functions of variables `vs` and `vSolvedFixed`.
The optional arguments `eSolvedFixed, vSolvedFixed`
define the starting Directed Acyclic Graph (solving equations eSolvedFixed[1] for variable vSolvedFixed[1],
eSolvedFixed[2] for variable vSolvedFixed[2] etc.) starting at `vSolvedFixed[1]`.

The function returns
the teared equations so that if `vTear` are given, `vSolved` can be computed from `eSolved`
in a forward sequence (so solving `eSolved[1]` for `vSolved[1]`, `eSolved[2]` for `vSolved[2]`,
and so on). `vTear` must be selected, so that the equations `eResidues` are fulfilled.
Equations `eSolved` and `eResidue` are the (reordered) union of `es` and `eSolvedFixed`.
Variables vSolved` and `vTear` are the (reordered) union of `vs` and `vSolvedFixed`.

This means that an algebraic equation system `0 = g(w)`
(`g` are equations `es` and `eSolvedFixed`;
 `w` are unknowns `vs` and `vSolvedFixed`) is solved as
much as possible explicitly for the unknowns resulting in

```
w_e := g_e(w_t)
  0  = g_r(w_t, w_e)
```

where

- `w` (= vs and vSolvedFixed) consists of all elements of `w_e, w_t`;
- equations `g_e` of `g` are explicitly solved for `w_e`;
- equations `g_r` are the equations of `g` that cannot be explicitly solved (= residual equations).


# Required arguments

- `GorTs`: Either a bi-partite graph `G::Vector{Vector{Int}}` or `ts::TearingSetup`
           generated with constructor `ts = TearingSetup(G)`.
           `ts` is re-initialized in `tearEquations!` for any call of the function.
           `TearingSetup` is useful to reduce the amount of memory to be allocated,
           if several equation systems shall be teared from `G`.

- `isSolvable(e,v)`: Function that returns true, if equation `e`
           can be solved for variable `v` without influencing the solution space
           (= rank preserving operation).

- `es::Vector{Int}`: Vector of equations that shall be solved with respect to variables `vs`.
           `es` must be equations from `G`.

- `vs`: Either of type `Vector{Int}` or of type `Vector{Vector{Int}}`.
            `vs` are the unknown variables that shall be solved from `es`.
            If `vs::Vector{Vector{Int}}`, it is first tried to solve `es`
            for `vs[1]`, then for `vs[2]` etc.
            (so `vs` defines a priority to solve for variables).

# Optional arguments

- `eSolvedFixed::Vector{Int}`: Equations of `G` that are already defined to be solved for `vSolvedFixed`.
- `vSolvedFixed::Vector{Int]`: Variables of `G` that are explicitly solved from `eSolvedFixed`.
- `check::Bool`: = true, if various checks shall be performed, for example
            that eSolvedFixed/vSolvedFixed and eSolved/vSolved are a DAG respectively.

# Return arguments

- `eSolved::Vector{Int}`: Equations that are explicitly solved in the order `eSolved[1], eSolved[2], ...`.
- `vSolved::Vector{Int}`: Equation `eSolved[i]` is explicitly solved for variable `vSolved[i]`.
- `eResdiue::Vector{Int}`: Residual equations that are not explicitly solved.
- `vTear::Vector{Int}`: Tearing variables, so variables that are assumed to be known, when solving
                        equations `eSolved`.

# Example

```
using ModiaBase

G = Vector{Int}[ [1,2,4],  # equation 1 depends on variables 1,2,4
                 [1,7],
                 [3,4],
                 [3,7],
                 [6,7],
                 [2] ]

es = [3,4,2,1]    # Solve equations 3,4,2,1
vs = [3,7,1,4]    # for variables 3,7,1,4
isSolvable(e,v) = true  # Assume that every equation is solvable for its unknowns

(eSolved,vSolved,eResidue,vTear) = tearEquations!(G, isSolvable, es, vs)

# eSolved  = [3,4,2]
# vSolved  = [3,7,1]
# eResidue = [1]
# vTear    = [4]
```

# Algorithm

The algorithm used in this function is sketched in the paper:

- Otter, Elmqvist (2017):
  [Transformation of Differential Algebraic Array Equations to Index One Form](http://www.ep.liu.se/ecp/132/064/ecp17132565.pdf).
  Modelica'2017 Conference.

The function uses several extensions of the described basic
tearing algorithm that are important for transforming higher index Differential Algebraic Equations
to index one form. Note, the traversals in Directed-Acyclic-Graphs - the core operation of
the tearing algorithm - is **not** performed with recursive function calls but with
while loops and an explicit stack, in order to avoid function stack overflow
for large algebraic loops. Tests up to 1 million equations in 1 million unknowns have been
performed.

For improving efficiency, algorithm N of the following paper is used as utility algorithm:

- Bender, Fineman, Gilbert, Tarjan (2016):
  [A New Approach to Incremental Cycle Detection and Related Problems](http://dl.acm.org/citation.cfm?id=2756553).
  ACM Transactions on Algorithms, Volume 12, Issue 2, Feb.


# Main developer

[Martin Otter](https://rmc.dlr.de/sr/en/staff/martin.otter/),
[DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en)
"""
tearEquations!(G               , isSolvable, es, vs             ; kwargs...) = tearEquations!(TearingSetup(G), isSolvable, es, vs; kwargs...)
tearEquations!(td::TearingSetup, isSolvable, es, vs::Vector{Int}; kwargs...) = tearEquations!(td             , isSolvable, es, fill(vs,1); kwargs...)
function tearEquations!(td::TearingSetup, isSolvable::Function, es::Vector{Int}, vs::Vector{ Vector{Int} };
        eSolvedFixed::Vector{Int}=Int[], vSolvedFixed::Vector{Int}=Int[],
        check::Bool=true, log::Bool=false)
    # Reinitialize td
    vs2 = reInitializeTearingSetup!(td, es, vs, eSolvedFixed, vSolvedFixed, check)

    # Try vs in order of priority
    for vse in vs
        if length(vse) > 0
            # Define variables that can be assigned
            for vsi in vse
                td.vAssignable[ vsi ] = true
            end

            # Perform tearing
            if log
                println("    Try to solve for variable(s) ", vse)
            end
            tearEquationsCore!(vs2, td, isSolvable, es, vse, log)

            # Undefine variables that can be assigned
            for vsi in vse
                td.vAssignable[ vsi ] = false
            end
        end
    end

    # Determine solved equations and variables
    (eSolved, vSolved) = sortDAG!(td, vs2)

    eResidue = setdiff(es, eSolved)
    if length(vs) == 1
        vTear = setdiff(vs[1], vSolved)
    else
        vs_vector = deepcopy(vs[1])
        for i in 2:length(vs)
            append!(vs_vector, vs[i])
        end
        vTear = setdiff(vs_vector, vSolved)
    end

    # Check result
    if check
        checkDAG!(td, eSolved, vSolved)
    end

    return (eSolved, vSolved, eResidue, vTear)
end
