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
    G                      # Vector{ Vector{Int} }
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
