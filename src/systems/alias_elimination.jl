using SymbolicUtils: Rewriters

const KEEP = typemin(Int)

function alias_eliminate_graph!(state::TransformationState)
    mm = linear_subsys_adjmat(state)
    size(mm, 1) == 0 && return AliasGraph(ndsts(state.structure.graph)), mm, BitSet() # No linear subsystems

    @unpack graph, var_to_diff = state.structure

    return alias_eliminate_graph!(complete(graph), complete(var_to_diff), mm)
end

# For debug purposes
function aag_bareiss(sys::AbstractSystem)
    state = TearingState(sys)
    mm = linear_subsys_adjmat(state)
    return aag_bareiss!(state.structure.graph, complete(state.structure.var_to_diff), mm)
end

function extreme_var(var_to_diff, v, level = nothing, ::Val{descend} = Val(true);
                     callback = _ -> nothing) where {descend}
    g = descend ? invview(var_to_diff) : var_to_diff
    callback(v)
    while (v′ = g[v]) !== nothing
        v::Int = v′
        callback(v)
        if level !== nothing
            descend ? (level -= 1) : (level += 1)
        end
    end
    level === nothing ? v : (v => level)
end

alias_elimination(sys) = alias_elimination!(TearingState(sys; quick_cancel = true))
function alias_elimination!(state::TearingState)
    Main._state[] = deepcopy(state)
    sys = state.sys
    complete!(state.structure)
    ag, mm, updated_diff_vars = alias_eliminate_graph!(state)
    isempty(ag) && return sys

    fullvars = state.fullvars
    @unpack var_to_diff, graph = state.structure

    if !isempty(updated_diff_vars)
        has_iv(sys) ||
            error(InvalidSystemException("The system has no independent variable!"))
        D = Differential(get_iv(sys))
        for v in updated_diff_vars
            dv = var_to_diff[v]
            fullvars[dv] = D(fullvars[v])
        end
    end

    subs = Dict()
    obs = Equation[]
    # If we encounter y = -D(x), then we need to expand the derivative when
    # D(y) appears in the equation, so that D(-D(x)) becomes -D(D(x)).
    to_expand = Int[]
    diff_to_var = invview(var_to_diff)
    for (v, (coeff, alias)) in pairs(ag)
        lhs = fullvars[v]
        rhs = iszero(coeff) ? 0 : coeff * fullvars[alias]
        subs[lhs] = rhs
        v != alias && push!(obs, lhs ~ rhs)
        if coeff == -1
            # if `alias` is like -D(x)
            diff_to_var[alias] === nothing && continue
            # if `v` is like y, and D(y) also exists
            (dv = var_to_diff[v]) === nothing && continue
            # all equations that contains D(y) needs to be expanded.
            append!(to_expand, 𝑑neighbors(graph, dv))
        end
    end

    dels = Int[]
    eqs = collect(equations(state))
    for (ei, e) in enumerate(mm.nzrows)
        vs = 𝑠neighbors(graph, e)
        if isempty(vs)
            # remove empty equations
            push!(dels, e)
        else
            rhs = mapfoldl(+, pairs(nonzerosmap(@view mm[ei, :]))) do (var, coeff)
                iszero(coeff) && return 0
                return coeff * fullvars[var]
            end
            eqs[e] = 0 ~ rhs
        end
    end
    deleteat!(eqs, sort!(dels))
    old_to_new = Vector{Int}(undef, length(var_to_diff))
    idx = 0
    cursor = 1
    ndels = length(dels)
    for i in eachindex(old_to_new)
        if cursor <= ndels && i == dels[cursor]
            cursor += 1
            old_to_new[i] = -1
            continue
        end
        idx += 1
        old_to_new[i] = idx
    end

    lineqs = BitSet(old_to_new[e] for e in mm.nzrows)
    for (ieq, eq) in enumerate(eqs)
        ieq in lineqs && continue
        eqs[ieq] = substitute(eq, subs)
    end

    for old_ieq in to_expand
        ieq = old_to_new[old_ieq]
        eqs[ieq] = expand_derivatives(eqs[ieq])
    end

    newstates = []
    diff_to_var = invview(var_to_diff)
    for j in eachindex(fullvars)
        if !(j in keys(ag))
            diff_to_var[j] === nothing && push!(newstates, fullvars[j])
        end
    end

    sys = state.sys
    @set! sys.eqs = eqs
    @set! sys.states = newstates
    @set! sys.observed = [observed(sys); obs]
    return invalidate_cache!(sys)
end

"""
$(SIGNATURES)

Find the first linear variable such that `𝑠neighbors(adj, i)[j]` is true given
the `constraint`.
"""
@inline function find_first_linear_variable(M::SparseMatrixCLIL,
                                            range,
                                            mask,
                                            constraint)
    eadj = M.row_cols
    for i in range
        vertices = eadj[i]
        if constraint(length(vertices))
            for (j, v) in enumerate(vertices)
                (mask === nothing || mask[v]) &&
                    return (CartesianIndex(i, v), M.row_vals[i][j])
            end
        end
    end
    return nothing
end

@inline function find_first_linear_variable(M::AbstractMatrix,
                                            range,
                                            mask,
                                            constraint)
    for i in range
        row = @view M[i, :]
        if constraint(count(!iszero, row))
            for (v, val) in enumerate(row)
                iszero(val) && continue
                if mask === nothing || mask[v]
                    return CartesianIndex(i, v), val
                end
            end
        end
    end
    return nothing
end

function find_masked_pivot(variables, M, k)
    r = find_first_linear_variable(M, k:size(M, 1), variables, isequal(1))
    r !== nothing && return r
    r = find_first_linear_variable(M, k:size(M, 1), variables, isequal(2))
    r !== nothing && return r
    r = find_first_linear_variable(M, k:size(M, 1), variables, _ -> true)
    return r
end

"""
    AliasGraph

When eliminating variables, keeps track of which variables where eliminated in
favor of which others.

Currently only supports elimination as direct aliases (+- 1).

We represent this as a dict from eliminated variables to a (coeff, var) pair
representing the variable that it was aliased to.
"""
struct AliasGraph <: AbstractDict{Int, Pair{Int, Int}}
    aliasto::Vector{Union{Int, Nothing}}
    eliminated::Vector{Int}
    function AliasGraph(nvars::Int)
        new(fill(nothing, nvars), Int[])
    end
end

Base.length(ag::AliasGraph) = length(ag.eliminated)

function Base.getindex(ag::AliasGraph, i::Integer)
    r = ag.aliasto[i]
    r === nothing && throw(KeyError(i))
    coeff, var = (sign(r), abs(r))
    nc = coeff
    av = var
    # We support `x -> -x` as an alias.
    if var != i && var in keys(ag)
        # Amortized lookup. Check if since we last looked this up, our alias was
        # itself aliased. If so, just adjust the alias table.
        ac, av = ag[var]
        nc = ac * coeff
        ag.aliasto[i] = nc > 0 ? av : -av
    end
    return (nc, av)
end

function Base.iterate(ag::AliasGraph, state...)
    r = Base.iterate(ag.eliminated, state...)
    r === nothing && return nothing
    c = ag.aliasto[r[1]]
    return (r[1] => (c == 0 ? 0 :
                     c >= 0 ? 1 :
                     -1, abs(c))), r[2]
end

function Base.setindex!(ag::AliasGraph, ::Nothing, i::Integer)
    if ag.aliasto[i] !== nothing
        ag.aliasto[i] = nothing
        deleteat!(ag.eliminated, findfirst(isequal(i), ag.eliminated))
    end
end
function Base.setindex!(ag::AliasGraph, v::Integer, i::Integer)
    @assert v == 0
    if ag.aliasto[i] === nothing
        push!(ag.eliminated, i)
    end
    ag.aliasto[i] = 0
    return 0 => 0
end

function Base.setindex!(ag::AliasGraph, p::Union{Pair{Int, Int}, Tuple{Int, Int}}, i::Integer)
    (c, v) = p
    if c == 0 || v == 0
        ag[i] = 0
        return p
    end
    @assert v != 0 && c in (-1, 1)
    if ag.aliasto[i] === nothing
        push!(ag.eliminated, i)
    end
    ag.aliasto[i] = c > 0 ? v : -v
    return p
end

function Base.get(ag::AliasGraph, i::Integer, default)
    i in keys(ag) || return default
    return ag[i]
end

struct AliasGraphKeySet <: AbstractSet{Int}
    ag::AliasGraph
end
Base.keys(ag::AliasGraph) = AliasGraphKeySet(ag)
Base.iterate(agk::AliasGraphKeySet, state...) = Base.iterate(agk.ag.eliminated, state...)
function Base.in(i::Int, agk::AliasGraphKeySet)
    aliasto = agk.ag.aliasto
    1 <= i <= length(aliasto) && aliasto[i] !== nothing
end

function reduce!(mm::SparseMatrixCLIL, ag::AliasGraph)
    for i in 1:size(mm, 1)
        adj_row = @view mm[i, :]
        locally_structure_simplify!(adj_row, nothing, ag)
    end
    mm
end

struct InducedAliasGraph
    ag::AliasGraph
    invag::SimpleDiGraph{Int}
    var_to_diff::DiffGraph
    visited::BitSet
end

function tograph(ag::AliasGraph, var_to_diff::DiffGraph)
    g = SimpleDiGraph{Int}(length(var_to_diff))
    zero_vars = Int[]
    for (v, (_, a)) in ag
        if iszero(a)
            push!(zero_vars, v)
            continue
        end
        add_edge!(g, v, a)
        add_edge!(g, a, v)
    end
    transitiveclosure!(g)
    zero_vars_set = BitSet(zero_vars)
    for v in zero_vars
        for a in outneighbors(g, v)
            push!(zero_vars_set, a)
        end
    end
    # Compute the largest transitive closure that doesn't include any diff
    # edges.
    og = g
    newg = SimpleDiGraph{Int}(length(var_to_diff))
    for e in Graphs.edges(og)
        s, d = src(e), dst(e)
        (var_to_diff[s] == d || var_to_diff[d] == s) && continue
        oldg = copy(newg)
        add_edge!(newg, s, d)
        add_edge!(newg, d, s)
        transitiveclosure!(newg)
        if any(e->(var_to_diff[src(e)] == dst(e) || var_to_diff[dst(e)] == src(e)), edges(newg))
            newg = oldg
        end
    end
    g = newg

    c = "green"
    edge_styles = Dict{Tuple{Int, Int}, String}()
    for (v, dv) in enumerate(var_to_diff)
        dv isa Int || continue
        edge_styles[(v, dv)] = c
        add_edge!(g, v, dv)
        add_edge!(g, dv, v)
    end
    g, zero_vars_set, edge_styles
end

using Graphs.Experimental.Traversals
struct DiffLevelState <: Traversals.AbstractTraversalState
    dists::Vector{Int}
    var_to_diff::DiffGraph
    visited::BitSet
end

DiffLevelState(g::SimpleDiGraph, var_to_diff) = DiffLevelState(fill(typemax(Int), nv(g)), var_to_diff, BitSet())

@inline function Traversals.initfn!(s::DiffLevelState, u)
    push!(s.visited, u)
    s.dists[u] = 0
    return true
end

@inline function Traversals.newvisitfn!(s::DiffLevelState, u, v)
    push!(s.visited, v)
    w = s.var_to_diff[u] == v ? 1 : s.var_to_diff[v] == u ? -1 : 0
    s.dists[v] = s.dists[u] + w
    return true
end

function find_root!(ss::DiffLevelState, g, s)
    Traversals.traverse_graph!(g, s, Traversals.BFS(), ss)
    argmin(Base.Fix1(getindex, ss.dists), ss.visited)
end

function get_levels(g, var_to_diff, s)
    ss = DiffLevelState(g, var_to_diff)
    Traversals.traverse_graph!(g, s, Traversals.BFS(), ss)
    return dists
end

function InducedAliasGraph(ag, invag, var_to_diff)
    InducedAliasGraph(ag, invag, var_to_diff, BitSet())
end

struct IAGNeighbors
    iag::InducedAliasGraph
    v::Int
end

function Base.iterate(it::IAGNeighbors, state = nothing)
    @unpack ag, invag, var_to_diff, visited = it.iag
    callback! = Base.Fix1(push!, visited)
    if state === nothing
        v, lv = extreme_var(var_to_diff, it.v, 0)
        used_ag = false
        nb = neighbors(invag, v)
        nit = iterate(nb)
        state = (v, lv, used_ag, nb, nit)
    end

    v, level, used_ag, nb, nit = state
    v in visited && return nothing
    while true
        @label TRYAGIN
        if used_ag
            if nit !== nothing
                n, ns = nit
                if !(n in visited)
                    n, lv = extreme_var(var_to_diff, n, level)
                    extreme_var(var_to_diff, n, nothing, Val(false), callback = callback!)
                    nit = iterate(nb, ns)
                    return n => lv, (v, level, used_ag, nb, nit)
                end
            end
        else
            used_ag = true
            # We don't care about the alising value because we only use this to
            # find the root of the tree.
            if (_n = get(ag, v, nothing)) !== nothing && (n = _n[2]) > 0
                if !(n in visited)
                    n, lv = extreme_var(var_to_diff, n, level)
                    extreme_var(var_to_diff, n, nothing, Val(false), callback = callback!)
                    return n => lv, (v, level, used_ag, nb, nit)
                end
            else
                @goto TRYAGIN
            end
        end
        push!(visited, v)
        (v = var_to_diff[v]) === nothing && return nothing
        level += 1
        used_ag = false
    end
end

Graphs.neighbors(iag::InducedAliasGraph, v::Integer) = IAGNeighbors(iag, v)

function _find_root!(iag::InducedAliasGraph, v::Integer, level = 0)
    brs = neighbors(iag, v)
    min_var_level = v => level
    for (x, lv′) in brs
        lv = lv′ + level
        x, lv = _find_root!(iag, x, lv)
        if min_var_level[2] > lv
            min_var_level = x => lv
        end
    end
    x, lv = extreme_var(iag.var_to_diff, min_var_level...)
    return x => lv
end

function find_root!(iag::InducedAliasGraph, v::Integer)
    clear_visited!(iag)
    _find_root!(iag, v)
end

clear_visited!(iag::InducedAliasGraph) = (empty!(iag.visited); iag)

struct RootedAliasTree
    iag::InducedAliasGraph
    root::Int
end

if Base.isdefined(AbstractTrees, :childtype)
    AbstractTrees.childtype(::Type{<:RootedAliasTree}) = Union{RootedAliasTree, Int}
else
    childtype(::Type{<:RootedAliasTree}) = Union{RootedAliasTree, Int}
end
AbstractTrees.children(rat::RootedAliasTree) = RootedAliasChildren{false}(rat)
AbstractTrees.children(rat::RootedAliasTree, ::Val{C}) where C = RootedAliasChildren{C}(rat)
AbstractTrees.nodetype(::Type{<:RootedAliasTree}) = Int
if Base.isdefined(AbstractTrees, :nodevalue)
    AbstractTrees.nodevalue(rat::RootedAliasTree) = rat.root
else
    nodevalue(rat::RootedAliasTree) = rat.root
    nodevalue(a) = a
end
if Base.isdefined(AbstractTrees, :shouldprintkeys)
    AbstractTrees.shouldprintkeys(rat::RootedAliasTree) = false
else
    shouldprintkeys(rat::RootedAliasTree) = false
end
has_fast_reverse(::Type{<:AbstractSimpleTreeIter{<:RootedAliasTree}}) = false

struct StatefulAliasBFS{T} <: AbstractSimpleTreeIter{T}
    t::T
end
# alias coefficient, depth, children
Base.eltype(::Type{<:StatefulAliasBFS{T}}) where {T} = Tuple{Int, Int, childtype(T)}
function Base.iterate(it::StatefulAliasBFS, queue = (eltype(it)[(1, 0, it.t)]))
    isempty(queue) && return nothing
    coeff, lv, t = popfirst!(queue)
    nextlv = lv + 1
    for (coeff′, c) in children(t)
        # -1 <= coeff <= 1
        push!(queue, (coeff * coeff′, nextlv, c))
    end
    return (coeff, lv, t), queue
end

struct RootedAliasChildren{C}
    t::RootedAliasTree
end

function Base.iterate(c::RootedAliasChildren, s = nothing)
    rat = c.t
    @unpack iag, root = rat
    @unpack visited = iag
    push!(visited, root)
    it = _iterate(c, s)
    it === nothing && return nothing
    while true
        node = nodevalue(it[1][2])
        if node in visited
            it = _iterate(c, it[2])
            it === nothing && return nothing
        else
            push!(visited, node)
            return it
        end
    end
end

@inline function _iterate(c::RootedAliasChildren{C}, s = nothing) where C
    rat = c.t
    @unpack iag, root = rat
    @unpack ag, invag, var_to_diff = iag
    iszero(root) && return nothing
    if !C
        root = var_to_diff[root]
    end
    root === nothing && return nothing
    if s === nothing
        stage = 1
        it = iterate(neighbors(invag, root))
        s = (stage, it)
    end
    (stage, it) = s
    if stage == 1 # root
        stage += 1
        return (1, root), (stage, it)
    elseif stage == 2 # ag
        stage += 1
        cv = get(ag, root, nothing)
        if cv !== nothing
            return (cv[1], RootedAliasTree(iag, cv[2])), (stage, it)
        end
    end
    # invag (stage 3)
    it === nothing && return nothing
    e, ns = it
    # c * a = b <=> a = c * b when -1 <= c <= 1
    return (ag[e][1], RootedAliasTree(iag, e)), (stage,
                                                 iterate(neighbors(invag, root), ns))
end

count_nonzeros(a::AbstractArray) = count(!iszero, a)

# N.B.: Ordinarily sparse vectors allow zero stored elements.
# Here we have a guarantee that they won't, so we can make this identification
count_nonzeros(a::SparseVector) = nnz(a)

# Linear variables are highest order differentiated variables that only appear
# in linear equations with only linear variables. Also, if a variable's any
# derivaitves is nonlinear, then all of them are not linear variables.
function find_linear_variables(graph, linear_equations, var_to_diff, irreducibles)
    stack = Int[]
    linear_variables = falses(length(var_to_diff))
    var_to_lineq = Dict{Int, BitSet}()
    mark_not_linear! = let linear_variables = linear_variables, stack = stack,
        var_to_lineq = var_to_lineq

        v -> begin
            linear_variables[v] = false
            push!(stack, v)
            while !isempty(stack)
                v = pop!(stack)
                eqs = get(var_to_lineq, v, nothing)
                eqs === nothing && continue
                for eq in eqs, v′ in 𝑠neighbors(graph, eq)
                    if linear_variables[v′]
                        linear_variables[v′] = false
                        push!(stack, v′)
                    end
                end
            end
        end
    end
    for eq in linear_equations, v in 𝑠neighbors(graph, eq)
        linear_variables[v] = true
        vlineqs = get!(() -> BitSet(), var_to_lineq, v)
        push!(vlineqs, eq)
    end
    for v in irreducibles
        lv = extreme_var(var_to_diff, v)
        while true
            mark_not_linear!(lv)
            lv = var_to_diff[lv]
            lv === nothing && break
        end
    end

    linear_equations_set = BitSet(linear_equations)
    for (v, islinear) in enumerate(linear_variables)
        islinear || continue
        lv = extreme_var(var_to_diff, v)
        oldlv = lv
        remove = invview(var_to_diff)[v] !== nothing
        while !remove
            for eq in 𝑑neighbors(graph, lv)
                if !(eq in linear_equations_set)
                    remove = true
                end
            end
            lv = var_to_diff[lv]
            lv === nothing && break
        end
        lv = oldlv
        remove && while true
            mark_not_linear!(lv)
            lv = var_to_diff[lv]
            lv === nothing && break
        end
    end

    return linear_variables
end

function aag_bareiss!(graph, var_to_diff, mm_orig::SparseMatrixCLIL, irreducibles = ())
    mm = copy(mm_orig)
    linear_equations = mm_orig.nzrows

    # If linear highest differentiated variables cannot be assigned to a pivot,
    # then we can set it to zero. We use `rank1` to track this.
    #
    # We only use alias graph to record reducible variables. We use `rank2` to
    # track this.
    #
    # For all the other variables, we can update the original system with
    # Bareiss'ed coefficients as Gaussian elimination is nullspace perserving
    # and we are only working on linear homogeneous subsystem.
    is_reducible = trues(length(var_to_diff))
    for v in irreducibles
        is_reducible[v] = false
    end
    # TODO/FIXME: This needs a proper recursion to compute the transitive
    # closure.
    is_linear_variables = find_linear_variables(graph, linear_equations, var_to_diff,
                                                irreducibles)
    solvable_variables = findall(is_linear_variables)

    function do_bareiss!(M, Mold = nothing)
        rank1 = rank2 = nothing
        pivots = Int[]
        function find_pivot(M, k)
            if rank1 === nothing
                r = find_masked_pivot(is_linear_variables, M, k)
                r !== nothing && return r
                rank1 = k - 1
            end
            if rank2 === nothing
                r = find_masked_pivot(is_reducible, M, k)
                r !== nothing && return r
                rank2 = k - 1
            end
            # TODO: It would be better to sort the variables by
            # derivative order here to enable more elimination
            # opportunities.
            return find_masked_pivot(nothing, M, k)
        end
        function find_and_record_pivot(M, k)
            r = find_pivot(M, k)
            r === nothing && return nothing
            push!(pivots, r[1][2])
            return r
        end
        function myswaprows!(M, i, j)
            Mold !== nothing && swaprows!(Mold, i, j)
            swaprows!(M, i, j)
        end
        bareiss_ops = ((M, i, j) -> nothing, myswaprows!,
                       bareiss_update_virtual_colswap_mtk!, bareiss_zero!)
        rank3, = bareiss!(M, bareiss_ops; find_pivot = find_and_record_pivot)
        rank2 = something(rank2, rank3)
        rank1 = something(rank1, rank2)
        (rank1, rank2, rank3, pivots)
    end

    return mm, solvable_variables, do_bareiss!(mm, mm_orig)
end

# Kind of like the backward substitution, but we don't actually rely on it
# being lower triangular. We eliminate a variable if there are at most 2
# variables left after the substitution.
function lss(mm, pivots, ag)
    ei -> let mm = mm, pivots = pivots, ag = ag
        vi = pivots[ei]
        locally_structure_simplify!((@view mm[ei, :]), vi, ag)
    end
end

function simple_aliases!(ag, graph, var_to_diff, mm_orig, irreducibles = ())
    mm, solvable_variables, (rank1, rank2, rank3, pivots) = aag_bareiss!(graph, var_to_diff,
                                                                         mm_orig,
                                                                         irreducibles)

    # Step 2: Simplify the system using the Bareiss factorization
    rk1vars = BitSet(@view pivots[1:rank1])
    for v in solvable_variables
        v in rk1vars && continue
        ag[v] = 0
    end

    echelon_mm = copy(mm)
    lss! = lss(mm, pivots, ag)
    # Step 2.1: Go backwards, collecting eliminated variables and substituting
    #         alias as we go.
    foreach(lss!, reverse(1:rank2))

    # Step 2.2: Sometimes Bareiss can make the equations more complicated.
    #         Go back and check the original matrix. If this happened,
    #         Replace the equation by the one from the original system,
    #         but be sure to also run `lss!` again, since we only ran that
    #         on the Bareiss'd matrix, not the original one.
    reduced = mapreduce(|, 1:rank2; init = false) do ei
        if count_nonzeros(@view mm_orig[ei, :]) < count_nonzeros(@view mm[ei, :])
            mm[ei, :] = @view mm_orig[ei, :]
            return lss!(ei)
        end
        return false
    end

    # Step 2.3: Iterate to convergence.
    #         N.B.: `lss!` modifies the array.
    # TODO: We know exactly what variable we eliminated. Starting over at the
    #       start is wasteful. We can lookup which equations have this variable
    #       using the graph.
    reduced && while any(lss!, 1:rank2)
    end

    return mm, echelon_mm
end

function mark_processed!(processed, var_to_diff, v)
    diff_to_var = invview(var_to_diff)
    processed[v] = true
    while (v = diff_to_var[v]) !== nothing
        processed[v] = true
    end
    return nothing
end

function is_self_aliasing((v, (_, a)), var_to_diff)
    iszero(a) && return false
    v = extreme_var(var_to_diff, v)
    while true
        v == a && return true
        v = var_to_diff[v]
        v === nothing && break
    end
    return false
end

function Base.filter(f, ag::AliasGraph)
    newag = AliasGraph(length(ag.aliasto))
    for (v, ca) in ag
        if f(v => ca)
            newag[v] = ca
        end
    end
    newag
end

function alias_eliminate_graph!(graph, var_to_diff, mm_orig::SparseMatrixCLIL)
    # Step 1: Perform Bareiss factorization on the adjacency matrix of the linear
    #         subsystem of the system we're interested in.
    #
    nvars = ndsts(graph)
    ag = AliasGraph(nvars)
    mm, echelon_mm = simple_aliases!(ag, graph, var_to_diff, mm_orig)
    fullvars = Main._state[].fullvars

    # Step 3: Handle differentiated variables
    # At this point, `var_to_diff` and `ag` form a tree structure like the
    # following:
    #
    #         x   -->   D(x)
    #         ⇓          ⇑
    #         ⇓         x_t   -->   D(x_t)
    #         ⇓               |---------------|
    # z --> D(z)  --> D(D(z))  |--> D(D(D(z))) |
    #         ⇑               |---------------|
    # k --> D(k)
    #
    # where `-->` is an edge in `var_to_diff`, `⇒` is an edge in `ag`, and the
    # part in the box are purely conceptual, i.e. `D(D(D(z)))` doesn't appear in
    # the system. We call the variables in the box "virtual" variables.
    #
    # To finish the algorithm, we backtrack to the root differentiation chain.
    # If the variable already exists in the chain, then we alias them
    # (e.g. `x_t ⇒ D(D(z))`), else, we substitute and update `var_to_diff`.
    #
    # Note that since we always prefer the higher differentiated variable and
    # with a tie breaking strategy, the root variable (in this case `z`) is
    # always uniquely determined. Thus, the result is well-defined.
    dag = AliasGraph(nvars) # alias graph for differentiated variables
    updated_diff_vars = Int[]
    diff_to_var = invview(var_to_diff)
    processed = falses(nvars)
    g, zero_vars_set = tograph(ag, var_to_diff)
    dls = DiffLevelState(g, var_to_diff)
    is_diff_edge = let var_to_diff = var_to_diff
        (v, w) -> var_to_diff[v] == w || var_to_diff[w] == v
    end
    diff_aliases = Vector{Pair{Int, Int}}[]
    for (v, dv) in enumerate(var_to_diff)
        processed[v] && continue
        (dv === nothing && diff_to_var[v] === nothing) && continue
        r = find_root!(dls, g, v)
        @show fullvars[r]
        level_to_var = Int[]
        extreme_var(var_to_diff, r, nothing, Val(false),
                    callback = Base.Fix1(push!, level_to_var))
        nlevels = length(level_to_var)
        prev_r = -1
        for _ in 1:10_000 # just to make sure that we don't stuck in an infinite loop
            reach₌ = Pair{Int, Int}[]
            r === nothing || for n in neighbors(g, r)
                (n == r || is_diff_edge(r, n)) && continue
                c = 1
                push!(reach₌, c => n)
            end
            if (n = length(diff_aliases)) >= 2
                as = diff_aliases[n-1]
                for (c, a) in as
                    (da = var_to_diff[a]) === nothing && continue
                    da === r && continue
                    push!(reach₌, c => da)
                end
            end
            for (c, a) in reach₌
                @info fullvars[r] => c * fullvars[a]
            end
            if r === nothing
                # TODO: updated_diff_vars check
                isempty(reach₌) && break
                dr = first(reach₌)
                var_to_diff[prev_r] = dr
                push!(updated_diff_vars, prev_r)
                prev_r = dr
            else
                prev_r = r
                r = var_to_diff[r]
            end
            for (c, v) in reach₌
                v == prev_r && continue
                dag[v] = c => prev_r
            end
            push!(diff_aliases, reach₌)
        end
        for v in zero_vars_set
            dag[v] = 0
        end
        @show nlevels
        display(diff_aliases)
        @assert length(diff_aliases) == nlevels
        @show zero_vars_set

        # clean up
        for v in dls.visited
            dls.dists[v] = typemax(Int)
            processed[v] = true
        end
        empty!(dls.visited)
        empty!(diff_aliases)
    end
    @show dag

    #=
    processed = falses(nvars)
    invag = SimpleDiGraph(nvars)
    for (v, (coeff, alias)) in pairs(ag)
        iszero(coeff) && continue
        add_edge!(invag, alias, v)
    end
    iag = InducedAliasGraph(ag, invag, var_to_diff)
    dag = AliasGraph(nvars) # alias graph for differentiated variables
    newinvag = SimpleDiGraph(nvars)
    updated_diff_vars = Int[]
    for (v, dv) in enumerate(var_to_diff)
        processed[v] && continue
        (dv === nothing && diff_to_var[v] === nothing) && continue

        r, _ = find_root!(iag, v)
           sv = fullvars[v]
           root = fullvars[r]
           @info "Found root $r" sv=>root
        level_to_var = Int[]
        extreme_var(var_to_diff, r, nothing, Val(false),
                    callback = Base.Fix1(push!, level_to_var))
        nlevels = length(level_to_var)
        current_coeff_level = Ref((0, 0))
        add_alias! = let current_coeff_level = current_coeff_level,
            level_to_var = level_to_var, dag = dag, newinvag = newinvag,
            processed = processed

            v -> begin
                coeff, level = current_coeff_level[]
                if level + 1 <= length(level_to_var)
                    av = level_to_var[level + 1]
                    if v != av # if the level_to_var isn't from the root branch
                        dag[v] = coeff => av
                        add_edge!(newinvag, av, v)

                        a = iszero(av) ? 0 : coeff * fullvars[av]
                        @info "dag $r" fullvars[v] => a
                    end
                else
                    @assert length(level_to_var) == level
                    if coeff != 1
                        dag[v] = coeff => v
                    end
                    push!(level_to_var, v)
                end
                mark_processed!(processed, var_to_diff, v)
                current_coeff_level[] = (coeff, level + 1)
            end
        end
        max_lv = 0
        clear_visited!(iag)
        Main._a[] = RootedAliasTree(iag, r)
        for (coeff, t) in children(RootedAliasTree(iag, r), Val(true))
            lv = 0
            max_lv = max(max_lv, lv)
            v = nodevalue(t)
            @info v
            iszero(v) && continue
            mark_processed!(processed, var_to_diff, v)
            v == r && continue
            if lv < length(level_to_var)
                if level_to_var[lv + 1] == v
                    continue
                end
            end
            current_coeff_level[] = coeff, lv
            extreme_var(var_to_diff, v, nothing, Val(false), callback = add_alias!)
        end
        @warn "after first"


        for (coeff, lv, t) in StatefulAliasBFS(RootedAliasTree(iag, r))
            max_lv = max(max_lv, lv)
            v = nodevalue(t)
            iszero(v) && continue
            mark_processed!(processed, var_to_diff, v)
            v == r && continue
            if lv < length(level_to_var)
                if level_to_var[lv + 1] == v
                    continue
                end
            end
            current_coeff_level[] = coeff, lv
            extreme_var(var_to_diff, v, nothing, Val(false), callback = add_alias!)
        end

        len = length(level_to_var)
        set_v_zero! = let dag = dag
            v -> dag[v] = 0
        end
        zero_av_idx = 0
        for (i, av) in enumerate(level_to_var)
            has_zero = iszero(get(ag, av, (1, 0))[1])
            for v in neighbors(newinvag, av)
                has_zero = has_zero || iszero(get(ag, v, (1, 0))[1])
            end
            if zero_av_idx == 0 && has_zero
                zero_av_idx = i
            end
        end
        # If a chain starts to equal to zero, then all its derivatives must be
        # zero. Irreducible variables are highest differentiated variables (with
        # order >= 1) that are not zero.
        if zero_av_idx > 0
            extreme_var(var_to_diff, level_to_var[zero_av_idx], nothing, Val(false),
                        callback = set_v_zero!)
        end
        # Handle virtual variables
        if nlevels < len
            for i in (nlevels + 1):len
                li = level_to_var[i]
                var_to_diff[level_to_var[i - 1]] = li
                push!(updated_diff_vars, level_to_var[i - 1])
            end
        end
    end
    =#

    for (v, (c, a)) in dag
        a = iszero(a) ? 0 : c * fullvars[a]
        @info "dag" fullvars[v] => a
    end

    # Step 4: Merge dag and ag
    removed_aliases = BitSet()
    freshag = AliasGraph(nvars)
    for (v, (c, a)) in dag
        # D(x) ~ D(y) cannot be removed if x and y are not aliases
        if v != a && !iszero(a)
            vv = v
            aa = a
            while true
                vv′ = vv
                vv = diff_to_var[vv]
                vv === nothing && break
                if !(haskey(dag, vv) && dag[vv][2] == diff_to_var[aa])
                    push!(removed_aliases, vv′)
                    @goto SKIP_FRESHAG
                end
            end
        end
        freshag[v] = c => a
        @label SKIP_FRESHAG
        push!(removed_aliases, a)
    end
    for (v, (c, a)) in ag
        (processed[v] || processed[a]) && continue
        v in removed_aliases && continue
        freshag[v] = c => a
    end
    if freshag != ag
        ag = freshag
        @show ag
        @warn "" echelon_mm
        mm = reduce!(copy(echelon_mm), ag)
        @warn "wow" mm
    end
    for (v, (c, a)) in ag
        a = iszero(a) ? 0 : c * fullvars[a]
        @info "ag" fullvars[v] => a
    end

    # Step 5: Reflect our update decisions back into the graph, and make sure
    # that the RHS of observable variables are defined.
    for (ei, e) in enumerate(mm.nzrows)
        set_neighbors!(graph, e, mm.row_cols[ei])
    end
    update_graph_neighbors!(graph, ag)
    finalag = AliasGraph(nvars)
    # RHS must still exist in the system to be valid aliases.
    needs_update = false
    for (v, (c, a)) in ag
        if iszero(a) || !isempty(𝑑neighbors(graph, a))
            finalag[v] = c => a
        else
            needs_update = true
        end
    end
    ag = finalag

    if needs_update
        mm = reduce!(copy(echelon_mm), ag)
        for (ei, e) in enumerate(mm.nzrows)
            set_neighbors!(graph, e, mm.row_cols[ei])
        end
        update_graph_neighbors!(graph, ag)
    end

    return ag, mm, updated_diff_vars
end

function update_graph_neighbors!(graph, ag)
    for eq in 1:nsrcs(graph)
        set_neighbors!(graph, eq,
                       [get(ag, n, (1, n))[2]
                        for n in 𝑠neighbors(graph, eq)
                        if !haskey(ag, n) || ag[n][2] != 0])
    end
    return graph
end

function exactdiv(a::Integer, b)
    d, r = divrem(a, b)
    @assert r == 0
    return d
end

function locally_structure_simplify!(adj_row, pivot_var, ag)
    if pivot_var === nothing
        pivot_val = nothing
    else
        pivot_val = adj_row[pivot_var]
        iszero(pivot_val) && return false
    end

    nirreducible = 0
    # When this row only as the pivot element, the pivot is zero by homogeneity
    # of the linear system.
    alias_candidate::Union{Int, Pair{Int, Int}} = 0

    # N.B.: Assumes that the non-zeros iterator is robust to modification
    # of the underlying array datastructure.
    for (var, val) in pairs(nonzerosmap(adj_row))
        # Go through every variable/coefficient in this row and apply all aliases
        # that we have so far accumulated in `ag`, updating the adj_row as
        # we go along.
        var == pivot_var && continue
        iszero(val) && continue
        alias = get(ag, var, nothing)
        if alias === nothing
            nirreducible += 1
            alias_candidate = val => var
            continue
        end
        (coeff, alias_var) = alias
        # `var = coeff * alias_var`, so we eliminate this var.
        adj_row[var] = 0
        if alias_var != 0
            # val * var = val * (coeff * alias_var) = (val * coeff) * alias_var
            val *= coeff
            # val * var + c * alias_var + ... = (val * coeff + c) * alias_var + ...
            new_coeff = (adj_row[alias_var] += val)
            if alias_var < var
                # If this adds to a coeff that was not previously accounted for,
                # and we've already passed it, make sure to count it here. We
                # need to know if there are at most 2 terms left after this
                # loop.
                #
                # We're relying on `var` being produced in sorted order here.
                nirreducible += !(alias_candidate isa Pair) ||
                                alias_var != alias_candidate[2]
                alias_candidate = new_coeff => alias_var
            end
        end
    end

    if pivot_var === nothing
        if iszero(nirreducible)
            zero!(adj_row)
        else
            dropzeros!(adj_row)
        end
        return true
    end
    # If there were only one or two terms left in the equation (including the
    # pivot variable). We can eliminate the pivot variable. Note that when
    # `nirreducible <= 1`, `alias_candidate` is uniquely determined.
    if nirreducible > 1
        dropzeros!(adj_row)
        return false
    end

    if alias_candidate isa Pair
        alias_val, alias_var = alias_candidate

        # `p` is the pivot variable, `a` is the alias variable, `v` and `c` are
        # their coefficients.
        # v * p + c * a = 0
        # v * p = -c * a
        # p = -(c / v) * a
        if iszero(alias_val)
            alias_candidate = 0
        else
            d, r = divrem(alias_val, pivot_val)
            if r == 0 && (d == 1 || d == -1)
                alias_candidate = -d => alias_var
            else
                return false
            end
        end
    end

    ag[pivot_var] = alias_candidate
    zero!(adj_row)
    return true
end

swap!(v, i, j) = v[i], v[j] = v[j], v[i]

function getcoeff(vars, coeffs, var)
    for (vj, v) in enumerate(vars)
        v == var && return coeffs[vj]
    end
    return 0
end

"""
$(SIGNATURES)

Use Kahn's algorithm to topologically sort observed equations.

Example:
```julia
julia> @variables t x(t) y(t) z(t) k(t)
(t, x(t), y(t), z(t), k(t))

julia> eqs = [
           x ~ y + z
           z ~ 2
           y ~ 2z + k
       ];

julia> ModelingToolkit.topsort_equations(eqs, [x, y, z, k])
3-element Vector{Equation}:
 Equation(z(t), 2)
 Equation(y(t), k(t) + 2z(t))
 Equation(x(t), y(t) + z(t))
```
"""
function topsort_equations(eqs, states; check = true)
    graph, assigns = observed2graph(eqs, states)
    neqs = length(eqs)
    degrees = zeros(Int, neqs)

    for 𝑠eq in 1:length(eqs)
        var = assigns[𝑠eq]
        for 𝑑eq in 𝑑neighbors(graph, var)
            # 𝑠eq => 𝑑eq
            degrees[𝑑eq] += 1
        end
    end

    q = Queue{Int}(neqs)
    for (i, d) in enumerate(degrees)
        d == 0 && enqueue!(q, i)
    end

    idx = 0
    ordered_eqs = similar(eqs, 0)
    sizehint!(ordered_eqs, neqs)
    while !isempty(q)
        𝑠eq = dequeue!(q)
        idx += 1
        push!(ordered_eqs, eqs[𝑠eq])
        var = assigns[𝑠eq]
        for 𝑑eq in 𝑑neighbors(graph, var)
            degree = degrees[𝑑eq] = degrees[𝑑eq] - 1
            degree == 0 && enqueue!(q, 𝑑eq)
        end
    end

    (check && idx != neqs) && throw(ArgumentError("The equations have at least one cycle."))

    return ordered_eqs
end

function observed2graph(eqs, states)
    graph = BipartiteGraph(length(eqs), length(states))
    v2j = Dict(states .=> 1:length(states))

    # `assigns: eq -> var`, `eq` defines `var`
    assigns = similar(eqs, Int)

    for (i, eq) in enumerate(eqs)
        lhs_j = get(v2j, eq.lhs, nothing)
        lhs_j === nothing &&
            throw(ArgumentError("The lhs $(eq.lhs) of $eq, doesn't appear in states."))
        assigns[i] = lhs_j
        vs = vars(eq.rhs)
        for v in vs
            j = get(v2j, v, nothing)
            j !== nothing && add_edge!(graph, i, j)
        end
    end

    return graph, assigns
end

function fixpoint_sub(x, dict)
    y = substitute(x, dict)
    while !isequal(x, y)
        y = x
        x = substitute(y, dict)
    end

    return x
end

function substitute_aliases(eqs, dict)
    sub = Base.Fix2(fixpoint_sub, dict)
    map(eq -> eq.lhs ~ sub(eq.rhs), eqs)
end
