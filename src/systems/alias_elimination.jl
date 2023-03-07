using SymbolicUtils: Rewriters
using Graphs.Experimental.Traversals

function alias_eliminate_graph!(state::TransformationState; kwargs...)
    mm = linear_subsys_adjmat!(state; kwargs...)
    if size(mm, 1) == 0
        ag = AliasGraph(ndsts(state.structure.graph))
        return ag, mm, ag, mm, BitSet() # No linear subsystems
    end

    @unpack graph, var_to_diff, solvable_graph = state.structure
    ag, mm, complete_ag, complete_mm = alias_eliminate_graph!(state, mm)
    if solvable_graph !== nothing
        for (ei, e) in enumerate(mm.nzrows)
            set_neighbors!(solvable_graph, e, mm.row_cols[ei])
        end
        update_graph_neighbors!(solvable_graph, ag)
    end

    return ag, mm, complete_ag, complete_mm
end

# For debug purposes
function aag_bareiss(sys::AbstractSystem)
    state = TearingState(sys)
    complete!(state.structure)
    mm = linear_subsys_adjmat!(state)
    return aag_bareiss!(state.structure.graph, state.structure.var_to_diff, mm)
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

alias_elimination(sys) = alias_elimination!(TearingState(sys))[1]
function alias_elimination!(state::TearingState; kwargs...)
    sys = state.sys
    complete!(state.structure)
    graph_orig = copy(state.structure.graph)
    ag, mm, complete_ag, complete_mm = alias_eliminate_graph!(state; kwargs...)
    isempty(ag) && return sys, ag

    fullvars = state.fullvars
    @unpack var_to_diff, graph, solvable_graph = state.structure

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
    old_to_new_eq = Vector{Int}(undef, nsrcs(graph))
    idx = 0
    cursor = 1
    ndels = length(dels)
    for i in eachindex(old_to_new_eq)
        if cursor <= ndels && i == dels[cursor]
            cursor += 1
            old_to_new_eq[i] = -1
            continue
        end
        idx += 1
        old_to_new_eq[i] = idx
    end
    n_new_eqs = idx

    old_to_new_var = Vector{Int}(undef, ndsts(graph))
    idx = 0
    for i in eachindex(old_to_new_var)
        if haskey(ag, i)
            old_to_new_var[i] = -1
        else
            idx += 1
            old_to_new_var[i] = idx
        end
    end
    n_new_vars = idx
    #for d in dels
    #    set_neighbors!(graph, d, ())
    #    set_neighbors!(solvable_graph, d, ())
    #end

    lineqs = BitSet(mm.nzrows)
    eqs_to_update = BitSet()
    nvs_orig = ndsts(graph_orig)
    for k in keys(ag)
        # We need to update `D(D(x))` when we subsitute `D(x)` as well.
        while true
            k > nvs_orig && break
            for ieq in 𝑑neighbors(graph_orig, k)
                ieq in lineqs && continue
                new_eq = old_to_new_eq[ieq]
                new_eq < 1 && continue
                push!(eqs_to_update, new_eq)
            end
            k = var_to_diff[k]
            k === nothing && break
        end
    end
    for ieq in eqs_to_update
        eq = eqs[ieq]
        eqs[ieq] = fast_substitute(eq, subs)
    end

    for old_ieq in to_expand
        ieq = old_to_new_eq[old_ieq]
        eqs[ieq] = expand_derivatives(eqs[ieq])
    end

    newstates = []
    diff_to_var = invview(var_to_diff)
    for j in eachindex(fullvars)
        if !(j in keys(ag))
            diff_to_var[j] === nothing && push!(newstates, fullvars[j])
        end
    end
    new_graph = BipartiteGraph(n_new_eqs, ndsts(graph))
    new_solvable_graph = BipartiteGraph(n_new_eqs, ndsts(graph))
    new_eq_to_diff = DiffGraph(n_new_eqs)
    eq_to_diff = state.structure.eq_to_diff
    for (i, ieq) in enumerate(old_to_new_eq)
        ieq > 0 || continue
        set_neighbors!(new_graph, ieq, 𝑠neighbors(graph, i))
        set_neighbors!(new_solvable_graph, ieq, 𝑠neighbors(solvable_graph, i))
        new_eq_to_diff[ieq] = eq_to_diff[i]
    end
    # update DiffGraph
    new_var_to_diff = DiffGraph(length(var_to_diff))
    for v in 1:length(var_to_diff)
        (haskey(ag, v)) && continue
        new_var_to_diff[v] = var_to_diff[v]
    end
    state.structure.graph = new_graph
    state.structure.solvable_graph = new_solvable_graph
    state.structure.eq_to_diff = new_eq_to_diff
    state.structure.var_to_diff = new_var_to_diff

    #=
    new_graph = BipartiteGraph(n_new_eqs, n_new_vars)
    new_solvable_graph = BipartiteGraph(n_new_eqs, n_new_vars)
    new_eq_to_diff = DiffGraph(n_new_eqs)
    eq_to_diff = state.structure.eq_to_diff
    new_var_to_diff = DiffGraph(n_new_vars)
    var_to_diff = state.structure.var_to_diff
    for (i, ieq) in enumerate(old_to_new_eq)
        ieq > 0 || continue
        set_neighbors!(new_graph, ieq, [old_to_new_var[v] for v in 𝑠neighbors(graph, i) if old_to_new_var[v] > 0])
        set_neighbors!(new_solvable_graph, ieq, [old_to_new_var[v] for v in 𝑠neighbors(solvable_graph, i) if old_to_new_var[v] > 0])
        new_eq_to_diff[ieq] = eq_to_diff[i]
    end
    new_fullvars = Vector{Any}(undef, n_new_vars)
    for (i, iv) in enumerate(old_to_new_var)
        iv > 0 || continue
        new_var_to_diff[iv] = var_to_diff[i]
        new_fullvars[iv] = fullvars[i]
    end
    state.structure.graph = new_graph
    state.structure.solvable_graph = new_solvable_graph
    state.structure.eq_to_diff = complete(new_eq_to_diff)
    state.structure.var_to_diff = complete(new_var_to_diff)
    state.fullvars = new_fullvars
    =#

    sys = state.sys
    @set! sys.eqs = eqs
    @set! sys.states = newstates
    @set! sys.observed = [observed(sys); obs]
    state.sys = sys
    return invalidate_cache!(sys), ag
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
                if (mask === nothing || mask[v])
                    return (CartesianIndex(i, v), M.row_vals[i][j])
                end
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
    if i > length(ag.aliasto)
        resize!(ag.aliasto, i)
    end
    if ag.aliasto[i] === nothing
        push!(ag.eliminated, i)
    end
    ag.aliasto[i] = 0
    return 0 => 0
end

function Base.setindex!(ag::AliasGraph,
                        p::Union{Pair{<:Integer, Int}, Tuple{<:Integer, Int}},
                        i::Integer)
    (c, v) = p
    if c == 0 || v == 0
        ag[i] = 0
        return p
    end
    @assert v != 0 && c in (-1, 1)
    if i > length(ag.aliasto)
        resize!(ag.aliasto, i)
    end
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

canonicalize(a, b) = a <= b ? (a, b) : (b, a)
struct WeightedGraph{T, W} <: AbstractGraph{T}
    graph::SimpleGraph{T}
    dict::Dict{Tuple{T, T}, W}
end
function WeightedGraph{T, W}(n) where {T, W}
    WeightedGraph{T, W}(SimpleGraph{T}(n), Dict{Tuple{T, T}, W}())
end

function Graphs.add_edge!(g::WeightedGraph, u, v, w)
    r = add_edge!(g.graph, u, v)
    r && (g.dict[canonicalize(u, v)] = w)
    r
end
Graphs.add_vertex!(g::WeightedGraph) = add_vertex!(g.graph)
Graphs.has_edge(g::WeightedGraph, u, v) = has_edge(g.graph, u, v)
Graphs.ne(g::WeightedGraph) = ne(g.graph)
Graphs.nv(g::WeightedGraph) = nv(g.graph)
get_weight(g::WeightedGraph, u, v) = g.dict[canonicalize(u, v)]
Graphs.is_directed(::Type{<:WeightedGraph}) = false
Graphs.inneighbors(g::WeightedGraph, v) = inneighbors(g.graph, v)
Graphs.outneighbors(g::WeightedGraph, v) = outneighbors(g.graph, v)
Graphs.vertices(g::WeightedGraph) = vertices(g.graph)
Graphs.edges(g::WeightedGraph) = vertices(g.graph)

function equality_diff_graph(ag::AliasGraph, var_to_diff::DiffGraph)
    g = SimpleDiGraph{Int}(length(var_to_diff))
    eqg = WeightedGraph{Int, Int}(length(var_to_diff))
    zero_vars = Int[]
    for (v, (c, a)) in ag
        if iszero(a)
            push!(zero_vars, v)
            continue
        end
        add_edge!(g, v, a)
        add_edge!(g, a, v)

        add_edge!(eqg, v, a, c)
    end
    transitiveclosure!(g)
    weighted_transitiveclosure!(eqg)

    for (v, dv) in enumerate(var_to_diff)
        dv isa Int || continue
        add_edge!(g, v, dv)
        add_edge!(g, dv, v)
    end
    g, eqg, zero_vars
end

function weighted_transitiveclosure!(g)
    cps = connected_components(g)
    for cp in cps
        n = length(cp)
        for k in cp
            for i′ in 1:n, j′ in (i′ + 1):n
                i = cp[i′]
                j = cp[j′]
                (has_edge(g, i, k) && has_edge(g, k, j)) || continue
                add_edge!(g, i, j, get_weight(g, i, k) * get_weight(g, k, j))
            end
        end
    end
    return g
end

struct DiffLevelState <: Traversals.AbstractTraversalState
    dists::Vector{Int}
    var_to_diff::DiffGraph
    visited::BitSet
end

function DiffLevelState(g::SimpleDiGraph, var_to_diff)
    DiffLevelState(fill(typemax(Int), nv(g)), var_to_diff, BitSet())
end

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

count_nonzeros(a::AbstractArray) = count(!iszero, a)

# N.B.: Ordinarily sparse vectors allow zero stored elements.
# Here we have a guarantee that they won't, so we can make this identification
count_nonzeros(a::CLILVector) = nnz(a)

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

function aag_bareiss!(graph, var_to_diff, mm_orig::SparseMatrixCLIL{T, Ti}) where {T, Ti}
    mm = copy(mm_orig)
    linear_equations_set = BitSet(mm_orig.nzrows)

    # All unassigned (not a pivot) algebraic variables that only appears in
    # linear algebraic equations can be set to 0.
    #
    # For all the other variables, we can update the original system with
    # Bareiss'ed coefficients as Gaussian elimination is nullspace perserving
    # and we are only working on linear homogeneous subsystem.

    is_algebraic = let var_to_diff = var_to_diff
        v -> var_to_diff[v] === nothing === invview(var_to_diff)[v]
    end
    is_linear_variables = is_algebraic.(1:length(var_to_diff))
    for i in 𝑠vertices(graph)
        # only consider linear algebraic equations
        (i in linear_equations_set && all(is_algebraic, 𝑠neighbors(graph, i))) &&
            continue
        for j in 𝑠neighbors(graph, i)
            is_linear_variables[j] = false
        end
    end
    solvable_variables = findall(is_linear_variables)

    local bar
    try
        bar = do_bareiss!(mm, mm_orig, is_linear_variables)
    catch e
        e isa OverflowError || rethrow(e)
        mm = convert(SparseMatrixCLIL{BigInt, Ti}, mm_orig)
        bar = do_bareiss!(mm, mm_orig, is_linear_variables)
    end

    return mm, solvable_variables, bar
end

function do_bareiss!(M, Mold, is_linear_variables)
    rank1r = Ref{Union{Nothing, Int}}(nothing)
    find_pivot = let rank1r = rank1r
        (M, k) -> begin
            if rank1r[] === nothing
                r = find_masked_pivot(is_linear_variables, M, k)
                r !== nothing && return r
                rank1r[] = k - 1
            end
            # TODO: It would be better to sort the variables by
            # derivative order here to enable more elimination
            # opportunities.
            return find_masked_pivot(nothing, M, k)
        end
    end
    pivots = Int[]
    find_and_record_pivot = let pivots = pivots
        (M, k) -> begin
            r = find_pivot(M, k)
            r === nothing && return nothing
            push!(pivots, r[1][2])
            return r
        end
    end
    myswaprows! = let Mold = Mold
        (M, i, j) -> begin
            Mold !== nothing && swaprows!(Mold, i, j)
            swaprows!(M, i, j)
        end
    end
    bareiss_ops = ((M, i, j) -> nothing, myswaprows!,
                   bareiss_update_virtual_colswap_mtk!, bareiss_zero!)

    rank2, = bareiss!(M, bareiss_ops; find_pivot = find_and_record_pivot)
    rank1 = something(rank1r[], rank2)
    (rank1, rank2, pivots)
end

# Kind of like the backward substitution, but we don't actually rely on it
# being lower triangular. We eliminate a variable if there are at most 2
# variables left after the substitution.
function lss(mm, ag, pivots)
    ei -> let mm = mm, pivots = pivots, ag = ag
        vi = pivots === nothing ? nothing : pivots[ei]
        locally_structure_simplify!((@view mm[ei, :]), vi, ag)
    end
end

function reduce!(mm, mm_orig, ag, rank2, pivots = nothing)
    lss! = lss(mm, ag, pivots)
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

    return mm
end

function simple_aliases!(ag, graph, var_to_diff, mm_orig)
    echelon_mm, solvable_variables, (rank1, rank2, pivots) = aag_bareiss!(graph,
                                                                          var_to_diff,
                                                                          mm_orig)

    # Step 2: Simplify the system using the Bareiss factorization
    rk1vars = BitSet(@view pivots[1:rank1])
    for v in solvable_variables
        v in rk1vars && continue
        ag[v] = 0
    end

    mm = reduce!(copy(echelon_mm), mm_orig, ag, rank2, pivots)
    return mm, echelon_mm
end

function var_derivative_here!(state, processed, g, eqg, dls, diff_var)
    newvar = var_derivative!(state, diff_var)
    @assert newvar == length(processed) + 1
    push!(processed, true)
    add_vertex!(g)
    add_vertex!(eqg)
    add_edge!(g, diff_var, newvar)
    add_edge!(g, newvar, diff_var)
    push!(dls.dists, typemax(Int))
    return newvar
end

function alias_eliminate_graph!(state::TransformationState, mm_orig::SparseMatrixCLIL)
    @unpack graph, var_to_diff = state.structure
    # Step 1: Perform Bareiss factorization on the adjacency matrix of the linear
    #         subsystem of the system we're interested in.
    #
    nvars = ndsts(graph)
    ag = AliasGraph(nvars)
    complete_ag = AliasGraph(nvars)
    mm, echelon_mm = simple_aliases!(ag, graph, var_to_diff, mm_orig)

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
    diff_to_var = invview(var_to_diff)
    processed = falses(nvars)
    g, eqg, zero_vars = equality_diff_graph(ag, var_to_diff)
    dls = DiffLevelState(g, var_to_diff)
    original_nvars = length(var_to_diff)

    is_diff_edge = let var_to_diff = var_to_diff
        (v, w) -> var_to_diff[v] == w || var_to_diff[w] == v
    end
    diff_aliases = Vector{Pair{Int, Int}}[]
    stems = Vector{Int}[]
    stem_set = BitSet()
    for (v, dv) in enumerate(var_to_diff)
        processed[v] && continue
        (dv === nothing && diff_to_var[v] === nothing) && continue
        stem = Int[]
        r = find_root!(dls, g, v)
        prev_r = -1
        for _ in 1:10_000 # just to make sure that we don't stuck in an infinite loop
            reach₌ = Pair{Int, Int}[]
            # `r` is aliased to its equality aliases
            if r !== nothing
                for n in neighbors(eqg, r)
                    (n == r || is_diff_edge(r, n)) && continue
                    c = get_weight(eqg, r, n)
                    push!(reach₌, c => n)
                end
            end
            # `r` is aliased to its previous differentiation level's aliases'
            # derivative
            if (n = length(diff_aliases)) >= 1
                as = diff_aliases[n]
                for (c, a) in as
                    (da = var_to_diff[a]) === nothing && continue
                    da === r && continue
                    push!(reach₌, c => da)
                    # `r` is aliased to its previous differentiation level's
                    # aliases' derivative's equality aliases
                    r === nothing || for n in neighbors(eqg, da)
                        (n == da || n == prev_r || is_diff_edge(prev_r, n)) && continue
                        c′ = get_weight(eqg, da, n)
                        push!(reach₌, c * c′ => n)
                    end
                end
            end
            if r === nothing
                isempty(reach₌) && break
                let stem_set = stem_set
                    any(x -> x[2] in stem_set, reach₌) && break
                end
                # See example in the box above where D(D(D(z))) doesn't appear
                # in the original system and needs to added, so we can alias to it.
                # We do that here.
                @assert prev_r !== -1
                prev_r = var_derivative_here!(state, processed, g, eqg, dls, prev_r)
                r = nothing
            else
                prev_r = r
                r = var_to_diff[r]
            end
            prev_r in stem_set && break
            push!(stem_set, prev_r)
            push!(stem, prev_r)
            push!(diff_aliases, reach₌)
            for (c, v) in reach₌
                v == prev_r && continue
                add_edge!(eqg, v, prev_r, c)
            end
        end

        @assert length(stem) == length(diff_aliases)
        for i in eachindex(stem)
            a = stem[i]
            for (c, v) in diff_aliases[i]
                # alias edges that coincide with diff edges are handled later
                v in stem_set && continue
                dag[v] = c => a
            end
        end
        push!(stems, stem)

        # clean up
        for v in dls.visited
            dls.dists[v] = typemax(Int)
            processed[v] = true
        end
        empty!(dls.visited)
        empty!(diff_aliases)
        empty!(stem_set)
    end

    # Obtain transitive closure after completing the alias edges from diff
    # edges. As a performance optimization, we only compute the transitive
    # closure once at the very end.
    weighted_transitiveclosure!(eqg)
    zero_vars_set = BitSet()
    for stem in stems
        # Canonicalize by preferring the lower differentiated variable
        # If we have the system
        # ```
        # D(x) ~ x
        # D(x) + D(y) ~ 0
        # ```
        # preferring the lower variable would lead to
        # ```
        # D(x) ~ x      <== added back because `x := D(x)` removes `D(x)`
        # D(y) ~ -x
        # ```
        # while preferring the higher variable would lead to
        # ```
        # D(x) + D(y) ~ 0
        # ```
        # which is not correct.
        for i in 1:(length(stem) - 1)
            r = stem[i]
            for dr in @view stem[(i + 1):end]
                # We cannot reduce newly introduced variables like `D(D(D(z)))`
                # in the example box above.
                dr > original_nvars && continue
                if has_edge(eqg, r, dr)
                    c = get_weight(eqg, r, dr)
                    dag[dr] = c => r
                end
            end
        end
        # If a non-differentiated variable equals to 0, then we can eliminate
        # the whole differentiation chain. Otherwise, we will still have to keep
        # the lowest differentiated variable in the differentiation chain.
        # E.g.
        # ```
        # D(x) ~ 0
        # D(D(x)) ~ y
        # ```
        # reduces to
        # ```
        # D(x) ~ 0
        # y := 0
        # ```
        # but
        # ```
        # x ~ 0
        # D(x) ~ y
        # ```
        # reduces to
        # ```
        # x := 0
        # y := 0
        # ```
        for v in zero_vars
            for a in Iterators.flatten((v, outneighbors(eqg, v)))
                while true
                    push!(zero_vars_set, a)
                    a = var_to_diff[a]
                    a === nothing && break
                end
            end
        end
        for v in zero_vars_set
            while (iv = diff_to_var[v]) in zero_vars_set
                v = iv
            end
            complete_ag[v] = 0
            if diff_to_var[v] === nothing # `v` is reducible
                dag[v] = 0
            end
            # reducible after v
            while (v = var_to_diff[v]) !== nothing
                complete_ag[v] = 0
                dag[v] = 0
            end
        end
        empty!(zero_vars_set)
    end

    # update `dag`
    for k in keys(dag)
        dag[k]
    end

    # Step 4: Merge dag and ag
    removed_aliases = BitSet()
    merged_ag = AliasGraph(nvars)
    for (v, (c, a)) in dag
        complete_ag[v] = c => a
    end
    for (v, (c, a)) in ag
        (processed[v] || (!iszero(a) && processed[a])) && continue
        complete_ag[v] = c => a
    end
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
                    @goto SKIP_merged_ag
                end
            end
        end
        merged_ag[v] = c => a
        @label SKIP_merged_ag
        push!(removed_aliases, a)
    end
    for (v, (c, a)) in ag
        (processed[v] || (!iszero(a) && processed[a])) && continue
        v in removed_aliases && continue
        merged_ag[v] = c => a
    end
    ag = merged_ag
    @set! echelon_mm.ncols = length(var_to_diff)
    @set! mm_orig.ncols = length(var_to_diff)
    mm = reduce!(copy(echelon_mm), mm_orig, ag, size(echelon_mm, 1))

    # Step 5: Reflect our update decisions back into the graph, and make sure
    # that the RHS of observable variables are defined.
    for (ei, e) in enumerate(mm.nzrows)
        set_neighbors!(graph, e, mm.row_cols[ei])
    end
    update_graph_neighbors!(graph, ag)
    finalag = AliasGraph(nvars)
    # RHS or its derivaitves must still exist in the system to be valid aliases.
    needs_update = false
    function contains_v_or_dv(var_to_diff, graph, v)
        counter = 0
        while true
            isempty(𝑑neighbors(graph, v)) || return true
            v = var_to_diff[v]
            v === nothing && return false
            counter += 1
            counter > 10_000 &&
                error("Internal error: there's an infinite loop in the `var_to_diff` graph.")
        end
    end
    for (v, (c, a)) in ag
        if iszero(a) || contains_v_or_dv(var_to_diff, graph, a)
            finalag[v] = c => a
        else
            needs_update = true
        end
    end
    ag = finalag

    if needs_update
        mm = reduce!(copy(echelon_mm), mm_orig, ag, size(echelon_mm, 1))
    end
    # applying new `ag` to `mm` might lead to linear dependence, so we have to
    # re-run Bareiss.
    mm, = aag_bareiss!(graph, var_to_diff, mm)
    for (ei, e) in enumerate(mm.nzrows)
        set_neighbors!(graph, e, mm.row_cols[ei])
    end
    update_graph_neighbors!(graph, ag)

    complete_mm = reduce!(copy(echelon_mm), mm_orig, complete_ag, size(echelon_mm, 1))
    return ag, mm, complete_ag, complete_mm
end

function update_graph_neighbors!(graph, ag)
    for eq in 1:nsrcs(graph)
        set_neighbors!(graph, eq,
                       Int[get(ag, n, (1, n))[2]
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
    # If `pivot_var === nothing`, then we only apply `ag` to `adj_row`
    if pivot_var === nothing
        pivot_val = nothing
    else
        pivot_val = adj_row[pivot_var]
        iszero(pivot_val) && return false
    end

    nirreducible = 0
    # When this row only as the pivot element, the pivot is zero by homogeneity
    # of the linear system.
    alias_candidate::Union{Int, Pair{eltype(adj_row), Int}} = 0

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
        return false
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
