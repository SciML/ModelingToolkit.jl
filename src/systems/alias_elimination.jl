using SymbolicUtils: Rewriters

const KEEP = typemin(Int)

function alias_eliminate_graph!(state::TransformationState)
    mm = linear_subsys_adjmat(state)
    size(mm, 1) == 0 && return nothing, mm, BitSet() # No linear subsystems

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
    sys = state.sys
    ag, mm, updated_diff_vars = alias_eliminate_graph!(state)
    ag === nothing && return sys

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
    for (v, (coeff, alias)) in pairs(ag)
        subs[fullvars[v]] = iszero(coeff) ? 0 : coeff * fullvars[alias]
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

    for (ieq, eq) in enumerate(eqs)
        eqs[ieq] = substitute(eq, subs)
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
    @set! sys.observed = [observed(sys); [lhs ~ rhs for (lhs, rhs) in pairs(subs)]]
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
    if var in keys(ag)
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

function Base.setindex!(ag::AliasGraph, p::Pair{Int, Int}, i::Integer)
    (c, v) = p
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
    dels = Int[]
    for (i, rs) in enumerate(mm.row_cols)
        rvals = mm.row_vals[i]
        j = 1
        while j <= length(rs)
            c = rs[j]
            _alias = get(ag, c, nothing)
            if _alias !== nothing
                push!(dels, j)
                coeff, alias = _alias
                iszero(coeff) && (j += 1; continue)
                inc = coeff * rvals[j]
                i = searchsortedfirst(rs, alias)
                if i > length(rs) || rs[i] != alias
                    # if we add a variable to what we already visited, make sure
                    # to bump the cursor.
                    j += i <= j
                    for (i, e) in enumerate(dels)
                        e >= i && (dels[i] += 1)
                    end
                    insert!(rs, i, alias)
                    insert!(rvals, i, inc)
                else
                    rvals[i] += inc
                end
            end
            j += 1
        end
        deleteat!(rs, dels)
        deleteat!(rvals, dels)
        empty!(dels)
        for (j, v) in enumerate(rvals)
            iszero(v) && push!(dels, j)
        end
        deleteat!(rs, dels)
        deleteat!(rvals, dels)
        empty!(dels)
    end
    mm
end

struct InducedAliasGraph
    ag::AliasGraph
    invag::SimpleDiGraph{Int}
    var_to_diff::DiffGraph
    visited::BitSet
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
AbstractTrees.children(rat::RootedAliasTree) = RootedAliasChildren(rat)
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

struct RootedAliasChildren
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

@inline function _iterate(c::RootedAliasChildren, s = nothing)
    rat = c.t
    @unpack iag, root = rat
    @unpack ag, invag, var_to_diff = iag
    (iszero(root) || (root = var_to_diff[root]) === nothing) && return nothing
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
    return (ag[e][1], RootedAliasTree(iag, e)), (stage, iterate(it, ns))
end

count_nonzeros(a::AbstractArray) = count(!iszero, a)

# N.B.: Ordinarily sparse vectors allow zero stored elements.
# Here we have a guarantee that they won't, so we can make this identification
count_nonzeros(a::SparseVector) = nnz(a)

function aag_bareiss!(graph, var_to_diff, mm_orig::SparseMatrixCLIL, irreducibles = ())
    mm = copy(mm_orig)
    is_linear_equations = falses(size(AsSubMatrix(mm_orig), 1))
    for e in mm_orig.nzrows
        is_linear_equations[e] = true
    end

    # If linear highest differentiated variables cannot be assigned to a pivot,
    # then we can set it to zero. We use `rank1` to track this.
    #
    # We only use alias graph to record reducible variables. We use `rank2` to
    # track this.
    #
    # For all the other variables, we can update the original system with
    # Bareiss'ed coefficients as Gaussian elimination is nullspace perserving
    # and we are only working on linear homogeneous subsystem.
    is_linear_variables = isnothing.(var_to_diff)
    is_reducible = trues(length(var_to_diff))
    for v in irreducibles
        is_linear_variables[v] = false
        is_reducible[v] = false
    end
    for i in 𝑠vertices(graph)
        is_linear_equations[i] && continue
        for j in 𝑠neighbors(graph, i)
            is_linear_variables[j] = false
        end
    end
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
    for v in setdiff(solvable_variables, @view pivots[1:rank1])
        ag[v] = 0
    end

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

    return mm
end

function alias_eliminate_graph!(graph, var_to_diff, mm_orig::SparseMatrixCLIL)
    # Step 1: Perform bareiss factorization on the adjacency matrix of the linear
    #         subsystem of the system we're interested in.
    #
    nvars = ndsts(graph)
    ag = AliasGraph(nvars)
    mm = simple_aliases!(ag, graph, var_to_diff, mm_orig)

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
    # the system.
    #
    # To finish the algorithm, we backtrack to the root differentiation chain.
    # If the variable already exists in the chain, then we alias them
    # (e.g. `x_t ⇒ D(D(z))`), else, we substitute and update `var_to_diff`.
    #
    # Note that since we always prefer the higher differentiated variable and
    # with a tie breaking strategy. The root variable (in this case `z`) is
    # always uniquely determined. Thus, the result is well-defined.
    diff_to_var = invview(var_to_diff)
    invag = SimpleDiGraph(nvars)
    for (v, (coeff, alias)) in pairs(ag)
        iszero(coeff) && continue
        add_edge!(invag, alias, v)
    end
    processed = falses(nvars)
    iag = InducedAliasGraph(ag, invag, var_to_diff)
    newag = AliasGraph(nvars)
    irreducibles = BitSet()
    updated_diff_vars = Int[]
    for (v, dv) in enumerate(var_to_diff)
        processed[v] && continue
        (dv === nothing && diff_to_var[v] === nothing) && continue

        r, _ = find_root!(iag, v)
        #   sv = fullvars[v]
        #   root = fullvars[r]
        #   @info "Found root $r" sv=>root
        level_to_var = Int[]
        extreme_var(var_to_diff, r, nothing, Val(false),
                    callback = Base.Fix1(push!, level_to_var))
        nlevels = length(level_to_var)
        current_coeff_level = Ref((0, 0))
        add_alias! = let current_coeff_level = current_coeff_level,
            level_to_var = level_to_var, newag = newag, processed = processed

            v -> begin
                coeff, level = current_coeff_level[]
                if level + 1 <= length(level_to_var)
                    av = level_to_var[level + 1]
                    if v != av # if the level_to_var isn't from the root branch
                        newag[v] = coeff => av
                    end
                else
                    @assert length(level_to_var) == level
                    push!(level_to_var, v)
                end
                processed[v] = true
                current_coeff_level[] = (coeff, level + 1)
            end
        end
        max_lv = 0
        clear_visited!(iag)
        for (coeff, lv, t) in StatefulAliasBFS(RootedAliasTree(iag, r))
            max_lv = max(max_lv, lv)
            v = nodevalue(t)
            iszero(v) && continue
            processed[v] = true
            v == r && continue
            if lv < length(level_to_var)
                if level_to_var[lv + 1] == v
                    continue
                end
            end
            current_coeff_level[] = coeff, lv
            extreme_var(var_to_diff, v, nothing, Val(false), callback = add_alias!)
        end
        max_lv > 0 || continue

        set_v_zero! = let newag = newag
            v -> newag[v] = 0
        end
        for (i, v) in enumerate(level_to_var)
            _alias = get(ag, v, nothing)
            push!(irreducibles, v)
            if _alias !== nothing && iszero(_alias[1]) && i < length(level_to_var)
                # we have `x = 0`
                v = level_to_var[i + 1]
                extreme_var(var_to_diff, v, nothing, Val(false), callback = set_v_zero!)
                break
            end
        end
        if nlevels < (new_nlevels = length(level_to_var))
            for i in (nlevels + 1):new_nlevels
                var_to_diff[level_to_var[i - 1]] = level_to_var[i]
                push!(updated_diff_vars, level_to_var[i - 1])
            end
        end
    end

    if !isempty(irreducibles)
        ag = newag
        mm_orig2 = isempty(ag) ? mm_orig : reduce!(copy(mm_orig), ag)
        mm = simple_aliases!(ag, graph, var_to_diff, mm_orig2, irreducibles)
    end

    # for (v, (c, a)) in ag
    #     va = iszero(a) ? a : fullvars[a]
    #     @info "new alias" fullvars[v]=>(c, va)
    # end

    # Step 4: Reflect our update decisions back into the graph
    for (ei, e) in enumerate(mm.nzrows)
        set_neighbors!(graph, e, mm.row_cols[ei])
    end

    # because of `irreducibles`, `mm` cannot always be trusted.
    return ag, mm, updated_diff_vars
end

function exactdiv(a::Integer, b)
    d, r = divrem(a, b)
    @assert r == 0
    return d
end

function locally_structure_simplify!(adj_row, pivot_var, ag)
    pivot_val = adj_row[pivot_var]
    iszero(pivot_val) && return false

    nirreducible = 0
    alias_candidate::Pair{Int, Int} = 0 => 0

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

    # If there were only one or two terms left in the equation (including the
    # pivot variable). We can eliminate the pivot variable. Note that when
    # `nirreducible <= 1`, `alias_candidate` is uniquely determined.
    nirreducible <= 1 || return false

    if alias_candidate isa Pair
        alias_val, alias_var = alias_candidate

        # `p` is the pivot variable, `a` is the alias variable, `v` and `c` are
        # their coefficients.
        # v * p + c * a = 0
        # v * p = -c * a
        # p = -(c / v) * a
        d, r = divrem(alias_val, pivot_val)
        if r == 0 && (d == 1 || d == -1)
            alias_candidate = -d => alias_var
        else
            return false
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
