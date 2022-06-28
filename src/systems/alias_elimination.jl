using SymbolicUtils: Rewriters

const KEEP = typemin(Int)

function alias_eliminate_graph!(state::TransformationState)
    mm = linear_subsys_adjmat(state)
    size(mm, 1) == 0 && return nothing, mm # No linear subsystems

    @unpack graph, var_to_diff = state.structure

    ag, mm = alias_eliminate_graph!(graph, complete(var_to_diff), mm)
    return ag, mm
end

# For debug purposes
function aag_bareiss(sys::AbstractSystem)
    state = TearingState(sys)
    mm = linear_subsys_adjmat(state)
    return aag_bareiss!(state.structure.graph, complete(state.structure.var_to_diff), mm)
end

function walk_to_root(ag, var_to_diff, v::Integer)
    diff_to_var = invview(var_to_diff)

    v′::Union{Nothing, Int} = v
    @label HAS_BRANCH
    while (v′ = diff_to_var[v]) !== nothing
        v = v′
    end
    # `v` is now not differentiated in the current chain.
    # Now we recursively walk to root variable's chain.
    while true
        next_v = get(ag, v, nothing)
        next_v === nothing || (v = next_v[2]; @goto HAS_BRANCH)
        (v′ = var_to_diff[v]) === nothing && break
        v = v′
    end

    # Descend to the root from the chain
    while (v′ = diff_to_var[v]) !== nothing
        v = v′
    end
    v
end

function visit_differential_aliases!(ag, level_to_var, processed, invag, var_to_diff, v, level=0)
    processed[v] && return nothing
    for n in neighbors(invag, v)
        # TODO: we currently only handle `coeff == 1`
        if isone(ag[n][1])
            visit_differential_aliases!(ag, level_to_var, processed, invag, var_to_diff, n, level)
        end
    end
    # Note that we don't need to update `invag`
    if 1 <= level + 1 <= length(level_to_var)
        root_var = level_to_var[level + 1]
        if v != root_var
            ag[v] = 1 => root_var
        end
    else
        @assert length(level_to_var) == level
        push!(level_to_var, v)
    end
    processed[v] = true
    if (dv = var_to_diff[v]) !== nothing
        visit_differential_aliases!(ag, level_to_var, processed, invag, var_to_diff, dv, level + 1)
    end
    if (iv = invview(var_to_diff)[v]) !== nothing
        visit_differential_aliases!(ag, level_to_var, processed, invag, var_to_diff, iv, level - 1)
    end
    return nothing
end

function alias_elimination(sys)
    state = TearingState(sys; quick_cancel = true)
    ag, mm = alias_eliminate_graph!(state)
    ag === nothing && return sys

    fullvars = state.fullvars
    @unpack var_to_diff, graph = state.structure

    # After `alias_eliminate_graph!`, `var_to_diff` and `ag` form a tree
    # structure like the following:
    #
    #         x   -->   D(x)
    #         ⇓          ⇑
    #         ⇓         x_t   -->   D(x_t)
    #         ⇓               |---------------|
    # z --> D(z)  --> D(D(z)) |--> D(D(D(z))) |
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
    D = has_iv(sys) ? Differential(get_iv(sys)) : nothing
    diff_to_var = invview(var_to_diff)
    invag = SimpleDiGraph(length(fullvars))
    for (v, (coeff, alias)) in pairs(ag)
        iszero(coeff) && continue
        add_edge!(invag, alias, v)
    end
    processed = falses(length(var_to_diff))
    for (v, dv) in enumerate(var_to_diff)
        processed[v] && continue
        (dv === nothing && diff_to_var[v] === nothing) && continue

        r = walk_to_root(ag, var_to_diff, v)
        level_to_var = Int[r]
        v′′::Union{Nothing, Int} = v′::Int = r
        while (v′′ = var_to_diff[v′]) !== nothing
            v′ = v′′
            push!(level_to_var, v′)
        end
        nlevels = length(level_to_var)
        visit_differential_aliases!(ag, level_to_var, processed, invag, var_to_diff, r)
        if nlevels < (new_nlevels = length(level_to_var))
            @assert !(D isa Nothing)
            for i in (nlevels + 1):new_nlevels
                var_to_diff[level_to_var[i-1]] = level_to_var[i]
                fullvars[level_to_var[i]] = D(fullvars[level_to_var[i - 1]])
            end
        end
    end

    subs = Dict()
    for (v, (coeff, alias)) in pairs(ag)
        subs[fullvars[v]] = iszero(coeff) ? 0 : coeff * fullvars[alias]
    end

    dels = Set{Int}()
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
    dels = sort(collect(dels))
    deleteat!(eqs, dels)

    for (ieq, eq) in enumerate(eqs)
        eqs[ieq] = substitute(eq, subs)
    end

    newstates = []
   for j in eachindex(fullvars)
        if j in keys(ag)
            _, var = ag[j]
            iszero(var) && continue
            # Put back equations for alias eliminated dervars
            if isdervar(state.structure, var)
                has_higher_order = false
                v = var
                while (v = var_to_diff[v]) !== nothing
                    if !(v in keys(ag))
                        has_higher_order = true
                        break
                    end
                end
                if !has_higher_order
                    rhs = fullvars[j]
                    push!(eqs, subs[fullvars[j]] ~ rhs)
                    diff_to_var[j] === nothing && push!(newstates, rhs)
                end
            end
        else
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
Base.length(agk::AliasGraphKeySet) = Base.length(agk.ag.eliminated)
function Base.in(i::Int, agk::AliasGraphKeySet)
    aliasto = agk.ag.aliasto
    1 <= i <= length(aliasto) && aliasto[i] !== nothing
end

count_nonzeros(a::AbstractArray) = count(!iszero, a)

# N.B.: Ordinarily sparse vectors allow zero stored elements.
# Here we have a guarantee that they won't, so we can make this identification
count_nonzeros(a::SparseVector) = nnz(a)

function aag_bareiss!(graph, var_to_diff, mm_orig::SparseMatrixCLIL)
    mm = copy(mm_orig)
    is_linear_equations = falses(size(AsSubMatrix(mm_orig), 1))
    for e in mm_orig.nzrows
        is_linear_equations[e] = true
    end

    # For now, only consider variables linear that are not differentiated.
    # We could potentially apply the same logic to variables whose derivative
    # is also linear, but that's a TODO.
    diff_to_var = invview(var_to_diff)
    is_linear_variables = .&(isnothing.(var_to_diff), isnothing.(diff_to_var))
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
        rank2, = bareiss!(M, bareiss_ops; find_pivot = find_and_record_pivot)
        rank1 = something(rank1, rank2)
        (rank1, rank2, pivots)
    end

    return mm, solvable_variables, do_bareiss!(mm, mm_orig)
end

function alias_eliminate_graph!(graph, var_to_diff, mm_orig::SparseMatrixCLIL)
    # Step 1: Perform bareiss factorization on the adjacency matrix of the linear
    #         subsystem of the system we're interested in.
    #
    # Let `m = the number of linear equations` and `n = the number of
    # variables`.
    #
    # `do_bareiss` conceptually gives us this system:
    # rank1 | [ M₁₁  M₁₂ | M₁₃ ]   [v₁] = [0]
    # rank2 | [ 0    M₂₂ | M₂₃ ] P [v₂] = [0]
    # -------------------|------------------------
    #         [ 0    0   | 0   ]   [v₃] = [0]
    #
    # Where `v₁` are the purely linear variables (i.e. those that only appear in linear equations),
    # `v₂` are the variables that may be potentially solved by the linear system and v₃ are the variables
    # that contribute to the equations, but are not solved by the linear system. Note
    # that the complete system may be larger than the linear subsystem and include variables
    # that do not appear here.
    mm, solvable_variables, (rank1, rank2, pivots) = aag_bareiss!(graph, var_to_diff,
                                                                  mm_orig)

    # Step 2: Simplify the system using the Bareiss factorization

    ag = AliasGraph(size(mm, 2))

    # First, eliminate variables that only appear in linear equations and were removed
    # completely from the coefficient matrix. These are technically singularities in
    # the matrix, but assigning them to 0 is a feasible assignment and works well in
    # practice.
    for v in setdiff(solvable_variables, @view pivots[1:rank1])
        ag[v] = 0
    end

    # Kind of like the backward substitution, but we don't actually rely on it
    # being lower triangular. We eliminate a variable if there are at most 2
    # variables left after the substitution.
    diff_to_var = invview(var_to_diff)
    function lss!(ei::Integer)
        vi = pivots[ei]
        locally_structure_simplify!((@view mm[ei, :]), vi, ag, var_to_diff)
    end

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

    # Step 3: Reflect our update decisions back into the graph
    for (ei, e) in enumerate(mm.nzrows)
        set_neighbors!(graph, e, mm.row_cols[ei])
    end

    return ag, mm
end

function exactdiv(a::Integer, b)
    d, r = divrem(a, b)
    @assert r == 0
    return d
end

function locally_structure_simplify!(adj_row, pivot_var, ag, var_to_diff)
    pivot_val = adj_row[pivot_var]
    iszero(pivot_val) && return false

    nirreducible = 0
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
                nirreducible += !(alias_candidate isa Pair) || alias_var != alias_candidate[2]
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
        #preferred_var = pivot_var
        switch = false # we prefer `alias_var` by default, unless we switch
        diff_to_var = invview(var_to_diff)
        pivot_var′′::Union{Nothing, Int} = pivot_var′::Int = pivot_var
        alias_var′′::Union{Nothing, Int} = alias_var′::Int = alias_var
        # We prefer the higher differenitated variable. Note that `{⋅}′′` vars
        # could be `nothing` while `{⋅}′` vars are always `Int`.
        while (pivot_var′′ = diff_to_var[pivot_var′]) !== nothing
            pivot_var′ = pivot_var′′
            if (alias_var′′ = diff_to_var[alias_var′]) === nothing
                switch = true
                break
            end
            pivot_var′ = pivot_var′′
        end
        # If we have a tie, then we prefer the lower variable.
        if alias_var′′ === pivot_var′′ === nothing
            @assert pivot_var′ != alias_var′
            switch = pivot_var′ < alias_var′
        end
        if switch
            pivot_var, alias_var = alias_var, pivot_var
            pivot_val, alias_val = alias_val, pivot_val
        end

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
