using SymbolicUtils: Rewriters

const KEEP = typemin(Int)

include("compat/bareiss.jl")

function alias_elimination(sys)
    sys = initialize_system_structure(sys; quick_cancel=true)
    s = structure(sys)

    mm = linear_subsys_adjmat(sys)
    size(mm, 1) == 0 && return sys # No linear subsystems

    ag, mm = alias_eliminate_graph!(s.graph, s.varassoc, mm)

    @unpack fullvars, graph = s

    subs = OrderedDict()
    for (v, (coeff, alias)) in pairs(ag)
        subs[fullvars[v]] = iszero(coeff) ? 0 : coeff * fullvars[alias]
    end

    dels = Set{Int}()
    eqs = copy(equations(sys))
    for (ei, e) in enumerate(mm.nzrows)
        vs = 𝑠neighbors(graph, e)
        if isempty(vs)
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

    dict = Dict(subs)
    for (ieq, eq) in enumerate(eqs)
        eqs[ieq] = eq.lhs ~ fixpoint_sub(eq.rhs, dict)
    end

    newstates = []
    sts = states(sys)
    for j in eachindex(fullvars)
        if !(j in keys(ag))
            isdervar(s, j) || push!(newstates, fullvars[j])
        end
    end

    @set! sys.eqs = eqs
    @set! sys.states = newstates
    @set! sys.observed = [observed(sys); [lhs ~ rhs for (lhs, rhs) in pairs(subs)]]
    @set! sys.structure = nothing
    return sys
end

"""
$(SIGNATURES)

Find the first linear variable such that `𝑠neighbors(adj, i)[j]` is true given
the `constraint`.
"""
@inline function find_first_linear_variable(
        M::SparseMatrixCLIL,
        range,
        mask,
        constraint,
    )
    eadj = M.row_cols
    for i in range
        vertices = eadj[i]
        if constraint(length(vertices))
            for (j, v) in enumerate(vertices)
                (mask === nothing || mask[v]) && return (CartesianIndex(i, v), M.row_vals[i][j])
            end
        end
    end
    return nothing
end

@inline function find_first_linear_variable(
        M::AbstractMatrix,
        range,
        mask,
        constraint,
    )
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
    r = find_first_linear_variable(M, k:size(M,1), variables, isequal(1))
    r !== nothing && return r
    r = find_first_linear_variable(M, k:size(M,1), variables, isequal(2))
    r !== nothing && return r
    r = find_first_linear_variable(M, k:size(M,1), variables, _->true)
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

function Base.getindex(ag::AliasGraph, i::Integer)
    r = ag.aliasto[i]
    r === nothing && throw(KeyError(i))
    coeff, var = (sign(r), abs(r))
    if var in keys(ag)
        # Amortized lookup. Check if since we last looked this up, our alias was
        # itself aliased. If so, adjust adjust the alias table.
        ac, av = ag[var]
        nc = ac * coeff
        ag.aliasto[var] = nc > 0 ? av : -av
    end
    return (coeff, var)
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
    return 0=>0
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
Base.in(i::Int, agk::AliasGraphKeySet) = agk.ag.aliasto[i] !== nothing

count_nonzeros(a::AbstractArray) = count(!iszero, a)

# N.B.: Ordinarily sparse vectors allow zero stored elements.
# Here we have a guarantee that they won't, so we can make this identification
count_nonzeros(a::SparseVector) = nnz(a)

function alias_eliminate_graph!(graph, varassoc, mm_orig::SparseMatrixCLIL)
    invvarassoc = inverse_mapping(varassoc)

    mm = copy(mm_orig)
    is_linear_equations = falses(size(AsSubMatrix(mm_orig), 1))
    for e in mm_orig.nzrows
        is_linear_equations[e] = true
    end

    is_not_potential_state = iszero.(varassoc)
    is_linear_variables = copy(is_not_potential_state)
    for i in 𝑠vertices(graph); is_linear_equations[i] && continue
        for j in 𝑠neighbors(graph, i)
            is_linear_variables[j] = false
        end
    end
    solvable_variables = findall(is_linear_variables)

    function do_bareiss!(M, Mold=nothing)
        rank1 = rank2 = nothing
        pivots = Int[]
        function find_pivot(M, k)
            if rank1 === nothing
                r = find_masked_pivot(is_linear_variables, M, k)
                r !== nothing && return r
                rank1 = k - 1
            end
            if rank2 === nothing
                r = find_masked_pivot(is_not_potential_state, M, k)
                r !== nothing && return r
                rank2 = k - 1
            end
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
        bareiss_ops = ((M,i,j)->nothing, myswaprows!, bareiss_update_virtual_colswap_mtk!, bareiss_zero!)
        rank3 = bareiss!(M, bareiss_ops; find_pivot=find_and_record_pivot)
        rank1 = something(rank1, rank3)
        rank2 = something(rank2, rank3)
        (rank1, rank2, rank3, pivots)
    end

    # mm2 = Array(copy(mm))
    # @show do_bareiss!(mm2)
    # display(mm2)

    # Step 1: Perform bareiss factorization on the adjacency matrix of the linear
    #         subsystem of the system we're interested in.
    (rank1, rank2, rank3, pivots) = do_bareiss!(mm, mm_orig)

    # Step 2: Simplify the system using the bareiss factorization
    ag = AliasGraph(size(mm, 2))
    for v in setdiff(solvable_variables, @view pivots[1:rank1])
        ag[v] = 0
    end

    # kind of like the backward substitution
    lss!(ei::Integer) = locally_structure_simplify!((@view mm[ei, :]), pivots[ei], ag, invvarassoc[pivots[ei]] == 0)

    # Step 2.1: Go backwards, collecting eliminated variables and substituting
    #         alias as we go.
    foreach(lss!, reverse(1:rank2))

    # Step 2.2: Sometimes bareiss can make the equations more complicated.
    #         Go back and check the original matrix. If this happened,
    #         Replace the equation by the one from the original system,
    #         but be sure to also run lss! again, since we only ran that
    #         on the bareiss'd matrix, not the original one.
    reduced = mapreduce(|, 1:rank2; init=false) do ei
        if count_nonzeros(@view mm_orig[ei, :]) < count_nonzeros(@view mm[ei, :])
            mm[ei, :] = @view mm_orig[ei, :]
            return lss!(ei)
        end
        return false
    end

    # Step 2.3: Iterate to convergance.
    #         N.B.: `lss!` modifies the array.
    # TODO: We know exactly what variable we eliminated. Starting over at the
    #       start is wasteful. We can lookup which equations have this variable
    #       using the graph.
    reduced && while any(lss!, 1:rank2); end

    # Step 3: Reflect our update decitions back into the graph
    for (ei, e) in enumerate(mm.nzrows)
        graph.fadjlist[e] = mm.row_cols[ei]
    end

    return ag, mm
end

iszeroterm(v_types, v) = v_types[v] == 0
isirreducible(v_types, v) = v_types[v] == KEEP
isalias(v_types, v) = v_types[v] > 0 && !isirreducible(v_types, v)
alias(v_types, v) = v_types[v]
negalias(v_types, v) = -v_types[v]

function exactdiv(a::Integer, b::Integer)
    d, r = divrem(a, b)
    @assert r == 0
    return d
end

function locally_structure_simplify!(adj_row, pivot_col, ag, may_eliminate)
    pivot_val = adj_row[pivot_col]
    iszero(pivot_val) && return false

    nirreducible = 0
    alias_candidate = 0

    # N.B.: Assumes that the non-zeros iterator is robust to modification
    # of the underlying array datastructure.
    for (var, val) in pairs(nonzerosmap(adj_row))
        # Go through every variable/coefficient in this row and apply all aliases
        # that we have so far accumulated in `ag`, updating the adj_row as
        # we go along.
        var == pivot_col && continue
        iszero(val) && continue
        alias = get(ag, var, nothing)
        if alias === nothing
            nirreducible += 1
            alias_candidate = val => var
            continue
        end
        (coeff, alias_var) = alias
        adj_row[var] = 0
        if alias_var != 0
            val *= coeff
            new_coeff = (adj_row[alias_var] += val)
            if alias_var < var
                # If this adds to a coeff that was not previously accounted for, and
                # we've already passed it, make sure to count it here. We're
                # relying on `var` being produced in sorted order here.
                nirreducible += 1
                alias_candidate = new_coeff => alias_var
            end
        end
    end

    if may_eliminate && nirreducible <= 1
        # There were only one or two terms left in the equation (including the pivot variable).
        # We can eliminate the pivot variable.
        if alias_candidate !== 0
            alias_candidate = -exactdiv(alias_candidate[1], pivot_val) => alias_candidate[2]
        end
        ag[pivot_col] = alias_candidate
        zero!(adj_row)
        return true
    end
    return false
end

swap!(v, i, j) = v[i], v[j] = v[j], v[i]

function getcoeff(vars, coeffs, var)
    for (vj, v) in enumerate(vars)
        v == var && return coeffs[vj]
    end
    return 0
end

function inverse_mapping(assoc)
    invassoc = zeros(Int, length(assoc))
    for (i, v) in enumerate(assoc)
        v <= 0 && continue
        invassoc[v] = i
    end
    return invassoc
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
function topsort_equations(eqs, states; check=true)
    graph, assigns = observed2graph(eqs, states)
    neqs = length(eqs)
    degrees = zeros(Int, neqs)

    for 𝑠eq in 1:length(eqs); var = assigns[𝑠eq]
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
    ordered_eqs = similar(eqs, 0); sizehint!(ordered_eqs, neqs)
    while !isempty(q)
        𝑠eq = dequeue!(q)
        idx+=1
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
        lhs_j === nothing && throw(ArgumentError("The lhs $(eq.lhs) of $eq, doesn't appear in states."))
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
    map(eq->eq.lhs ~ sub(eq.rhs), eqs)
end
