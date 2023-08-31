# This code is derived from the Modia project and is licensed as follows:
# https://github.com/ModiaSim/Modia.jl/blob/b61daad643ef7edd0c1ccce6bf462c6acfb4ad1a/LICENSE

struct OrderedBitSet <: AbstractSet{Int}
    bitset::BitSet
    order::Vector{Int}
end
OrderedBitSet() = OrderedBitSet(BitSet(), Int[])
Base.iterate(o::OrderedBitSet) = Base.iterate(o.order)
Base.iterate(o::OrderedBitSet, state) = Base.iterate(o.order, state)
Base.in(a::Int, o::OrderedBitSet) = a in o.bitset
function Base.push!(o::OrderedBitSet, a::Int)
    if !(a in o.bitset)
        push!(o.bitset, a)
        push!(o.order, a)
    end
    o
end
function Base.empty!(o::OrderedBitSet)
    empty!(o.bitset)
    empty!(o.order)
    o
end
function Base.sort!(o::OrderedBitSet; kw...)
    sort!(o.order; kw...)
    o
end

function try_assign_eq!(ict::IncrementalCycleTracker, vj::Integer, eq::Integer)
    G = ict.graph
    add_edge_checked!(ict, Iterators.filter(!=(vj), 𝑠neighbors(G.graph, eq)), vj) do G
        G.matching[vj] = eq
        G.ne += length(𝑠neighbors(G.graph, eq)) - 1
    end
end

function try_assign_eq!(ict::IncrementalCycleTracker, vars, v_active, eq::Integer,
    condition::F = _ -> true) where {F}
    G = ict.graph
    for vj in vars
        (vj in v_active && G.matching[vj] === unassigned && condition(vj)) || continue
        try_assign_eq!(ict, vj, eq) && return true
    end
    return false
end

function tearEquations!(ict::IncrementalCycleTracker, Gsolvable, es::Vector{Int},
    v_active::OrderedBitSet, isder′::F) where {F}
    check_der = isder′ !== nothing
    if check_der
        has_der = Ref(false)
        isder = let has_der = has_der, isder′ = isder′
            v -> begin
                r = isder′(v)
                has_der[] |= r
                r
            end
        end
    end
    # Heuristic: As a first pass, try to assign any equations that only have one
    # solvable variable.
    for only_single_solvable in (true, false)
        for eq in es  # iterate only over equations that are not in eSolvedFixed
            vs = Gsolvable[eq]
            ((length(vs) == 1) ⊻ only_single_solvable) && continue
            if check_der
                # if there're differentiated variables, then only consider them
                try_assign_eq!(ict, vs, v_active, eq, isder)
                if has_der[]
                    has_der[] = false
                    continue
                end
            end
            try_assign_eq!(ict, vs, v_active, eq)
        end
    end

    return ict
end

function tear_graph_block_modia!(var_eq_matching, ict, solvable_graph, eqs, vars,
    isder::F, solved_eq) where {F}
    tearEquations!(ict, solvable_graph.fadjlist, eqs, vars, isder)
    for var in vars
        eq = var_eq_matching[var] = ict.graph.matching[var]
        if eq isa Int
            push!(solved_eq, eq)
        end
    end
    return nothing
end

function tear_graph_modia(structure::SystemStructure, isder::F = nothing,
    ::Type{U} = Unassigned;
    varfilter::F2 = v -> true,
    eqfilter::F3 = eq -> true,
    state_priority::F4 = nothing) where {F, U, F2, F3, F4}
    # It would be possible here to simply iterate over all variables and attempt to
    # use tearEquations! to produce a matching that greedily selects the minimal
    # number of torn variables. However, we can do this process faster if we first
    # compute the strongly connected components. In the absence of cycles and
    # non-solvability, a maximal matching on the original graph will give us an
    # optimal assignment. However, even with cycles, we can use the maximal matching
    # to give us a good starting point for a good matching and then proceed to
    # reverse edges in each scc to improve the solution. Note that it is possible
    # to have optimal solutions that cannot be found by this process. We will not
    # find them here [TODO: It would be good to have an explicit example of this.]

    @unpack graph, solvable_graph = structure
    var_eq_matching = maximal_matching(graph, eqfilter, varfilter, U)
    var_eq_matching = complete(var_eq_matching,
        max(length(var_eq_matching),
            maximum(x -> x isa Int ? x : 0, var_eq_matching)))
    full_var_eq_matching = copy(var_eq_matching)
    var_sccs::Vector{Union{Vector{Int}, Int}} = find_var_sccs(graph, var_eq_matching)
    vargraph = DiCMOBiGraph{true}(graph)
    ict = IncrementalCycleTracker(vargraph; dir = :in)

    ieqs = Int[]
    filtered_vars = OrderedBitSet()
    solved_eq = BitSet()
    if state_priority !== nothing
        sort!(var_sccs, by = Base.Fix1(sum, state_priority))
    end
    for vars in var_sccs
        for var in vars
            if varfilter(var)
                push!(filtered_vars, var)
                for eq in 𝑑neighbors(graph, var)
                    if !(eq in solved_eq)
                        push!(ieqs, eq)
                    end
                end
            end
            var_eq_matching[var] = unassigned
        end
        if state_priority !== nothing
            sort!(filtered_vars, by = state_priority)
        end
        tear_graph_block_modia!(var_eq_matching, ict, solvable_graph, ieqs,
            filtered_vars,
            isder, solved_eq)

        # clear cache
        vargraph.ne = 0
        for var in vars
            vargraph.matching[var] = unassigned
        end
        empty!(ieqs)
        empty!(filtered_vars)
    end
    return var_eq_matching, full_var_eq_matching
end
