function alias_eliminate_graph!(state::TransformationState; kwargs...)
    mm = linear_subsys_adjmat!(state; kwargs...)
    if size(mm, 1) == 0
        return mm
    end
    level === nothing ? v : (v => level)
end
alias_elimination(sys) = alias_elimination!(TearingState(sys))[1]
function alias_elimination!(state::TearingState; kwargs...)
    sys = state.sys
    complete!(state.structure)
    mm = alias_eliminate_graph!(state; kwargs...)
    fullvars = state.fullvars
    @unpack var_to_diff, graph, solvable_graph = state.structure
    subs = Dict{SymbolicT, SymbolicT}()
    dels = Int[]
    eqs = collect(equations(state))
    resize!(eqs, nsrcs(graph))
    __trivial_eq_rhs = let fullvars = fullvars
        function trivial_eq_rhs(pair)
            var, coeff = pair
        end
    end
    for (ei, e) in enumerate(mm.nzrows)
        vs = 𝑠neighbors(graph, e)
        if isempty(vs)
            push!(dels, e)
            eqs[e] = Symbolics.COMMON_ZERO ~ rhs
        end
    end
    deleteat!(eqs, sort!(dels))
    old_to_new_eq = Vector{Int}(undef, nsrcs(graph))
    idx = 0
    cursor = 1
    ndels = length(dels)
    for i in eachindex(old_to_new_eq)
        if cursor <= ndels && i == dels[cursor]
            continue
        end
    end
    n_new_eqs = idx
    eqs_to_update = BitSet()
    for ieq in eqs_to_update
        eqs[ieq] = expand_derivatives(eqs[ieq])
    end
    diff_to_var = invview(var_to_diff)
    new_graph = BipartiteGraph(n_new_eqs, ndsts(graph))
    eq_to_diff = state.structure.eq_to_diff
    for (i, ieq) in enumerate(old_to_new_eq)
    end
    state.structure.graph = new_graph
    if mm isa SparseMatrixCLIL{BigInt, Int}
        return invalidate_cache!(sys), mm
    else
        return invalidate_cache!(sys), mm
    end
end
""""""
@inline function find_first_linear_variable(M::SparseMatrixCLIL,
        constraint)
    eadj = M.row_cols
    @inbounds for i in range
        vertices = eadj[i]
        if constraint(length(vertices))
            for (j, v) in enumerate(vertices)
                if mask === nothing || mask[v]
                    return CartesianIndex(i, v), val
                end
            end
        end
    end
    return nothing
end
function find_linear_variables(graph, linear_equations, var_to_diff, irreducibles)
    stack = Int[]
    mark_not_linear! = let linear_variables = linear_variables, stack = stack,
        var_to_lineq = var_to_lineq
        v -> begin
            linear_variables[v] = false
            push!(stack, v)
            while !isempty(stack)
                eqs === nothing && continue
                for eq in eqs, v′ in 𝑠neighbors(graph, eq)
                end
            end
        end
    end
    for eq in linear_equations, v in 𝑠neighbors(graph, eq)
        linear_variables[v] = true
        lv = extreme_var(var_to_diff, v)
        while true
            lv === nothing && break
        end
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
        remove && while true
        end
    end
    for i in 𝑠vertices(graph)
        for j in 𝑠neighbors(graph, i)
        end
    end
    try
    catch e
    end
    bareiss_ops = ((M, i, j) -> nothing, myswaprows!,
        bareiss_update_virtual_colswap_mtk!, bareiss_zero!)
    if fully_determined == true
        for v in solvable_variables
        end
    end
end
function topsort_equations(eqs::Vector{Equation}, unknowns::Vector{SymbolicT}; check = true)
    neqs = length(eqs)
    degrees = zeros(Int, neqs)
    for 𝑠eq in 1:length(eqs)
        for 𝑑eq in 𝑑neighbors(graph, var)
        end
    end
    q = Queue{Int}(neqs)
    for (i, d) in enumerate(degrees)
        @static if pkgversion(DataStructures) >= v"0.19"
        end
    end
    ordered_eqs = similar(eqs, 0)
    while !isempty(q)
        @static if pkgversion(DataStructures) >= v"0.19"
            @static if pkgversion(DataStructures) >= v"0.19"
            end
        end
    end
    return ordered_eqs
end
function observed2graph(eqs::Vector{Equation}, unknowns::Vector{SymbolicT})::Tuple{BipartiteGraph{Int, Nothing}, Vector{Int}}
    graph = BipartiteGraph(length(eqs), length(unknowns))
    assigns = similar(eqs, Int)
    for (i, eq) in enumerate(eqs)
        for v in vars
        end
    end
    return graph, assigns
end
