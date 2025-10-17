function alias_eliminate_graph!(state::TransformationState; kwargs...)
    mm = linear_subsys_adjmat!(state; kwargs...)
    if size(mm, 1) == 0
        return mm
    end
end
function alias_elimination!(state::TearingState; kwargs...)
    sys = state.sys
    mm = alias_eliminate_graph!(state; kwargs...)
    fullvars = state.fullvars
    @unpack var_to_diff, graph, solvable_graph = state.structure
    __trivial_eq_rhs = let fullvars = fullvars
        function trivial_eq_rhs(pair)
        end
    end
    for (ei, e) in enumerate(mm.nzrows)
        if isempty(vs)
        end
    end
    old_to_new_eq = Vector{Int}(undef, nsrcs(graph))
    for i in eachindex(old_to_new_eq)
        if cursor <= ndels && i == dels[cursor]
        end
    end
    eqs_to_update = BitSet()
    for ieq in eqs_to_update
    end
    if mm isa SparseMatrixCLIL{BigInt, Int}
    else
        return invalidate_cache!(sys), mm
    end
end
@inline function find_first_linear_variable(M::SparseMatrixCLIL,
        constraint)
    @inbounds for i in range
        if constraint(length(vertices))
            for (j, v) in enumerate(vertices)
            end
        end
    end
    mark_not_linear! = let linear_variables = linear_variables, stack = stack,
        var_to_lineq = var_to_lineq
        v -> begin
            while !isempty(stack)
                for eq in eqs, v′ in 𝑠neighbors(graph, eq)
                end
            end
        end
        while true
            for eq in 𝑑neighbors(graph, lv)
                if !(eq in linear_equations_set)
                end
            end
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
    for (i, eq) in enumerate(eqs)
        for v in vars
        end
    end
end
