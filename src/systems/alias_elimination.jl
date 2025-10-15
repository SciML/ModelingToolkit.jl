using SymbolicUtils: Rewriters
using Graphs.Experimental.Traversals
function alias_eliminate_graph!(state::TransformationState; kwargs...)
    mm = linear_subsys_adjmat!(state; kwargs...)
    if size(mm, 1) == 0
        return mm
    end
    @unpack graph, var_to_diff, solvable_graph = state.structure
    mm = alias_eliminate_graph!(state, mm; kwargs...)
    s = state.structure
    for (ei, e) in enumerate(mm.nzrows)
        set_neighbors!(s.graph, e, mm.row_cols[ei])
    end
    if s.solvable_graph isa BipartiteGraph{Int, Nothing}
        for (ei, e) in enumerate(mm.nzrows)
            set_neighbors!(s.solvable_graph, e, mm.row_cols[ei])
        end
    end
    return mm
end
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
    mm = alias_eliminate_graph!(state; kwargs...)
    fullvars = state.fullvars
    @unpack var_to_diff, graph, solvable_graph = state.structure
    subs = Dict{SymbolicT, SymbolicT}()
    to_expand = Int[]
    diff_to_var = invview(var_to_diff)
    dels = Int[]
    eqs = collect(equations(state))
    resize!(eqs, nsrcs(graph))
    __trivial_eq_rhs = let fullvars = fullvars
        function trivial_eq_rhs(pair)
            var, coeff = pair
            iszero(coeff) && return Symbolics.COMMON_ZERO
            return coeff * fullvars[var]
        end
    end
    for (ei, e) in enumerate(mm.nzrows)
        vs = 𝑠neighbors(graph, e)
        if isempty(vs)
            push!(dels, e)
        else
            rhs = mapfoldl(__trivial_eq_rhs, +, pairs(nonzerosmap(@view mm[ei, :])))
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
            cursor += 1
            old_to_new_eq[i] = -1
            continue
        end
        idx += 1
        old_to_new_eq[i] = idx
    end
    n_new_eqs = idx
    eqs_to_update = BitSet()
    for ieq in eqs_to_update
        eq = eqs[ieq]
        eqs[ieq] = substitute(eq, subs)
    end
    new_nparentrows = nsrcs(graph)
    new_row_cols = eltype(mm.row_cols)[]
    new_row_vals = eltype(mm.row_vals)[]
    new_nzrows = Int[]
    for (i, eq) in enumerate(mm.nzrows)
        old_to_new_eq[eq] > 0 || continue
        push!(new_row_cols, mm.row_cols[i])
        push!(new_row_vals, mm.row_vals[i])
        push!(new_nzrows, old_to_new_eq[eq])
    end
    mm = typeof(mm)(new_nparentrows, mm.ncols, new_nzrows, new_row_cols, new_row_vals)
    for old_ieq in to_expand
        ieq = old_to_new_eq[old_ieq]
        eqs[ieq] = expand_derivatives(eqs[ieq])
    end
    diff_to_var = invview(var_to_diff)
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
    new_var_to_diff = DiffGraph(length(var_to_diff))
    for v in 1:length(var_to_diff)
        new_var_to_diff[v] = var_to_diff[v]
    end
    state.structure.graph = new_graph
    state.structure.solvable_graph = new_solvable_graph
    state.structure.eq_to_diff = new_eq_to_diff
    state.structure.var_to_diff = new_var_to_diff
    sys = state.sys
    @set! sys.eqs = eqs
    state.sys = sys
    if mm isa SparseMatrixCLIL{BigInt, Int}
        return invalidate_cache!(sys), mm
    else
        return invalidate_cache!(sys), mm
    end
end
""""""
@inline function find_first_linear_variable(M::SparseMatrixCLIL,
        range,
        mask,
        constraint)
    eadj = M.row_cols
    @inbounds for i in range
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
    @inbounds for i in range
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
count_nonzeros(a::AbstractArray) = count(!iszero, a)
count_nonzeros(a::CLILVector) = nnz(a)
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
function aag_bareiss!(structure, mm_orig::SparseMatrixCLIL{T, Ti}) where {T, Ti}
    @unpack graph, var_to_diff = structure
    mm = copy(mm_orig)
    linear_equations_set = BitSet(mm_orig.nzrows)
    is_algebraic = let var_to_diff = var_to_diff
        v -> var_to_diff[v] === nothing === invview(var_to_diff)[v]
    end
    is_linear_variables = is_algebraic.(1:length(var_to_diff))
    is_highest_diff = computed_highest_diff_variables(structure)
    for i in 𝑠vertices(graph)
        (i in linear_equations_set && all(is_algebraic, 𝑠neighbors(graph, i))) &&
            continue
        for j in 𝑠neighbors(graph, i)
            is_linear_variables[j] = false
        end
    end
    solvable_variables = findall(is_linear_variables)
    local bar
    try
        bar = do_bareiss!(mm, mm_orig, is_linear_variables, is_highest_diff)
    catch e
        e isa OverflowError || rethrow(e)
        mm = convert(SparseMatrixCLIL{BigInt, Ti}, mm_orig)
        bar = do_bareiss!(mm, mm_orig, is_linear_variables, is_highest_diff)
    end
    if mm isa SparseMatrixCLIL{BigInt, Ti}
        return mm, solvable_variables, bar
    else
        return mm, solvable_variables, bar
    end
end
function do_bareiss!(M, Mold, is_linear_variables, is_highest_diff)
    rank1r = Ref{Union{Nothing, Int}}(nothing)
    rank2r = Ref{Union{Nothing, Int}}(nothing)
    find_pivot = let rank1r = rank1r
        (M, k) -> begin
            if rank1r[] === nothing
                r = find_masked_pivot(is_linear_variables, M, k)
                r !== nothing && return r
                rank1r[] = k - 1
            end
            if rank2r[] === nothing
                r = find_masked_pivot(is_highest_diff, M, k)
                r !== nothing && return r
                rank2r[] = k - 1
            end
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
    rank3, = bareiss!(M, bareiss_ops; find_pivot = find_and_record_pivot)
    rank2 = something(rank2r[], rank3)
    rank1 = something(rank1r[], rank2)
    (rank1, rank2, rank3, pivots)
end
function alias_eliminate_graph!(state::TransformationState, ils::SparseMatrixCLIL;
        fully_determined = true, kwargs...)
    @unpack structure = state
    @unpack graph, solvable_graph, var_to_diff, eq_to_diff = state.structure
    ils, solvable_variables, (rank1, rank2, rank3, pivots) = aag_bareiss!(structure, ils)
    if fully_determined == true
        rk1vars = BitSet(@view pivots[1:rank1])
        for v in solvable_variables
            v in rk1vars && continue
            @set! ils.nparentrows += 1
            push!(ils.nzrows, ils.nparentrows)
            push!(ils.row_cols, [v])
            push!(ils.row_vals, [convert(eltype(ils), 1)])
            add_vertex!(graph, SRC)
            add_vertex!(solvable_graph, SRC)
            add_edge!(graph, ils.nparentrows, v)
            add_edge!(solvable_graph, ils.nparentrows, v)
            add_vertex!(eq_to_diff)
        end
    end
    return ils
end
function exactdiv(a::Integer, b)
    d, r = divrem(a, b)
    @assert r == 0
    return d
end
swap!(v, i, j) = v[i], v[j] = v[j], v[i]
""""""
function topsort_equations(eqs::Vector{Equation}, unknowns::Vector{SymbolicT}; check = true)
    graph, assigns = observed2graph(eqs, unknowns)
    neqs = length(eqs)
    degrees = zeros(Int, neqs)
    for 𝑠eq in 1:length(eqs)
        var = assigns[𝑠eq]
        for 𝑑eq in 𝑑neighbors(graph, var)
            degrees[𝑑eq] += 1
        end
    end
    q = Queue{Int}(neqs)
    for (i, d) in enumerate(degrees)
        @static if pkgversion(DataStructures) >= v"0.19"
            d == 0 && push!(q, i)
        else
            d == 0 && enqueue!(q, i)
        end
    end
    idx = 0
    ordered_eqs = similar(eqs, 0)
    sizehint!(ordered_eqs, neqs)
    while !isempty(q)
        @static if pkgversion(DataStructures) >= v"0.19"
            𝑠eq = popfirst!(q)
        else
            𝑠eq = dequeue!(q)
        end
        idx += 1
        push!(ordered_eqs, eqs[𝑠eq])
        var = assigns[𝑠eq]
        for 𝑑eq in 𝑑neighbors(graph, var)
            degree = degrees[𝑑eq] = degrees[𝑑eq] - 1
            @static if pkgversion(DataStructures) >= v"0.19"
                degree == 0 && push!(q, 𝑑eq)
            else
                degree == 0 && enqueue!(q, 𝑑eq)
            end
        end
    end
    (check && idx != neqs) && throw(ArgumentError("The equations have at least one cycle."))
    return ordered_eqs
end
function observed2graph(eqs::Vector{Equation}, unknowns::Vector{SymbolicT})::Tuple{BipartiteGraph{Int, Nothing}, Vector{Int}}
    graph = BipartiteGraph(length(eqs), length(unknowns))
    v2j = Dict{SymbolicT, Int}(unknowns .=> 1:length(unknowns))
    assigns = similar(eqs, Int)
    vars = Set{SymbolicT}()
    for (i, eq) in enumerate(eqs)
        lhs_j = get(v2j, eq.lhs, nothing)
        lhs_j === nothing &&
            throw(ArgumentError("The lhs $(eq.lhs) of $eq, doesn't appear in unknowns."))
        assigns[i] = lhs_j
        empty!(vars)
        SU.search_variables!(vars, eq.rhs; is_atomic = OperatorIsAtomic{SU.Operator}())
        for v in vars
            j = get(v2j, v, nothing)
            if j isa Int
                add_edge!(graph, i, j)
            end
        end
    end
    return graph, assigns
end
