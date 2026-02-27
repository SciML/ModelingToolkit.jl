struct CacheWriter{F}
    fn::F
end

function (cw::CacheWriter)(p, sols)
    return cw.fn(p.caches, sols, p)
end

const SCCCacheVarsExprsElT = Dict{TypeT, Vector{SymbolicT}}

function CacheWriter(
        sys::AbstractSystem, buffer_types::Vector{TypeT},
        exprs::SCCCacheVarsExprsElT, solsyms, obseqs::Vector{Equation};
        eval_expression = false, eval_module = @__MODULE__, cse = true, sparse = false
    )
    ps = parameters(sys; initial_parameters = true)
    rps = reorder_parameters(sys, ps)
    obs_assigns = [eq.lhs ← eq.rhs for eq in obseqs]
    body = map(eachindex(buffer_types), buffer_types) do i, T
        Symbol(:tmp, i) ← SetArray(true, :(out[$i]), get(exprs, T, []))
    end

    function argument_name(i::Int)
        if i <= length(solsyms)
            return :($(generated_argument_name(1))[$i])
        end
        return generated_argument_name(i - length(solsyms))
    end
    array_assignments = array_variable_assignments(solsyms...; argument_name)
    fn = build_function_wrapper(
        sys, nothing, :out,
        DestructuredArgs(DestructuredArgs.(solsyms), generated_argument_name(1)),
        rps...; p_start = 3, p_end = length(rps) + 2,
        expression = Val{true}, add_observed = false, cse,
        extra_assignments = [array_assignments; obs_assigns; body]
    )
    fn = eval_or_rgf(fn; eval_expression, eval_module)
    fn = GeneratedFunctionWrapper{(3, 3, is_split(sys))}(fn, nothing)
    return CacheWriter(fn)
end

"""
    $TYPEDSIGNATURES

Subset a system to have only the given unknowns `vscc` and equations `escc`. Observed
equations are subset, delete accordingly. Requires that `sys` is complete and flattened.

# Keyword arguments

- `available_vars`: A list of variables that the subset system should assume are precomputed
  or already available. Will be mutated with the unknowns and observables of the subset. This
  is useful for SCC decomposition.
- `prevobsidxs`: Indices of observed equations that the subset system should assume are
  precomputed or already available. Will be appended with the indices of observed equations
  required by this subset.
"""
function subset_system(
        sys::System, vscc::Vector{Int}, escc::Vector{Int};
        available_vars = Set{SymbolicT}(), prevobsidxs = Int[]
    )
    check_complete(sys, "subset_system")
    @assert isempty(get_systems(sys)) "`subset_system` requires a flattened system"

    dvs = unknowns(sys)
    ps = parameters(sys)
    eqs = equations(sys)
    obs = observed(sys)

    # subset unknowns and equations
    _dvs = dvs[vscc]
    _eqs = eqs[escc]
    # get observed equations required by this SCC
    union!(available_vars, _dvs)
    obsidxs = observed_equations_used_by(sys, _eqs; available_vars)
    # the ones used by previous SCCs can be precomputed into the cache
    setdiff!(obsidxs, prevobsidxs)
    _obs = obs[obsidxs]
    _observables = SymbolicT[]
    for eq in _obs
        push!(_observables, eq.lhs)
    end
    union!(available_vars, _observables)
    append!(prevobsidxs, obsidxs)

    subsys = ConstructionBase.setproperties(
        sys; unknowns = _dvs, eqs = _eqs, observed = _obs,
        parameter_bindings_graph = get_parameter_bindings_graph(sys), complete = true
    )
    if get_index_cache(sys) !== nothing
        @set! subsys.index_cache = subset_unknowns_observed(
            get_index_cache(sys), sys, _dvs, getproperty.(_obs, (:lhs,))
        )
    end
    cached_param_arr_assigns = check_mutable_cache(
        sys, MTKBase.ParameterArrayAssignments, MTKBase.ParameterArrayAssignments, nothing
    )
    if cached_param_arr_assigns isa MTKBase.ParameterArrayAssignments
        store_to_mutable_cache!(
            subsys, MTKBase.ParameterArrayAssignments, cached_param_arr_assigns
        )
    end

    return subsys
end

const BlockIdxsT = typeof(BlockVector{Int}(undef_blocks, Int[]))

mutable struct SCCDecomposition
    const subsystems::Vector{System}
    const var_sccs::Vector{Vector{Int}}
    const eq_sccs::Vector{Vector{Int}}
    const islinear::BitVector
    const hints::Vector{StructuralHint.Type}
    # Cache buffer types and corresponding sizes. Stored as a pair of arrays instead of a
    # dict to maintain a consistent order of buffers across SCCs
    const cachetypes::Vector{TypeT}
    const cachesizes::Vector{Int}
    # explicitfun! related information for each SCC
    # We need to compute buffer sizes before doing any codegen
    const scc_cachevars::Vector{SCCCacheVarsExprsElT}
    const scc_cacheexprs::Vector{SCCCacheVarsExprsElT}
    obsidxs::BlockIdxsT
    const obsidxs_for_cacheexprs::Vector{Vector{Int}}
end

function SCCDecomposition()
    return SCCDecomposition(
        System[], Vector{Int}[], Vector{Int}[], BitVector(), StructuralHint.Type[],
        TypeT[], Int[], SCCCacheVarsExprsElT[], SCCCacheVarsExprsElT[],
        BlockVector{Int}(undef_blocks, Int[]), Vector{Int}[]
    )
end

function SCCDecomposition(sys::System, var_sccs::Vector{Vector{Int}}, eq_sccs::Vector{Vector{Int}})
    active_decomposition = SCCDecomposition()
    final_decomposition = SCCDecomposition()
    available_vars = Set{SymbolicT}()
    ts = get_tearing_state(sys)::TearingState
    icg = BipartiteGraphs.InducedCondensationGraph(ts.structure.graph, var_sccs)
    active = Set{Int}()
    for (i, (escc, vscc)) in enumerate(zip(eq_sccs, var_sccs))
        nbors = Set{Int}(collect(Graphs.inneighbors(icg, i)))
        intersect!(nbors, active)
        while !isempty(nbors)
            finalize_scc!(final_decomposition, active_decomposition, first(nbors), active, nbors)
        end
        blockpush!(active_decomposition.obsidxs, Int[])
        subsys = subset_system(sys, vscc, escc; available_vars, prevobsidxs = active_decomposition.obsidxs)
        push!(active_decomposition.subsystems, subsys)
        push!(active_decomposition.var_sccs, vscc)
        push!(active_decomposition.eq_sccs, escc)
        push!(active_decomposition.hints, StructuralHint.NoHint())
        push!(active_decomposition.islinear, calculate_A_b(subsys; throw = false) !== nothing)
        push!(active, i)
    end

    while !isempty(active)
        finalize_scc!(final_decomposition, active_decomposition, first(active), active, active)
    end

    return final_decomposition
end

function copy_scc!(dst::SCCDecomposition, src::SCCDecomposition, tgt::Int)
    blockpush!(dst.obsidxs, copy(view(src.obsidxs, Block(tgt))))
    push!(dst.subsystems, src.subsystems[tgt])
    push!(dst.var_sccs, src.var_sccs[tgt])
    push!(dst.eq_sccs, src.eq_sccs[tgt])
    push!(dst.hints, src.hints[tgt])
    push!(dst.islinear, src.islinear[tgt])
    return dst
end

function finalize_scc!(final_decomposition::SCCDecomposition, active_decomposition::SCCDecomposition, i::Int, active::Set{Int}, nbors::Set{Int}; linear_scc_combine_range::Int = 2)
    if !active_decomposition.islinear[i]
        copy_scc!(final_decomposition, active_decomposition, i)
        delete!(active, i)
        delete!(nbors, i)
        return
    end

    subsys = active_decomposition.subsystems[i]
    neqs = length(equations(subsys))
    if neqs == 1
        to_merge = Int[]
        for j in active
            if active_decomposition.islinear[j] && length(equations(active_decomposition.subsystems[j])) == 1
                push!(to_merge, j)
            end
        end
        if length(to_merge) > 1
            active_decomposition.hints[to_merge[1]] = StructuralHint.Diagonal()
        end
        # Don't remove old to avoid affecting ordering
        _collapse_into!(active_decomposition, to_merge[1], @view(to_merge[2:end]))
        copy_scc!(final_decomposition, active_decomposition, to_merge[1])
        setdiff!(active, to_merge)
        setdiff!(nbors, to_merge)
        # TODO: Merge this SCC again if it is small
        return
    end

    merge_candidates = Int[]
    for j in active
        active_decomposition.islinear[j] || continue
        neqs_j = length(equations(active_decomposition.subsystems[j]))
        neqs_j == 1 && continue
        push!(merge_candidates, j)
    end

    comparator = length ∘ equations ∘ Base.Fix1(getindex, active_decomposition.subsystems)
    sort!(merge_candidates; by = comparator)
    # mapping from size of SCC to number of times it occurs in the range `low..high`.
    # Effectively a sorted multiset.
    sizes_in_range = DataStructures.SortedDict{Int, Int}()
    while i in active
        empty!(sizes_in_range)
        low = 1
        high = 1
        best_low = 1
        best_high = 1
        while checkbounds(Bool, merge_candidates, high)
            neqs_high = comparator(high)
            sizes_in_range[neqs_high] = get(sizes_in_range, neqs_high, 0) + 1
            absdiff = first(last(sizes_in_range)) - first(first(sizes_in_range))
            while absdiff > linear_scc_combine_range
                neqs_low = comparator(low)
                low += 1
                n_low = sizes_in_range[neqs_low] -= 1
                if iszero(n_low)
                    delete!(sizes_in_range, neqs_low)
                end
                absdiff = first(last(sizes_in_range)) - first(first(sizes_in_range))
            end
            if (high - low) > (best_high - best_low)
                best_high = high
                best_low = low
            end
            high += 1
        end

        if best_low == best_high
            copy_scc!(final_decomposition, active_decomposition, merge_candidates[best_low])
            delete!(active, merge_candidates[best_low])
            delete!(nbors, merge_candidates[best_low])
            deleteat!(merge_candidates, best_low)
            continue
        end

        merge_target = merge_candidates[best_low]
        to_merge = view(merge_candidates, (best_low + 1):best_high)
        largest_collapsed_scc = max(maximum(comparator, to_merge), comparator(merge_target))
        band_size = largest_collapsed_scc - 1
        active_decomposition.hints[merge_target] = StructuralHint.Banded(band_size, band_size)
        _collapse_into!(active_decomposition, merge_target, to_merge)
        copy_scc!(final_decomposition, active_decomposition, merge_target)
        setdiff!(active, to_merge)
        setdiff!(nbors, to_merge)
        delete!(active, merge_target)
        delete!(nbors, merge_target)
        deleteat!(merge_candidates, best_low:best_high)
    end

    return nothing
end

function build_caches!(sys::System, decomposition::SCCDecomposition)
    banned_vars = Set{SymbolicT}()
    state = Dict()
    prev_obseqs = Equation[]
    # We use the unhacked system to avoid spurious array hack observed equations
    # from making it in when the elements of the array are spread across SCCs.
    unhacked_sys = unhack_system(sys)
    for i in eachindex(decomposition.subsystems)
        empty!(banned_vars)
        empty!(state)
        
        subsys = decomposition.subsystems[i]
        union!(banned_vars, unknowns(subsys))
        union!(banned_vars, observables(subsys))

        _obs = get_observed(subsys)
        for i in eachindex(_obs)
            _obs[i] = _obs[i].lhs ~ subexpressions_not_involving_vars!(
                _obs[i].rhs, banned_vars, state
            )
        end
        _eqs = get_eqs(subsys)
        for i in eachindex(_eqs)
            _eqs[i] = _eqs[i].lhs ~ subexpressions_not_involving_vars!(
                _eqs[i].rhs, banned_vars, state
            )
        end

        if decomposition.islinear[i]
            store_to_mutable_cache!(subsys, CachedLinearAb, nothing)
            # cached_ab = check_mutable_cache(
            #     subsys, CachedLinearAb, CachedLinearAb, nothing
            # )
            # if cached_ab isa CachedLinearAb
            #     subber = SU.Substituter{false}(state)
            #     I, J, V = findnz(cached_ab.A)
            #     map!(subber, V, V)
            #     map!(subber, cached_ab.b, cached_ab.b)
            # end
        end

        # map from symtype to cached variables and their expressions
        cachevars = SCCCacheVarsExprsElT()
        cacheexprs = SCCCacheVarsExprsElT()
        push!(decomposition.scc_cachevars, cachevars)
        push!(decomposition.scc_cacheexprs, cacheexprs)
        # observed of previous SCCs are in the cache
        # NOTE: When we get proper CSE, we can substitute these
        # and then use `subexpressions_not_involving_vars!`
        for eq in prev_obseqs
            T = symtype(eq.lhs)
            buf = get!(() -> SymbolicT[], cachevars, T)
            push!(buf, eq.lhs)

            buf = get!(() -> SymbolicT[], cacheexprs, T)
            push!(buf, eq.lhs)
        end

        append!(prev_obseqs, observed(subsys))

        for (k, v) in state
            k = unwrap(k)
            v = unwrap(v)
            T = symtype(k)
            buf = get!(() -> SymbolicT[], cachevars, T)
            push!(buf, v)
            buf = get!(() -> SymbolicT[], cacheexprs, T)
            push!(buf, k)
        end
        all_cacheexprs = reduce(vcat, values(cacheexprs); init = SymbolicT[])
        obsidxs_for_scc_cacheexprs = observed_equations_used_by(unhacked_sys, all_cacheexprs)
        push!(decomposition.obsidxs_for_cacheexprs, obsidxs_for_scc_cacheexprs)
        # update the sizes of cache buffers
        for (T, buf) in cachevars
            idx = findfirst(isequal(T), decomposition.cachetypes)
            if idx === nothing
                push!(decomposition.cachetypes, T)
                push!(decomposition.cachesizes, 0)
                idx = lastindex(decomposition.cachetypes)
            end
            decomposition.cachesizes[idx] = max(decomposition.cachesizes[idx], length(buf))
        end
    end
end

"""
    $TYPEDSIGNATURES

Make SCC `i` a combination of SCC `i` and SCCs in `js`, where `js` is an iterable of
integers. Only modifies SCC `i`. Does not change the number of SCCs stored.
"""
function _collapse_into!(decomposition::SCCDecomposition, i::Int, js)
    parent = decomposition.subsystems[i]
    new_eqs = copy(equations(parent))
    new_dvs = copy(unknowns(parent))
    new_obs = copy(observed(parent))
    cached_ab::Union{CachedLinearAb, Nothing} = if decomposition.islinear[i]
        calculate_A_b(parent)
        check_mutable_cache(parent, CachedLinearAb, CachedLinearAb, nothing)::CachedLinearAb
    else
        nothing
    end
    for j in js
        cur = decomposition.subsystems[j]
        append!(new_eqs, equations(cur))
        append!(new_dvs, unknowns(cur))
        append!(new_obs, observed(cur))
        decomposition.islinear[i] &= decomposition.islinear[j]
        if cached_ab isa CachedLinearAb && decomposition.islinear[j]
            A = cached_ab.A
            b = cached_ab.b
            calculate_A_b(cur)
            jcache = check_mutable_cache(cur, CachedLinearAb, CachedLinearAb, nothing)::CachedLinearAb
            A = blockdiag(A, jcache.A)
            b = vcat(b, jcache.b)
            cached_ab = CachedLinearAb(A, b)
        end

        append!(decomposition.var_sccs[i], decomposition.var_sccs[j])
        append!(decomposition.eq_sccs[i], decomposition.eq_sccs[j])
    end
    unique!(new_obs)
    new_parent = decomposition.subsystems[i] = ConstructionBase.setproperties(
        parent; eqs = new_eqs, unknowns = new_dvs, observed = new_obs
    )
    if cached_ab isa CachedLinearAb
        store_to_mutable_cache!(new_parent, CachedLinearAb, cached_ab)
    end

    parent_obsidxs = view(decomposition.obsidxs, Block(i))
    for j in js
        append!(parent_obsidxs, view(decomposition.obsidxs, Block(j)))
    end
    sort!(parent_obsidxs)
    unique!(parent_obsidxs)

    return nothing
end

struct SCCNonlinearFunction{iip} end

function SCCNonlinearFunction{iip}(
        decomposition::SCCDecomposition, i::Int, cachesyms, op; eval_expression = false,
        eval_module = @__MODULE__, cse = true, kwargs...
    ) where {iip}
    subsys = decomposition.subsystems[i]
    islin = decomposition.islinear[i]
    # generate linear problem instead
    if islin
        return LinearFunction{iip}(
            subsys; eval_expression, eval_module, cse, cachesyms,
            structural_hint = decomposition.hints[i], kwargs...
        )
    end
    rps = reorder_parameters(subsys)
    f = generate_rhs(
        subsys; expression = Val{false}, wrap_gfw = Val{true}, cachesyms,
        obsidxs_to_use = eachindex(observed(subsys))
    )

    return NonlinearFunction{iip}(f; sys = subsys)
end

function SciMLBase.SCCNonlinearProblem(sys::System, args...; kwargs...)
    return SCCNonlinearProblem{true}(sys, args...; kwargs...)
end

function SciMLBase.SCCNonlinearProblem{iip}(
        sys::System, op; eval_expression = false,
        eval_module = @__MODULE__, cse = true, u0_constructor = identity,
        missing_guess_value = default_missing_guess_value(), combine_sccs = true, kwargs...
    ) where {iip}
    if !iscomplete(sys) || get_tearing_state(sys) === nothing
        error("A simplified `System` is required. Call `mtkcompile` on the system before creating an `SCCNonlinearProblem`.")
    end

    if !is_split(sys)
        error("The system has been simplified with `split = false`. `SCCNonlinearProblem` is not compatible with this system. Pass `split = true` to `mtkcompile` to use `SCCNonlinearProblem`.")
    end

    ts = get_tearing_state(sys)
    sched = get_schedule(sys)
    if sched === nothing
        @warn "System is simplified but does not have a schedule. This should not happen."
        var_eq_matching, var_sccs = StructuralTransformations.algebraic_variables_scc(ts)
        condensed_graph = MatchedCondensationGraph(
            DiCMOBiGraph{true}(
                complete(ts.structure.graph),
                complete(var_eq_matching)
            ),
            var_sccs
        )
        toporder = topological_sort_by_dfs(condensed_graph)
        var_sccs = var_sccs[toporder]
        eq_sccs = map(Base.Fix1(getindex, var_eq_matching), var_sccs)
    else
        var_sccs = sched.var_sccs
        # Equations are already in the order of SCCs
        eq_sccs = length.(var_sccs)
        cumsum!(eq_sccs, eq_sccs)
        eq_sccs = map(enumerate(eq_sccs)) do (i, lasti)
            i == 1 ? collect(1:lasti) : collect((eq_sccs[i - 1] + 1):lasti)
        end
    end

    if length(var_sccs) == 1
        if calculate_A_b(sys; throw = false) !== nothing
            linprob = LinearProblem{iip}(
                sys, op; eval_expression, eval_module,
                u0_constructor, cse, kwargs...
            )
            # Required for filling missing parameter values when this is an initialization
            # problem
            if state_values(linprob) === nothing
                linprob = remake(
                    linprob;
                    u0 = u0_constructor(ones(eltype(linprob.A), size(linprob.A, 2)))
                )
            end
            return SCCNonlinearProblem((linprob,), (Returns(nothing),), parameter_values(linprob), true; sys)
        else
            return NonlinearProblem{iip}(
                sys, op; eval_expression, eval_module, u0_constructor, cse, kwargs...
            )
        end
    end

    dvs = unknowns(sys)
    ps = parameters(sys)
    eqs = equations(sys)
    obs = observed(sys)
    unhacked_sys = unhack_system(sys)
    unhacked_obs = observed(unhacked_sys)

    _, u0, p = process_SciMLProblem(
        EmptySciMLFunction{iip}, sys, op; eval_expression, eval_module, symbolic_u0 = true,
        missing_guess_value, kwargs...
    )
    op = calculate_op_from_u0_p(sys, u0, p)

    explicitfuns = []
    nlfuns = []
    decomposition = SCCDecomposition(sys, var_sccs, eq_sccs)


    # Invalidate the SCC information - `decomposition` is the source of truth now
    var_sccs = nothing
    eq_sccs = nothing

    build_caches!(sys, decomposition)

    for i in eachindex(decomposition.subsystems)
        cachevars = decomposition.scc_cachevars[i]
        cacheexprs = decomposition.scc_cacheexprs[i]
        subsys = decomposition.subsystems[i]
        _prevobsidxs = decomposition.obsidxs_for_cacheexprs[i]
        if isempty(cachevars)
            push!(explicitfuns, Returns(nothing))
        else
            solsyms = view.((dvs,), view(decomposition.var_sccs, 1:(i - 1)))
            push!(
                explicitfuns,
                CacheWriter(
                    sys, decomposition.cachetypes, cacheexprs, solsyms, unhacked_obs[_prevobsidxs];
                    eval_expression, eval_module, cse
                )
            )
        end
        cachebufsyms = Vector{SymbolicT}[]
        for T in decomposition.cachetypes
            push!(cachebufsyms, get(cachevars, T, SymbolicT[]))
        end
        f = SCCNonlinearFunction{iip}(
            decomposition, i, cachebufsyms, op;
            eval_expression, eval_module, cse, kwargs...
        )
        push!(nlfuns, f)
    end

    u0_eltype = Union{}
    for x in u0
        symbolic_type(x) == NotSymbolic() || continue
        u0_eltype = typeof(x)
        break
    end
    if u0_eltype === Union{} || u0_eltype === Nothing
        u0_eltype = Float64
    end
    u0_eltype = float(u0_eltype)

    if !isempty(decomposition.cachetypes)
        templates = map(decomposition.cachetypes, decomposition.cachesizes) do T, n
            # Real refers to `eltype(u0)`
            if T == Real
                T = u0_eltype
            elseif T <: Array && eltype(T) == Real
                T = Array{u0_eltype, ndims(T)}
            end
            BufferTemplate(T, n)
        end
        p = rebuild_with_caches(p, templates...)
    end

    # yes, `get_p_constructor` since this is only used for `LinearProblem` and
    # will retain the shape of `A`
    u0_constructor = get_p_constructor(u0_constructor, typeof(u0), u0_eltype)
    subprobs = []
    subber = Symbolics.FixpointSubstituter{true}(AtomicArrayDictSubstitutionWrapper(op))
    for (i, (f, vscc)) in enumerate(zip(nlfuns, decomposition.var_sccs))
        _u0 = SymbolicUtils.Code.create_array(
            typeof(u0), eltype(u0), Val(1), Val(length(vscc)), u0[vscc]...
        )
        symbolic_idxs = findall(x -> x === nothing || symbolic_type(x) !== NotSymbolic(), _u0)
        if f isa LinearFunction
            _u0 = isempty(symbolic_idxs) ? _u0 : zeros(u0_eltype, length(_u0))
            _u0 = u0_constructor(u0_eltype.(_u0))
            cachevars = decomposition.scc_cachevars[i]
            cacheexprs = decomposition.scc_cacheexprs[i]
            for T in keys(cachevars)
                for (var, expr) in zip(cachevars[T], cacheexprs[T])
                    isequal(var, expr) && continue
                    has_possibly_indexed_key(op, var) && continue
                    write_possibly_indexed_array!(op, var, expr, COMMON_NOTHING)
                end
            end
            symbolic_interface = f.interface
            A, b = get_A_b_from_LinearFunction(
                sys, f, subber; eval_expression, eval_module, u0_constructor, u0_eltype
            )
            for (j, val) in zip(vscc, _u0)
                write_possibly_indexed_array!(op, dvs[j], Symbolics.SConst(val), COMMON_NOTHING)
            end
            prob = LinearProblem{iip}(A, b, p; f = symbolic_interface, u0 = _u0)
        else
            if !isempty(symbolic_idxs)
                Moshi.Match.@match missing_guess_value begin
                    MissingGuessValue.Constant(val) => begin
                        _u0[symbolic_idxs] .= val
                        _u0 = unwrap_const.(_u0)
                        cval = Symbolics.SConst(val)
                        for j in symbolic_idxs
                            write_possibly_indexed_array!(op, dvs[vscc[j]], cval, COMMON_NOTHING)
                        end
                    end
                    MissingGuessValue.Random(rng) => begin
                        newval = rand(rng, length(symbolic_idxs))
                        _u0[symbolic_idxs] .= newval
                        for (idx, j) in enumerate(symbolic_idxs)
                            write_possibly_indexed_array!(
                                op, dvs[vscc[j]], Symbolics.SConst(newval[idx]), COMMON_NOTHING
                            )
                        end
                    end
                    MissingGuessValue.Error() => throw(MissingGuessError(dvs[vscc], _u0))
                end
            end
            _u0 = u0_constructor(u0_eltype.(_u0))
            prob = NonlinearProblem(f, _u0, p)
        end
        push!(subprobs, prob)
    end

    new_dvs = dvs[reduce(vcat, decomposition.var_sccs)]
    new_eqs = eqs[reduce(vcat, decomposition.eq_sccs)]
    sys = ConstructionBase.setproperties(
        sys; unknowns = new_dvs, eqs = new_eqs, index_cache = subset_unknowns_observed(
            get_index_cache(sys), sys, new_dvs, getproperty.(obs, (:lhs,))
        )
    )
    return SCCNonlinearProblem(Tuple(subprobs), Tuple(explicitfuns), p, true; sys)
end

function calculate_op_from_u0_p(sys::System, u0::Union{Nothing, AbstractVector}, p::MTKParameters)
    op = SymmapT()
    if u0 !== nothing
        for (var, val) in zip(unknowns(sys), u0)
            val === nothing && continue
            write_possibly_indexed_array!(op, var, Symbolics.SConst(val), COMMON_NOTHING)
        end
    end
    rps = reorder_parameters(sys)
    @assert length(rps) == length(p)

    for (i, pvars) in enumerate(rps)
        for (var, val) in zip(pvars, p[i])
            write_possibly_indexed_array!(op, var, Symbolics.SConst(val), COMMON_NOTHING)
        end
    end

    _ss = unhack_system(sys)
    for eq in observed(_ss)
        write_possibly_indexed_array!(op, eq.lhs, eq.rhs, COMMON_NOTHING)
    end
    merge!(op, bindings(sys))
    return op
end
