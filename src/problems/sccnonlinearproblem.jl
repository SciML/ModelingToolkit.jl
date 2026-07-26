struct CacheWriter{F}
    fn::F
end

function (cw::CacheWriter)(p::MTKParameters, sols)
    return cw.fn(p.caches, sols, p)
end

const SCCCacheVarsExprsElT = Dict{TypeT, Vector{SymbolicT}}

const SCC_EXPLICITFUN_CACHE_OUT = unwrap(only(@parameters __outₘₜₖ::Vector{Vector{Any}}))

function CacheWriter(
        sys::AbstractSystem, buffer_types::Vector{TypeT},
        exprs::SCCCacheVarsExprsElT, solsyms, opts::GeneratedFunctionOptions;
        sparse = false
    )
    (; eval_expression, eval_module) = opts
    rps = reorder_parameters(sys)  # 1 arg to use the cached version
    cache_writes = SymbolicT[]
    for (i, T) in enumerate(buffer_types)
        regions = SU.RegionsT()
        values = Symbolics.SArgsT()
        output = SCC_EXPLICITFUN_CACHE_OUT[i]
        cacheexprs = get(exprs, T, SymbolicT[])
        isempty(cacheexprs) && continue
        N = length(cacheexprs)
        allocator = Symbolics.STerm(
            Returns, Symbolics.SArgsT((output,));
            type = SU.FnType{Tuple, Vector{T}, Any}, shape = SU.ShapeVecT((1:N,))
        )
        for (j, expr) in enumerate(cacheexprs)
            push!(regions, SU.ShapeVecT((j:j,)))
            push!(values, Symbolics.SConst([expr]))
        end
        maker = SU.ArrayMaker{VartypeT}(regions, values; shape = SU.ShapeVecT((1:N,)))
        writer = Code.with_allocator(allocator, maker)
        push!(cache_writes, writer)
    end
    body = Symbolics.STerm(
        tuple, cache_writes;
        type = Vector{Any}, shape = SU.ShapeVecT((1:length(cache_writes),))
    )

    fn, _ = build_function_wrapper(
        sys, body, [Any[SCC_EXPLICITFUN_CACHE_OUT]; solsyms; rps],
        BuildFunctionWrapperOptions(;
            p_start = length(solsyms) + 2, p_end = length(rps) + length(solsyms) + 1,
            compress_args = [2:(length(solsyms) + 1)],
            codegen_function_options = ConstructionBase.setproperties(
                opts.codegen, (; iip_config = (true, false))
            )
        )
    )
    fn = eval_or_rgf(fn; eval_expression, eval_module)
    fn = GeneratedFunctionWrapper{(3, 3, is_split(sys))}(fn, nothing)
    return CacheWriter{Any}(fn)
end

# Backward-compatibility keyword method. The positional `opts::GeneratedFunctionOptions`
# method above is the primary; this wrapper preserves the historical keyword API.
function CacheWriter(
        sys::AbstractSystem, buffer_types::Vector{TypeT},
        exprs::SCCCacheVarsExprsElT, solsyms;
        eval_expression = false, eval_module = @__MODULE__, sparse = false
    )
    return CacheWriter(
        sys, buffer_types, exprs, solsyms,
        GeneratedFunctionOptions(; eval_expression, eval_module); sparse
    )
end

# This phrasing allows us to precompile the calls
@noinline function __explicitfun_copy_states_helper(buffer, sols)
    offset = 0
    for sol in sols
        u = sol.u
        copyto!(buffer, CartesianIndices(((offset + 1):(offset + length(u)),)), u, CartesianIndices((1:length(u),)))
        offset += length(u)
    end
    return nothing
end

function __explicitfun_copy_states(@nospecialize(p), sols)
    buffer_idx = findfirst(Base.Fix2(==, eltype(sols[1].u)) ∘ eltype, p.caches)::Int
    buffer = p.caches[buffer_idx]::Vector{eltype(sols[1].u)}
    __explicitfun_copy_states_helper(buffer, sols)
    return nothing
end

"""
    $TYPEDSIGNATURES

Subset a system to have only the given unknowns `vscc` and equations `escc`. Observed
equations are subset, delete accordingly. Requires that `sys` is complete and flattened.

# Keyword arguments

- `available_vars`: A list of variables that the subset system should assume are precomputed
  or already available. Will be mutated with the unknowns and observables of the subset. This
  is useful for SCC decomposition.
"""
function subset_system(
        sys::System, vscc::Vector{Int}, escc::Vector{Int};
        available_vars = Set{SymbolicT}()
    )
    check_complete(sys, "subset_system")
    @assert isempty(get_systems(sys)) "`subset_system` requires a flattened system"

    dvs = unknowns(sys)
    ps = parameters(sys)
    eqs = full_equations(sys)

    # subset unknowns and equations
    _dvs = dvs[vscc]
    _eqs = eqs[escc]
    union!(available_vars, _dvs)

    subsys = ConstructionBase.setproperties(
        sys; unknowns = _dvs, eqs = _eqs, observed = Equation[],
        parameter_bindings_graph = get_parameter_bindings_graph(sys), complete = true
    )
    if get_index_cache(sys) !== nothing
        @set! subsys.index_cache = subset_unknowns_observed(
            get_index_cache(sys), sys, _dvs, SymbolicT[],
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

struct SCCDecomposition
    subsystems::Vector{System}
    var_sccs::Vector{Vector{Int}}
    eq_sccs::Vector{Vector{Int}}
    islinear::BitVector
    hints::Vector{StructuralHint.Type}
    # Cache buffer types and corresponding sizes. Stored as a pair of arrays instead of a
    # dict to maintain a consistent order of buffers across SCCs
    cachetypes::Vector{TypeT}
    cachesizes::Vector{Int}
    # explicitfun! related information for each SCC
    # We need to compute buffer sizes before doing any codegen
    scc_cachevars::Vector{SCCCacheVarsExprsElT}
    scc_cacheexprs::Vector{SCCCacheVarsExprsElT}
end

function SCCDecomposition()
    return SCCDecomposition(
        System[], Vector{Int}[], Vector{Int}[], BitVector(), StructuralHint.Type[],
        TypeT[], Int[], SCCCacheVarsExprsElT[], SCCCacheVarsExprsElT[]
    )
end

function SCCDecomposition(
        sys::System, var_sccs::Vector{Vector{Int}}, eq_sccs::Vector{Vector{Int}};
        combine_sccs = true
    )
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
            finalize_scc!(
                final_decomposition, active_decomposition, first(nbors), active, nbors
            )
        end
        subsys = subset_system(
            sys, vscc, escc; available_vars
        )
        push!(active_decomposition.subsystems, subsys)
        push!(active_decomposition.var_sccs, vscc)
        push!(active_decomposition.eq_sccs, escc)
        push!(active_decomposition.hints, StructuralHint.NoHint())
        push!(active_decomposition.islinear, calculate_A_b(subsys; throw = false) !== nothing)
        push!(active, i)
        if !combine_sccs
            finalize_scc!(
                final_decomposition, active_decomposition, first(active), active, active
            )
        end
    end

    while !isempty(active)
        finalize_scc!(final_decomposition, active_decomposition, first(active), active, active)
    end

    return final_decomposition
end

function copy_scc!(dst::SCCDecomposition, src::SCCDecomposition, tgt::Int)
    push!(dst.subsystems, src.subsystems[tgt])
    push!(dst.var_sccs, src.var_sccs[tgt])
    push!(dst.eq_sccs, src.eq_sccs[tgt])
    push!(dst.hints, src.hints[tgt])
    push!(dst.islinear, src.islinear[tgt])
    return dst
end

function finalize_scc!(
        final_decomposition::SCCDecomposition, active_decomposition::SCCDecomposition, i::Int,
        active::Set{Int}, nbors::Set{Int}; linear_scc_combine_range::Int = 2
    )
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
        # The way we run `subset_system`, if a later SCC (later index in `.subsystems`)
        # needs an observed from a previous SCC, that equation won't be present in `obsidxs`
        # and it will rely on `explicitfun!` for the required value. This means that if we
        # don't merge SCCs in sorted order, the observed equations of the resultant bigger
        # system won't be topologically sorted.
        sort!(to_merge)
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
        # See the note in the diagonal SCC case above to know why sorting is necessary
        sort!(to_merge)
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
    state = Dict{SymbolicT, SymbolicT}()
    ir = get_irstructure(sys)
    previous_scc_unknowns = SymbolicT[]
    for i in eachindex(decomposition.subsystems)
        subsys = decomposition.subsystems[i]
        if decomposition.islinear[i]
            store_to_mutable_cache!(subsys, CachedLinearAb, nothing)
            # For linear systems, there is no point in a complicated explicitfun. Just copy
            # all previous states to the cache buffer.

            cachevars = SCCCacheVarsExprsElT()
            cacheexprs = SCCCacheVarsExprsElT()
            push!(decomposition.scc_cachevars, cachevars)
            push!(decomposition.scc_cacheexprs, cacheexprs)

            if isempty(previous_scc_unknowns)
                append!(previous_scc_unknowns, unknowns(subsys))
                continue
            end
            T = Real
            buffer = cachevars[T] = cacheexprs[T] = SymbolicT[]
            append!(buffer, previous_scc_unknowns)
            idx = findfirst(isequal(T), decomposition.cachetypes)
            if idx === nothing
                push!(decomposition.cachetypes, T)
                push!(decomposition.cachesizes, 0)
                idx = lastindex(decomposition.cachetypes)
            end
            decomposition.cachesizes[idx] = max(decomposition.cachesizes[idx], length(buffer))
            append!(previous_scc_unknowns, unknowns(subsys))
            continue
        end
        empty!(banned_vars)
        empty!(state)

        union!(banned_vars, unknowns(subsys))
        append!(previous_scc_unknowns, unknowns(subsys))
        for u in unknowns(subsys)
            push!(banned_vars, split_indexed_var(u)[1])
        end

        # While we own the system and so mutation _should_ be safe, `IRInfo` exists.
        # It stores the indices corresponding to equations in the `IRStructure`, which
        # would be incorrect if we mutate.
        _eqs = copy(get_eqs(subsys))
        exprs_to_search = SymbolicT[]
        for i in eachindex(_eqs)
            push!(exprs_to_search, _eqs[i].rhs)
        end
        subexpressions_not_involving_vars!(ir, exprs_to_search, banned_vars, state)
        subber = SU.IRSubstituter{false}(ir, state; filterer = !SU.default_is_atomic)
        for i in eachindex(_eqs)
            _eqs[i] = _eqs[i].lhs ~ subber(_eqs[i].rhs)
        end
        subsys = decomposition.subsystems[i] = ConstructionBase.setproperties(subsys; eqs = _eqs)

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

        # Iterate `state` in canonical (hash/objectid-independent) order so the cache
        # buffers, the generated cache-writer code and `cachetypes` are deterministic
        # rather than ordered by the `Dict` hash of the (symbolic) keys / (type) buckets.
        for (k, v) in sort!(collect(state); by = p -> MTKBase.canonical_sort_key(unwrap(first(p))))
            k = unwrap(k)
            v = unwrap(v)
            T = symtype(k)
            buf = get!(() -> SymbolicT[], cachevars, T)
            push!(buf, v)
            buf = get!(() -> SymbolicT[], cacheexprs, T)
            push!(buf, k)
        end
        sorted_cache_types = sort!(collect(keys(cachevars)); by = string)
        all_cacheexprs = reduce(
            vcat, (cacheexprs[T] for T in sorted_cache_types); init = SymbolicT[])
        # update the sizes of cache buffers
        for T in sorted_cache_types
            buf = cachevars[T]
            idx = findfirst(isequal(T), decomposition.cachetypes)
            if idx === nothing
                push!(decomposition.cachetypes, T)
                push!(decomposition.cachesizes, 0)
                idx = lastindex(decomposition.cachetypes)
            end
            decomposition.cachesizes[idx] = max(decomposition.cachesizes[idx], length(buf))
        end
    end
    return
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
    new_parent = decomposition.subsystems[i] = ConstructionBase.setproperties(
        parent; eqs = new_eqs, unknowns = new_dvs,
    )
    if cached_ab isa CachedLinearAb
        store_to_mutable_cache!(new_parent, CachedLinearAb, cached_ab)
    end

    return nothing
end

struct SCCNonlinearFunction{iip} end

function SCCNonlinearFunction{iip}(
        decomposition::SCCDecomposition, i::Int, cachesyms, op; eval_expression = false,
        eval_module = @__MODULE__, kwargs...
    ) where {iip}
    subsys = decomposition.subsystems[i]
    islin = decomposition.islinear[i]
    # generate linear problem instead
    if islin
        return LinearFunction{iip}(
            subsys; eval_expression, eval_module, cachesyms,
            structural_hint = decomposition.hints[i], kwargs...
        )
    end
    rps = reorder_parameters(subsys)

    if MTKBase.has_any_homotopy(subsys)
        # This block carries Modelica `homotopy(actual, simplified)` nodes, so compile the
        # λ-swept residual `f(u, p, λ)` instead of the opaque-`actual` one: the block is
        # built as a `SciMLBase.HomotopyProblem` below and solved by continuation. `λ` is a
        # trailing argument, never a parameter, so the block's cache buffers (which are
        # appended to the parameter list) are unaffected.
        shadow, λ = MTKBase.lower_homotopy(subsys)
        hf = MTKBase.generate_homotopy_residual(
            shadow, λ; expression = Val{false}, wrap_gfw = Val{true},
            eval_expression, eval_module, cachesyms
        )
        return NonlinearFunction{iip}(hf; sys = subsys)
    end

    f = generate_rhs(
        subsys,
        GeneratedFunctionOptions(;
            expression = Val{false}, wrap_gfw = Val{true}, eval_expression, eval_module
        );
        cachesyms
    )

    return NonlinearFunction{iip}(f; sys = subsys)
end

function SciMLBase.SCCNonlinearProblem(sys::System, op; kwargs...)
    return SCCNonlinearProblem{true}(sys, op; kwargs...)
end

function SciMLBase.SCCNonlinearProblem{iip}(
        sys::System, op; eval_expression = false,
        eval_module = @__MODULE__, u0_constructor = identity,
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
        var_sccs = map(copy, sched.var_sccs)
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
                u0_constructor, kwargs...
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
            # A single nonlinear block is solved directly rather than wrapped in an
            # `SCCNonlinearProblem`. When it carries Modelica `homotopy(actual, simplified)`
            # nodes this must be a `HomotopyProblem` so the block is solved by continuation
            # (the multi-block path below already builds `HomotopyProblem` blocks); a plain
            # `NonlinearProblem` would drop the λ-sweep and Newton-solve the target directly.
            TProb = MTKBase.get_nonlinear_problem_type(sys)
            return TProb{iip}(
                sys, op; eval_expression, eval_module, u0_constructor, missing_guess_value, kwargs...
            )
        end
    end

    dvs = unknowns(sys)
    ps = parameters(sys)
    eqs = equations(sys)

    _, u0, p = process_SciMLProblem(
        EmptySciMLFunction{iip}, sys, op; eval_expression, eval_module, symbolic_u0 = true,
        missing_guess_value, kwargs...
    )
    op = calculate_op_from_u0_p(sys, u0, p)

    explicitfuns = Union{Returns{Nothing}, CacheWriter{Any}, typeof(__explicitfun_copy_states)}[]
    nlfuns = []
    decomposition = SCCDecomposition(sys, var_sccs, eq_sccs; combine_sccs)


    # Invalidate the SCC information - `decomposition` is the source of truth now
    var_sccs = nothing
    eq_sccs = nothing

    build_caches!(sys, decomposition)

    # Every SCC subsystem is derived from `sys` via `subset_unknowns_observed`, which shares
    # the entire parameter portion of the index cache by reference, and none of the subsystem
    # constructions (`subset_system`, `_collapse_into!`) touch `ps`. The default parameter
    # reordering is therefore identical between `sys` and every (merged or unmerged) subsystem.
    if get_index_cache(sys) !== nothing
        rp = reorder_parameters(sys)
        cached_reorder = check_mutable_cache(
            sys, MTKBase.ReorderedDefaultParameters, MTKBase.ReorderedDefaultParameters, nothing
        )
        # The parameter array decomposition (`ParameterArrayAssignments`) is likewise a pure
        # function of the shared parameter layout. Compute it once (param-only) and propagate
        # to every subsystem; `build_function_wrapper` decomposes each subsystem's small
        # cachesym tail on top of it (see the `n_param_buffers` handling there).
        param_paa = MTKBase.ParameterArrayAssignments(
            MTKBase.compute_array_variable_buffer_idxs(rp)
        )
        for subsys in decomposition.subsystems
            if cached_reorder isa MTKBase.ReorderedDefaultParameters
                store_to_mutable_cache!(
                    subsys, MTKBase.ReorderedDefaultParameters, cached_reorder
                )
            end
            store_to_mutable_cache!(subsys, MTKBase.ParameterArrayAssignments, param_paa)
        end
    end

    for i in eachindex(decomposition.subsystems)
        cachevars = decomposition.scc_cachevars[i]
        cacheexprs = decomposition.scc_cacheexprs[i]
        subsys = decomposition.subsystems[i]
        if isempty(cachevars)
            push!(explicitfuns, Returns(nothing))
        elseif decomposition.islinear[i]
            push!(explicitfuns, __explicitfun_copy_states)
        else
            solsyms = view.((dvs,), view(decomposition.var_sccs, 1:(i - 1)))
            push!(
                explicitfuns,
                CacheWriter(
                    sys, decomposition.cachetypes, cacheexprs, solsyms,
                    GeneratedFunctionOptions(; eval_expression, eval_module)
                )
            )
        end
        cachebufsyms = Vector{SymbolicT}[]
        for T in decomposition.cachetypes
            push!(cachebufsyms, get(cachevars, T, SymbolicT[]))
        end
        f = SCCNonlinearFunction{iip}(
            decomposition, i, cachebufsyms, op;
            eval_expression, eval_module, kwargs...
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
            typeof(u0), Any, Val(1), Val(length(vscc)), u0[vscc]...
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
            symbolic_interface = MTKBase.wrap_symbolic_linear_interface(symbolic_interface, iip, A, b, p)
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
                    MissingGuessValue.HashedRandom() => begin
                        newval = [hash(dvs[vscc[j]]) for j in symbolic_idxs] ./ 0x1p64
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
            prob = if MTKBase.has_any_homotopy(decomposition.subsystems[i])
                # `f` is the λ-swept residual for this block (see `SCCNonlinearFunction`);
                # solving the block sweeps λ from the `simplified` form to `actual`. Blocks
                # are solved in dependency order, so each one reaches `λ = 1` before the
                # next begins.
                SciMLBase.HomotopyProblem(f, _u0, p; λspan = (0.0, 1.0))
            else
                NonlinearProblem(f, _u0, p)
            end
        end
        push!(subprobs, prob)
    end

    new_dvs = dvs[reduce(vcat, decomposition.var_sccs)]
    new_eqs = eqs[reduce(vcat, decomposition.eq_sccs)]
    sys = ConstructionBase.setproperties(
        sys; unknowns = new_dvs, eqs = new_eqs, index_cache = subset_unknowns_observed(
            get_index_cache(sys), sys, new_dvs, SymbolicT[]
        )
    )

    if length(subprobs) <= 5
        return SCCNonlinearProblem(Tuple(subprobs), Tuple(explicitfuns), p, true; sys)
    else
        return SCCNonlinearProblem(subprobs, SciMLBase.Void{Any}.(explicitfuns), p, true; sys)
    end
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

    _ss = MTKBase.reverse_all_default_reversible_transformations(sys)
    for eq in observed(_ss)
        write_possibly_indexed_array!(op, eq.lhs, eq.rhs, COMMON_NOTHING)
    end
    merge!(op, bindings(sys))
    return op
end
