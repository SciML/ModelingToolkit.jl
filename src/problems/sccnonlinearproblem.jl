const TypeT = Union{DataType, UnionAll}

struct CacheWriter{F}
    fn::F
end

function (cw::CacheWriter)(p, sols)
    cw.fn(p.caches, sols, p)
end

function CacheWriter(sys::AbstractSystem, buffer_types::Vector{TypeT},
        exprs::Dict{TypeT, Vector{Any}}, solsyms, obseqs::Vector{Equation};
        eval_expression = false, eval_module = @__MODULE__, cse = true, sparse = false)
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
        extra_assignments = [array_assignments; obs_assigns; body])
    fn = eval_or_rgf(fn; eval_expression, eval_module)
    fn = GeneratedFunctionWrapper{(3, 3, is_split(sys))}(fn, nothing)
    return CacheWriter(fn)
end

struct SCCNonlinearFunction{iip} end

function SCCNonlinearFunction{iip}(
        sys::System, _eqs, _dvs, _obs, cachesyms, op; eval_expression = false,
        eval_module = @__MODULE__, cse = true, kwargs...) where {iip}
    ps = parameters(sys; initial_parameters = true)
    subsys = System(
        _eqs, _dvs, ps; observed = _obs, name = nameof(sys), defaults = defaults(sys))
    @set! subsys.parameter_dependencies = parameter_dependencies(sys)
    if get_index_cache(sys) !== nothing
        @set! subsys.index_cache = subset_unknowns_observed(
            get_index_cache(sys), sys, _dvs, getproperty.(_obs, (:lhs,)))
        @set! subsys.complete = true
    end
    # generate linear problem instead
    if isaffine(subsys)
        return LinearFunction{iip}(
            subsys; eval_expression, eval_module, cse, cachesyms, kwargs...)
    end
    rps = reorder_parameters(sys, ps)

    obs_assignments = [eq.lhs ← eq.rhs for eq in _obs]

    rhss = [eq.rhs - eq.lhs for eq in _eqs]
    f_gen = build_function_wrapper(sys,
        rhss, _dvs, rps..., cachesyms...; p_start = 2,
        p_end = length(rps) + length(cachesyms) + 1, add_observed = false,
        extra_assignments = obs_assignments, expression = Val{true}, cse)
    f_oop, f_iip = eval_or_rgf.(f_gen; eval_expression, eval_module)
    f = GeneratedFunctionWrapper{(2, 2, is_split(sys))}(f_oop, f_iip)

    return NonlinearFunction{iip}(f; sys = subsys)
end

function SciMLBase.SCCNonlinearProblem(sys::System, args...; kwargs...)
    SCCNonlinearProblem{true}(sys, args...; kwargs...)
end

function SciMLBase.SCCNonlinearProblem{iip}(sys::System, op; eval_expression = false,
        eval_module = @__MODULE__, cse = true, u0_constructor = identity, kwargs...) where {iip}
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
            DiCMOBiGraph{true}(complete(ts.structure.graph),
                complete(var_eq_matching)),
            var_sccs)
        toporder = topological_sort_by_dfs(condensed_graph)
        var_sccs = var_sccs[toporder]
        eq_sccs = map(Base.Fix1(getindex, var_eq_matching), var_sccs)
    else
        var_sccs = sched.var_sccs
        # Equations are already in the order of SCCs
        eq_sccs = length.(var_sccs)
        cumsum!(eq_sccs, eq_sccs)
        eq_sccs = map(enumerate(eq_sccs)) do (i, lasti)
            i == 1 ? (1:lasti) : ((eq_sccs[i - 1] + 1):lasti)
        end
    end

    if length(var_sccs) == 1
        return NonlinearProblem{iip}(
            sys, op; eval_expression, eval_module, kwargs...)
    end

    dvs = unknowns(sys)
    ps = parameters(sys)
    eqs = equations(sys)
    obs = observed(sys)

    _, u0,
    p = process_SciMLProblem(
        EmptySciMLFunction{iip}, sys, op; eval_expression, eval_module, symbolic_u0 = true, kwargs...)

    explicitfuns = []
    nlfuns = []
    prevobsidxs = BlockArray(undef_blocks, Vector{Int}, Int[])
    # Cache buffer types and corresponding sizes. Stored as a pair of arrays instead of a
    # dict to maintain a consistent order of buffers across SCCs
    cachetypes = TypeT[]
    cachesizes = Int[]
    # explicitfun! related information for each SCC
    # We need to compute buffer sizes before doing any codegen
    scc_cachevars = Dict{TypeT, Vector{Any}}[]
    scc_cacheexprs = Dict{TypeT, Vector{Any}}[]
    scc_eqs = Vector{Equation}[]
    scc_obs = Vector{Equation}[]
    # variables solved in previous SCCs
    available_vars = Set()
    for (i, (escc, vscc)) in enumerate(zip(eq_sccs, var_sccs))
        # subset unknowns and equations
        _dvs = dvs[vscc]
        _eqs = eqs[escc]
        # get observed equations required by this SCC
        union!(available_vars, _dvs)
        obsidxs = observed_equations_used_by(sys, _eqs; available_vars)
        # the ones used by previous SCCs can be precomputed into the cache
        setdiff!(obsidxs, prevobsidxs)
        _obs = obs[obsidxs]
        union!(available_vars, getproperty.(_obs, (:lhs,)))

        # get all subexpressions in the RHS which we can precompute in the cache
        # precomputed subexpressions should not contain `banned_vars`
        banned_vars = Set{Any}(vcat(_dvs, getproperty.(_obs, (:lhs,))))
        state = Dict()
        for i in eachindex(_obs)
            _obs[i] = _obs[i].lhs ~ subexpressions_not_involving_vars!(
                _obs[i].rhs, banned_vars, state)
        end
        for i in eachindex(_eqs)
            _eqs[i] = _eqs[i].lhs ~ subexpressions_not_involving_vars!(
                _eqs[i].rhs, banned_vars, state)
        end

        # map from symtype to cached variables and their expressions
        cachevars = Dict{Union{DataType, UnionAll}, Vector{Any}}()
        cacheexprs = Dict{Union{DataType, UnionAll}, Vector{Any}}()
        # observed of previous SCCs are in the cache
        # NOTE: When we get proper CSE, we can substitute these
        # and then use `subexpressions_not_involving_vars!`
        for i in prevobsidxs
            T = symtype(obs[i].lhs)
            buf = get!(() -> Any[], cachevars, T)
            push!(buf, obs[i].lhs)

            buf = get!(() -> Any[], cacheexprs, T)
            push!(buf, obs[i].lhs)
        end

        for (k, v) in state
            k = unwrap(k)
            v = unwrap(v)
            T = symtype(k)
            buf = get!(() -> Any[], cachevars, T)
            push!(buf, v)
            buf = get!(() -> Any[], cacheexprs, T)
            push!(buf, k)
        end

        # update the sizes of cache buffers
        for (T, buf) in cachevars
            idx = findfirst(isequal(T), cachetypes)
            if idx === nothing
                push!(cachetypes, T)
                push!(cachesizes, 0)
                idx = lastindex(cachetypes)
            end
            cachesizes[idx] = max(cachesizes[idx], length(buf))
        end

        push!(scc_cachevars, cachevars)
        push!(scc_cacheexprs, cacheexprs)
        push!(scc_eqs, _eqs)
        push!(scc_obs, _obs)
        blockpush!(prevobsidxs, obsidxs)
    end

    for (i, (escc, vscc)) in enumerate(zip(eq_sccs, var_sccs))
        _dvs = dvs[vscc]
        _eqs = scc_eqs[i]
        _prevobsidxs = reduce(vcat, blocks(prevobsidxs)[1:(i - 1)]; init = Int[])
        _obs = scc_obs[i]
        cachevars = scc_cachevars[i]
        cacheexprs = scc_cacheexprs[i]
        available_vars = [dvs[reduce(vcat, var_sccs[1:(i - 1)]; init = Int[])];
                          getproperty.(
                              reduce(vcat, scc_obs[1:(i - 1)]; init = []), (:lhs,))]
        _prevobsidxs = vcat(_prevobsidxs,
            observed_equations_used_by(
                sys, reduce(vcat, values(cacheexprs); init = []); available_vars))
        if isempty(cachevars)
            push!(explicitfuns, Returns(nothing))
        else
            solsyms = getindex.((dvs,), view(var_sccs, 1:(i - 1)))
            push!(explicitfuns,
                CacheWriter(sys, cachetypes, cacheexprs, solsyms, obs[_prevobsidxs];
                    eval_expression, eval_module, cse))
        end

        cachebufsyms = Tuple(map(cachetypes) do T
            get(cachevars, T, [])
        end)
        f = SCCNonlinearFunction{iip}(
            sys, _eqs, _dvs, _obs, cachebufsyms, op;
            eval_expression, eval_module, cse, kwargs...)
        push!(nlfuns, f)
    end

    u0_eltype = Union{}
    for x in u0
        symbolic_type(x) == NotSymbolic() || continue
        u0_eltype = typeof(x)
        break
    end
    if u0_eltype == Union{}
        u0_eltype = Float64
    end
    u0_eltype = float(u0_eltype)

    if !isempty(cachetypes)
        templates = map(cachetypes, cachesizes) do T, n
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
    for (i, (f, vscc)) in enumerate(zip(nlfuns, var_sccs))
        _u0 = SymbolicUtils.Code.create_array(
            typeof(u0), eltype(u0), Val(1), Val(length(vscc)), u0[vscc]...)
        symbolic_idxs = findall(x -> symbolic_type(x) != NotSymbolic(), _u0)
        if f isa LinearFunction
            _u0 = isempty(symbolic_idxs) ? _u0 : zeros(u0_eltype, length(_u0))
            _u0 = u0_eltype.(_u0)
            symbolic_interface = f.interface
            A,
            b = get_A_b_from_LinearFunction(
                sys, f, p; eval_expression, eval_module, u0_constructor, u0_eltype)
            prob = LinearProblem{iip}(A, b, p; f = symbolic_interface, u0 = _u0)
        else
            isempty(symbolic_idxs) || throw(MissingGuessError(dvs[vscc], _u0))
            _u0 = u0_eltype.(_u0)
            prob = NonlinearProblem(f, _u0, p)
        end
        push!(subprobs, prob)
    end

    new_dvs = dvs[reduce(vcat, var_sccs)]
    new_eqs = eqs[reduce(vcat, eq_sccs)]
    @set! sys.unknowns = new_dvs
    @set! sys.eqs = new_eqs
    @set! sys.index_cache = subset_unknowns_observed(
        get_index_cache(sys), sys, new_dvs, getproperty.(obs, (:lhs,)))
    return SCCNonlinearProblem(Tuple(subprobs), Tuple(explicitfuns), p, true; sys)
end
