using LinearAlgebra

using ModelingToolkit: isdifferenceeq, has_continuous_events, generate_rootfinding_callback, generate_difference_cb, merge_cb

const MAX_INLINE_NLSOLVE_SIZE = 8

function torn_system_jacobian_sparsity(sys, var_eq_matching, var_sccs, nlsolve_scc_idxs, eqs_idxs, states_idxs)
    s = structure(sys)
    @unpack fullvars, graph = s

    # The sparsity pattern of `nlsolve(f, u, p)` w.r.t `p` is difficult to
    # determine in general. Consider the "simplest" case, a linear system. We
    # have
    #                   A u = p.
    # Clearly, the sparsity of `u` depends on the sparsity of both `p` and `A`
    # in a non-trivial way. However, in the generic case, `u` is dense even when
    # `A` and `p` are sparse. For instance
    #
    # ```julia
    # julia> using Random, SparseArrays
    #
    # julia> A = sprand(MersenneTwister(1234), 100, 100, 0.1);
    #
    # julia> p = sprand(MersenneTwister(12345), 100, 0.05);
    #
    # julia> count(x->abs(x) < 1e-5, A \ Vector(p))
    # 0
    # ```
    #
    # Let 𝑇 be the set of tearing variables and 𝑉 be the set of all *states* in
    # the residual equations. In the following code, we are going to assume the
    # connection between 𝑇 (the `u` in from above) and 𝑉 ∖ 𝑇 (the `p` in from
    # above) has full incidence.
    #
    # Note that as we are reducing algebraic equations numerically, it could be
    # the case that a later partition (a BLT block) contains tearing variables
    # from other partitions.
    #
    # We know that partitions are BLT ordered. Hence, the tearing variables in
    # each partition is unique, and all states in a partition must be
    # either differential variables or algebraic tearing variables that are
    # from previous partitions. Hence, we can build the dependency chain as we
    # traverse the partitions.

    var_rename = ones(Int64, ndsts(graph))
    nlsolve_vars = Int[]
    for i in nlsolve_scc_idxs, c in var_sccs[i]
        append!(nlsolve_vars, c)
        for v in c
            var_rename[v] = 0
        end
    end
    masked_cumsum!(var_rename)

    dig = DiCMOBiGraph{true}(graph, var_eq_matching)

    fused_var_deps = map(1:ndsts(graph)) do v
        BitSet(v′ for v′ in neighborhood(dig, v, Inf; dir=:in) if var_rename[v′] != 0)
    end

    for scc in var_sccs[nlsolve_scc_idxs]
        if length(scc) >= 2
            deps = fused_var_deps[scc[1]]
            for c in 2:length(scc)
                union!(deps, fused_var_deps[c])
                fused_var_deps[c] = deps
            end
        end
    end

    var2idx = Dict{Int,Int}(v => i for (i, v) in enumerate(states_idxs))
    eqs2idx = Dict{Int,Int}(v => i for (i, v) in enumerate(eqs_idxs))
    nlsolve_vars_set = BitSet(nlsolve_vars)

    I = Int[]; J = Int[]
    for ieq in 𝑠vertices(graph)
        nieq = get(eqs2idx, ieq, 0)
        nieq == 0 && continue
        for ivar in 𝑠neighbors(graph, ieq)
            isdervar(s, ivar) && continue
            if var_rename[ivar] != 0
                push!(I, nieq)
                push!(J, var2idx[ivar])
            else
                for dvar in fused_var_deps[ivar]
                    isdervar(s, dvar) && continue
                    niv = get(var2idx, dvar, 0)
                    niv == 0 && continue
                    push!(I, nieq)
                    push!(J, niv)
                end
            end
        end
    end
    sparse(I, J, true)
end

function gen_nlsolve!(is_not_prepended_assignment, eqs, vars, u0map::AbstractDict, assignments, (deps, invdeps), var2assignment; checkbounds=true)
    isempty(vars) && throw(ArgumentError("vars may not be empty"))
    length(eqs) == length(vars) || throw(ArgumentError("vars must be of the same length as the number of equations to find the roots of"))
    rhss = map(x->x.rhs, eqs)
    # We use `vars` instead of `graph` to capture parameters, too.
    paramset = ModelingToolkit.vars(r for r in rhss)

    # Compute necessary assignments for the nlsolve expr
    init_assignments = [var2assignment[p] for p in paramset if haskey(var2assignment, p)]
    tmp = [init_assignments]
    # `deps[init_assignments]` gives the dependency of `init_assignments`
    while true
        next_assignments = reduce(vcat, deps[init_assignments])
        isempty(next_assignments) && break
        init_assignments = next_assignments
        push!(tmp, init_assignments)
    end
    needed_assignments_idxs = reduce(vcat, unique(reverse(tmp)))
    needed_assignments = assignments[needed_assignments_idxs]

    # Compute `params`. They are like enclosed variables
    rhsvars = [ModelingToolkit.vars(r.rhs) for r in needed_assignments]
    vars_set = Set(vars)
    outer_set = BitSet()
    inner_set = BitSet()
    for (i, vs) in enumerate(rhsvars)
        j = needed_assignments_idxs[i]
        if isdisjoint(vars_set, vs)
            push!(outer_set, j)
        else
            push!(inner_set, j)
        end
    end
    init_refine = BitSet()
    for i in inner_set
        union!(init_refine, invdeps[i])
    end
    intersect!(init_refine, outer_set)
    setdiff!(outer_set, init_refine)
    union!(inner_set, init_refine)

    next_refine = BitSet()
    while true
        for i in init_refine
            id = invdeps[i]
            isempty(id) && break
            union!(next_refine, id)
        end
        intersect!(next_refine, outer_set)
        isempty(next_refine) && break
        setdiff!(outer_set, next_refine)
        union!(inner_set, next_refine)

        init_refine, next_refine = next_refine, init_refine
        empty!(next_refine)
    end
    global2local = Dict(j=>i for (i, j) in enumerate(needed_assignments_idxs))
    inner_idxs = [global2local[i] for i in collect(inner_set)]
    outer_idxs = [global2local[i] for i in collect(outer_set)]
    extravars = reduce(union!, rhsvars[inner_idxs], init=Set())
    union!(paramset, extravars)
    setdiff!(paramset, vars)
    setdiff!(paramset, [needed_assignments[i].lhs for i in inner_idxs])
    union!(paramset, [needed_assignments[i].lhs for i in outer_idxs])
    params = collect(paramset)

    # splatting to tighten the type
    u0 = []
    for v in vars
        v in keys(u0map) || (push!(u0, 1e-3); continue)
        u = substitute(v, u0map)
        for i in 1:length(u0map)
            u = substitute(u, u0map)
            u isa Number && (push!(u0, u); break)
        end
        u isa Number || error("$v doesn't have a default.")
    end
    u0 = [u0...]
    # specialize on the scalar case
    isscalar = length(u0) == 1
    u0 = isscalar ? u0[1] : SVector(u0...)

    fname = gensym("fun")
    # f is the function to find roots on
    f = Func(
        [
         DestructuredArgs(vars, inbounds=!checkbounds)
         DestructuredArgs(params, inbounds=!checkbounds)
        ],
        [],
        Let(
            needed_assignments[inner_idxs],
            isscalar ? rhss[1] : MakeArray(rhss, SVector)
           )
    ) |> SymbolicUtils.Code.toexpr

    # solver call contains code to call the root-finding solver on the function f
    solver_call = LiteralExpr(quote
            $numerical_nlsolve(
                               $fname,
                               # initial guess
                               $u0,
                               # "captured variables"
                               ($(params...),)
                              )
        end)

    preassignments = []
    for i in outer_idxs
        ii = needed_assignments_idxs[i]
        is_not_prepended_assignment[ii] || continue
        is_not_prepended_assignment[ii] = false
        push!(preassignments, assignments[ii])
    end

    nlsolve_expr = Assignment[
                    preassignments
                    fname ← @RuntimeGeneratedFunction(f)
                    DestructuredArgs(vars, inbounds=!checkbounds) ← solver_call
                   ]

    nlsolve_expr
end

function build_torn_function(
        sys;
        expression=false,
        jacobian_sparsity=true,
        checkbounds=false,
        max_inlining_size=nothing,
        kw...
    )

    max_inlining_size = something(max_inlining_size, MAX_INLINE_NLSOLVE_SIZE)
    rhss = []
    eqs = equations(sys)
    eqs_idxs = Int[]
    for (i, eq) in enumerate(eqs)
        isdiffeq(eq) || continue
        push!(eqs_idxs, i)
        push!(rhss, eq.rhs)
    end

    s = structure(sys)
    @unpack fullvars = s
    var_eq_matching, var_sccs = algebraic_variables_scc(sys)
    condensed_graph = MatchedCondensationGraph(
        DiCMOBiGraph{true}(complete(s.graph), complete(var_eq_matching)), var_sccs)
    toporder = topological_sort_by_dfs(condensed_graph)
    var_sccs = var_sccs[toporder]

    states_idxs = collect(diffvars_range(s))
    mass_matrix_diag = ones(length(states_idxs))

    assignments, deps, sol_states = tearing_assignments(sys)
    invdeps = map(_->BitSet(), deps)
    for (i, d) in enumerate(deps)
        for a in d
            push!(invdeps[a], i)
        end
    end
    var2assignment = Dict{Any,Int}(eq.lhs => i for (i, eq) in enumerate(assignments))
    is_not_prepended_assignment = trues(length(assignments))

    torn_expr = Assignment[]

    defs = defaults(sys)
    nlsolve_scc_idxs = Int[]

    needs_extending = false
    @views for (i, scc) in enumerate(var_sccs)
        torn_vars_idxs = Int[var for var in scc if var_eq_matching[var] !== unassigned]
        torn_eqs_idxs = [var_eq_matching[var] for var in torn_vars_idxs]
        isempty(torn_eqs_idxs) && continue
        if length(torn_eqs_idxs) <= max_inlining_size
            nlsolve_expr = gen_nlsolve!(is_not_prepended_assignment, eqs[torn_eqs_idxs], s.fullvars[torn_vars_idxs], defs, assignments, (deps, invdeps), var2assignment, checkbounds=checkbounds)
            append!(torn_expr, nlsolve_expr)
            push!(nlsolve_scc_idxs, i)
        else
            needs_extending = true
            append!(eqs_idxs, torn_eqs_idxs)
            append!(rhss, map(x->x.rhs, eqs[torn_eqs_idxs]))
            append!(states_idxs, torn_vars_idxs)
            append!(mass_matrix_diag, zeros(length(torn_eqs_idxs)))
        end
    end

    mass_matrix = needs_extending ? Diagonal(mass_matrix_diag) : I

    out = Sym{Any}(gensym("out"))
    funbody = SetArray(
        !checkbounds,
        out,
        rhss
    )

    states = s.fullvars[states_idxs]
    syms = map(Symbol, states_idxs)

    pre = get_postprocess_fbody(sys)

    expr = SymbolicUtils.Code.toexpr(
        Func(
             [
              out
              DestructuredArgs(states, inbounds=!checkbounds)
              DestructuredArgs(parameters(sys), inbounds=!checkbounds)
              independent_variables(sys)
             ],
             [],
             pre(Let(
                     [torn_expr; assignments[is_not_prepended_assignment]],
                 funbody
                ))
            ),
        sol_states
    )
    if expression
        expr, states
    else
        observedfun = let sys=sys, dict=Dict(), assignments=assignments, deps=(deps, invdeps), sol_states=sol_states, var2assignment=var2assignment
            function generated_observed(obsvar, u, p, t)
                obs = get!(dict, value(obsvar)) do
                    build_observed_function(sys, obsvar, var_eq_matching, var_sccs,
                                            assignments, deps, sol_states, var2assignment,
                                            checkbounds=checkbounds,
                                           )
                end
                obs(u, p, t)
            end
        end

        ODEFunction{true}(
                          @RuntimeGeneratedFunction(expr),
                          sparsity = jacobian_sparsity ? torn_system_jacobian_sparsity(sys, var_eq_matching, var_sccs, nlsolve_scc_idxs, eqs_idxs, states_idxs) : nothing,
                          syms = syms,
                          observed = observedfun,
                          mass_matrix = mass_matrix,
                         ), states
    end
end

"""
    find_solve_sequence(sccs, vars)

given a set of `vars`, find the groups of equations we need to solve for
to obtain the solution to `vars`
"""
function find_solve_sequence(sccs, vars)
    subset = filter(i -> !isdisjoint(sccs[i], vars), 1:length(sccs))
    isempty(subset) && return []
    vars′ = mapreduce(i->sccs[i], union, subset)
    if vars′ == vars
        return subset
    else
        return find_solve_sequence(sccs, vars′)
    end
end

function build_observed_function(
        sys, ts, var_eq_matching, var_sccs,
        assignments,
        deps,
        sol_states,
        var2assignment;
        expression=false,
        output_type=Array,
        checkbounds=true,
    )

    is_not_prepended_assignment = trues(length(assignments))
    if (isscalar = !(ts isa AbstractVector))
        ts = [ts]
    end
    ts = Symbolics.scalarize.(value.(ts))

    vars = Set()
    foreach(Base.Fix1(vars!, vars), ts)
    ivs = independent_variables(sys)
    dep_vars = collect(setdiff(vars, ivs))

    s = structure(sys)
    @unpack fullvars, graph = s
    diffvars = map(i->fullvars[i], diffvars_range(s))
    algvars = map(i->fullvars[i], algvars_range(s))

    required_algvars = Set(intersect(algvars, vars))
    obs = observed(sys)
    observed_idx = Dict(map(x->x.lhs, obs) .=> 1:length(obs))
    # FIXME: this is a rather rough estimate of dependencies.
    maxidx = 0
    sts = Set(states(sys))
    for (i, s) in enumerate(dep_vars)
        idx = get(observed_idx, s, nothing)
        if idx === nothing
            if !(s in sts)
                throw(ArgumentError("$s is either an observed nor a state variable."))
            end
            continue
        end
        idx > maxidx && (maxidx = idx)
    end
    vs = Set()
    for idx in 1:maxidx
        vars!(vs, obs[idx].rhs)
        union!(required_algvars, intersect(algvars, vs))
        empty!(vs)
    end

    varidxs = findall(x->x in required_algvars, fullvars)
    subset = find_solve_sequence(var_sccs, varidxs)
    if !isempty(subset)
        eqs = equations(sys)

        torn_eqs  = map(i->map(v->eqs[var_eq_matching[v]], var_sccs[i]), subset)
        torn_vars = map(i->map(v->fullvars[v], var_sccs[i]), subset)
        u0map = defaults(sys)
        assignments = copy(assignments)
        solves = map(zip(torn_eqs, torn_vars)) do (eqs, vars)
            gen_nlsolve!(is_not_prepended_assignment, eqs, vars,
                         u0map, assignments, deps, var2assignment; checkbounds=checkbounds)
        end
    else
        solves = []
    end

    subs = []
    for sym in vars
        eqidx = get(observed_idx, sym, nothing)
        eqidx === nothing && continue
        push!(subs, sym ← obs[eqidx].rhs)
    end
    pre = get_postprocess_fbody(sys)

    ex = Code.toexpr(Func(
        [
         DestructuredArgs(diffvars, inbounds=!checkbounds)
         DestructuredArgs(parameters(sys), inbounds=!checkbounds)
         independent_variables(sys)
        ],
        [],
        pre(Let(
            [
             collect(Iterators.flatten(solves))
             map(eq -> eq.lhs←eq.rhs, obs[1:maxidx])
             subs
             assignments[is_not_prepended_assignment]
            ],
            isscalar ? ts[1] : MakeArray(ts, output_type)
           ))
       ), sol_states)

    expression ? ex : @RuntimeGeneratedFunction(ex)
end

struct ODAEProblem{iip}
end

ODAEProblem(args...; kw...) = ODAEProblem{true}(args...; kw...)

"""
    ODAEProblem{iip}(sys, u0map, tspan, parammap = DiffEqBase.NullParameters(); kw...)

This constructor acts similar to the one for [`ODEProblem`](@ref) with the following changes:
`ODESystem`s can sometimes be further reduced if `structural_simplify` has
already been applied to them. This is done this constructor.
In these cases, the constructor uses the knowledge of the strongly connected
components calculated during the process of simplification as the basis for
building pre-simplified nonlinear systems in the implicit solving.
In summary: these problems are structurally modified, but could be
more efficient and more stable. Note, the returned object is still of type [`ODEProblem`](@ref).
"""
function ODAEProblem{iip}(
                          sys,
                          u0map,
                          tspan,
                          parammap=DiffEqBase.NullParameters();
                          callback = nothing,
                          kwargs...
                         ) where {iip}
    fun, dvs = build_torn_function(sys; kwargs...)
    ps = parameters(sys)
    defs = defaults(sys)

    u0 = ModelingToolkit.varmap_to_vars(u0map, dvs; defaults=defs)
    p = ModelingToolkit.varmap_to_vars(parammap, ps; defaults=defs)

    has_difference = any(isdifferenceeq, equations(sys))
    if has_continuous_events(sys)
        event_cb = generate_rootfinding_callback(sys; kwargs...)
    else
        event_cb = nothing
    end
    difference_cb = has_difference ? generate_difference_cb(sys; kwargs...) : nothing
    cb = merge_cb(event_cb, difference_cb)
    cb = merge_cb(cb, callback)

    if cb === nothing
        ODEProblem{iip}(fun, u0, tspan, p; kwargs...)
    else
        ODEProblem{iip}(fun, u0, tspan, p; callback=cb, kwargs...)
    end
end
