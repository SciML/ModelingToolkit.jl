using LinearAlgebra

using ModelingToolkit: isdifferenceeq, has_continuous_events, generate_rootfinding_callback, generate_difference_cb, merge_cb

const MAX_INLINE_NLSOLVE_SIZE = 8

function torn_system_jacobian_sparsity(sys, var_eq_matching, var_sccs, nlsolve_scc_idxs, states_idxs)
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
        BitSet(var_rename[v′] for v′ in neighborhood(dig, v, Inf; dir=:in) if var_rename[v′] != 0)
    end

    for scc in var_sccs
        if length(scc) >= 2
            deps = fused_var_deps[scc[1]]
            for c in 2:length(scc)
                union!(deps, fused_var_deps[c])
            end
        end
    end

    nlsolve_eqs = BitSet(var_eq_matching[c]::Int for c in nlsolve_vars if var_eq_matching[c] !== unassigned)

    var2idx = Dict(v => i for (i, v) in enumerate(states_idxs))
    nlsolve_vars_set = BitSet(nlsolve_vars)

    I = Int[]; J = Int[]
    eqidx = 0
    for ieq in 𝑠vertices(graph)
        ieq in nlsolve_eqs && continue
        eqidx += 1
        for ivar in 𝑠neighbors(graph, ieq)
            isdervar(s, ivar) && continue
            if var_rename[ivar] != 0
                push!(I, eqidx)
                push!(J, var2idx[ivar])
            else
                for dvar in fused_var_deps[ivar]
                    isdervar(s, dvar) && continue
                    dvar in nlsolve_vars_set && continue
                    push!(I, eqidx)
                    push!(J, var2idx[dvar])
                end
            end
        end
    end
    sparse(I, J, true)
end

"""
    exprs = gen_nlsolve(eqs::Vector{Equation}, vars::Vector, u0map::Dict; checkbounds = true)

Generate `SymbolicUtils` expressions for a root-finding function based on `eqs`,
as well as a call to the root-finding solver.

`exprs` is a two element vector
```
exprs = [fname = f, numerical_nlsolve(fname, ...)]
```

# Arguments:
- `eqs`: Equations to find roots of.
- `vars`: ???
- `u0map`: A `Dict` which maps variables in `eqs` to values, e.g., `defaults(sys)` if `eqs = equations(sys)`.
- `checkbounds`: Apply bounds checking in the generated code.
"""
function gen_nlsolve(eqs, vars, u0map::AbstractDict; checkbounds=true)
    isempty(vars) && throw(ArgumentError("vars may not be empty"))
    length(eqs) == length(vars) || throw(ArgumentError("vars must be of the same length as the number of equations to find the roots of"))
    rhss = map(x->x.rhs, eqs)
    # We use `vars` instead of `graph` to capture parameters, too.
    allvars = unique(collect(Iterators.flatten(map(ModelingToolkit.vars, rhss))))
    params = setdiff(allvars, vars) # these are not the subject of the root finding

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
        isscalar ? rhss[1] : MakeArray(rhss, SVector)
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

    [
     fname ← @RuntimeGeneratedFunction(f)
     DestructuredArgs(vars, inbounds=!checkbounds) ← solver_call
    ]
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
    for eq in eqs
        isdiffeq(eq) && push!(rhss, eq.rhs)
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
    torn_expr = []
    defs = defaults(sys)
    nlsolve_scc_idxs = Int[]

    needs_extending = false
    for (i, scc) in enumerate(var_sccs)
        #torn_vars = [s.fullvars[var] for var in scc if var_eq_matching[var] !== unassigned]
        torn_vars_idxs = Int[var for var in scc if var_eq_matching[var] !== unassigned]
        torn_eqs = [eqs[var_eq_matching[var]] for var in torn_vars_idxs]
        isempty(torn_eqs) && continue
        if length(torn_eqs) <= max_inlining_size
            append!(torn_expr, gen_nlsolve(torn_eqs, s.fullvars[torn_vars_idxs], defs, checkbounds=checkbounds))
            push!(nlsolve_scc_idxs, i)
        else
            needs_extending = true
            append!(rhss, map(x->x.rhs, torn_eqs))
            append!(states_idxs, torn_vars_idxs)
            append!(mass_matrix_diag, zeros(length(torn_eqs)))
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
                 torn_expr,
                 funbody
                ))
            )
    )
    if expression
        expr, states
    else
        observedfun = let sys = sys, dict = Dict()
            function generated_observed(obsvar, u, p, t)
                obs = get!(dict, value(obsvar)) do
                    build_observed_function(sys, obsvar, var_eq_matching, var_sccs, checkbounds=checkbounds)
                end
                obs(u, p, t)
            end
        end

        ODEFunction{true}(
                          @RuntimeGeneratedFunction(expr),
                          sparsity = jacobian_sparsity ? torn_system_jacobian_sparsity(sys, var_eq_matching, var_sccs, nlsolve_scc_idxs, states_idxs) : nothing,
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
        sys, ts, var_eq_matching, var_sccs;
        expression=false,
        output_type=Array,
        checkbounds=true
    )

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
        solves = gen_nlsolve.(torn_eqs, torn_vars, (u0map,); checkbounds=checkbounds)
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

    ex = Func(
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
            ],
            isscalar ? ts[1] : MakeArray(ts, output_type)
           ))
    ) |> Code.toexpr

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
