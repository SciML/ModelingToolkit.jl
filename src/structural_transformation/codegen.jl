function torn_system_jacobian_sparsity(sys)
    s = structure(sys)
    @unpack fullvars, graph, partitions = s

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

    # `avars2dvars` maps a algebraic variable to its differential variable
    # dependencies.
    avars2dvars = Dict{Int,Set{Int}}()
    c = 0
    for partition in partitions
        @unpack e_residual, v_residual = partition
        # initialization
        for tvar in v_residual
            avars2dvars[tvar] = Set{Int}()
        end
        for teq in e_residual
            c += 1
            for var in 𝑠neighbors(graph, teq)
                # Skip the tearing variables in the current partition, because
                # we are computing them from all the other states.
                Graphs.insorted(var, v_residual) && continue
                deps = get(avars2dvars, var, nothing)
                if deps === nothing # differential variable
                    @assert !isalgvar(s, var)
                    for tvar in v_residual
                        push!(avars2dvars[tvar], var)
                    end
                else # tearing variable from previous partitions
                    @assert isalgvar(s, var)
                    for tvar in v_residual
                        union!(avars2dvars[tvar], avars2dvars[var])
                    end
                end
            end
        end
    end

    dvrange = diffvars_range(s)
    dvar2idx = Dict(v=>i for (i, v) in enumerate(dvrange))
    I = Int[]; J = Int[]
    eqidx = 0
    for ieq in 𝑠vertices(graph)
        isalgeq(s, ieq) && continue
        eqidx += 1
        for ivar in 𝑠neighbors(graph, ieq)
            if isdiffvar(s, ivar)
                push!(I, eqidx)
                push!(J, dvar2idx[ivar])
            elseif isalgvar(s, ivar)
                for dvar in avars2dvars[ivar]
                    push!(I, eqidx)
                    push!(J, dvar2idx[dvar])
                end
            end
        end
    end
    sparse(I, J, true)
end

"""
    partitions_dag(s::SystemStructure)

Return a DAG (sparse matrix) of partitions to use for parallelism.
"""
function partitions_dag(s::SystemStructure)
    @unpack partitions, graph = s

    # `partvars[i]` contains all the states that appear in `partitions[i]`
    partvars = map(partitions) do partition
        ipartvars = Set{Int}()
        for req in partition.e_residual
            union!(ipartvars, 𝑠neighbors(graph, req))
        end
        ipartvars
    end

    I, J = Int[], Int[]
    n = length(partitions)
    for (i, partition) in enumerate(partitions)
        for j in i+1:n
            # The only way for a later partition `j` to depend on an earlier
            # partition `i` is when `partvars[j]` contains one of tearing
            # variables of partition `i`.
            if !isdisjoint(partvars[j], partition.v_residual)
                # j depends on i
                push!(I, i)
                push!(J, j)
            end
        end
    end

    sparse(I, J, true, n, n)
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
    u0 = [map(var->get(u0map, var, 1e-3), vars)...]
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

function get_torn_eqs_vars(sys; checkbounds=true)
    s = structure(sys)
    partitions = s.partitions
    vars = s.fullvars
    eqs = equations(sys)

    torn_eqs  = map(idxs-> eqs[idxs], map(x->x.e_residual, partitions))
    torn_vars = map(idxs->vars[idxs], map(x->x.v_residual, partitions))
    u0map = defaults(sys)

    gen_nlsolve.(torn_eqs, torn_vars, (u0map,), checkbounds=checkbounds)
end

function build_torn_function(
        sys;
        expression=false,
        jacobian_sparsity=true,
        checkbounds=false,
        kw...
    )

    rhss = []
    for eq in equations(sys)
        isdiffeq(eq) && push!(rhss, eq.rhs)
    end

    out = Sym{Any}(gensym("out"))
    odefunbody = SetArray(
        !checkbounds,
        out,
        rhss
    )

    s = structure(sys)
    states = map(i->s.fullvars[i], diffvars_range(s))
    syms = map(Symbol, states)
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
                 collect(Iterators.flatten(get_torn_eqs_vars(sys, checkbounds=checkbounds))),
                 odefunbody
                ))
            )
    )
    if expression
        expr
    else
        observedfun = let sys = sys, dict = Dict()
            function generated_observed(obsvar, u, p, t)
                obs = get!(dict, value(obsvar)) do
                    build_observed_function(sys, obsvar, checkbounds=checkbounds)
                end
                obs(u, p, t)
            end
        end

        ODEFunction{true}(
                          @RuntimeGeneratedFunction(expr),
                          sparsity = torn_system_jacobian_sparsity(sys),
                          syms = syms,
                          observed = observedfun,
                         )
    end
end

"""
    find_solve_sequence(partitions, vars)

given a set of `vars`, find the groups of equations we need to solve for
to obtain the solution to `vars`
"""
function find_solve_sequence(partitions, vars)
    subset = filter(x -> !isdisjoint(x.v_residual, vars), partitions)
    isempty(subset) && return []
    vars′ = mapreduce(x->x.v_residual, union, subset)
    if vars′ == vars
        return subset
    else
        return find_solve_sequence(partitions, vars′)
    end
end

function build_observed_function(
        sys, ts;
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
    @unpack partitions, fullvars, graph = s
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
    subset = find_solve_sequence(partitions, varidxs)
    if !isempty(subset)
        eqs = equations(sys)

        torn_eqs  = map(idxs-> eqs[idxs.e_residual], subset)
        torn_vars = map(idxs->fullvars[idxs.v_residual], subset)
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
function ODAEProblem{iip}(
                          sys,
                          u0map,
                          tspan,
                          parammap=DiffEqBase.NullParameters();
                          kw...
                         ) where {iip}
    s = structure(sys)
    @unpack fullvars = s
    dvs = map(i->fullvars[i], diffvars_range(s))
    ps = parameters(sys)
    defs = defaults(sys)

    u0 = ModelingToolkit.varmap_to_vars(u0map, dvs; defaults=defs)
    p = ModelingToolkit.varmap_to_vars(parammap, ps; defaults=defs)

    ODEProblem{iip}(build_torn_function(sys; kw...), u0, tspan, p; kw...)
end
