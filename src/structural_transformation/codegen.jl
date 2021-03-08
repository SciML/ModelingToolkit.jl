import Symbolics: ClosureExpr

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
    for (_, _, teqs, tvars) in partitions
        # initialization
        for tvar in tvars
            avars2dvars[tvar] = Set{Int}()
        end
        for teq in teqs
            c += 1
            for var in 𝑠neighbors(graph, teq)
                # Skip the tearing variables in the current partition, because
                # we are computing them from all the other states.
                LightGraphs.insorted(var, tvars) && continue
                deps = get(avars2dvars, var, nothing)
                if deps === nothing # differential variable
                    @assert !isalgvar(s, var)
                    for tvar in tvars
                        push!(avars2dvars[tvar], var)
                    end
                else # tearing variable from previous partitions
                    @assert isalgvar(s, var)
                    for tvar in tvars
                        union!(avars2dvars[tvar], avars2dvars[var])
                    end
                end
            end
        end
    end

    dvar2idx(idx) = idx # maps `dvar` to the index of the states
    I = Int[]; J = Int[]
    eqidx = 0
    for ieq in 𝑠vertices(graph)
        isalgeq(s, ieq) && continue
        eqidx += 1
        for ivar in 𝑠neighbors(graph, ieq)
            if isdiffvar(s, ivar)
                push!(I, eqidx)
                push!(J, dvar2idx(ivar))
            elseif isalgvar(s, ivar)
                for dvar in avars2dvars[ivar]
                    push!(I, eqidx)
                    push!(J, dvar2idx(dvar))
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
    partvars = map(partitions) do (_, _, reqs, tvars)
        ipartvars = Set{Int}()
        for req in reqs
            union!(ipartvars, 𝑠neighbors(graph, req))
        end
        ipartvars
    end

    I, J = Int[], Int[]
    n = length(partitions)
    for i in 1:n
        for j in i+1:n
            # The only way for a later partition `j` to depend on an earlier
            # partition `i` is when `partvars[j]` contains one of tearing
            # variables of partition `i`.
            if !isdisjoint(partvars[j], partitions[i][4])
                # j depends on i
                push!(I, i)
                push!(J, j)
            end
        end
    end

    sparse(I, J, true, n, n)
end

function gen_nlsolve(sys, eqs, vars, destructure=true)
    @assert !isempty(vars)
    @assert length(eqs) == length(vars)
    rhss = map(x->x.rhs, eqs)
    # We use `vars` instead of `graph` to capture parameters, too.
    allvars = unique(collect(Iterators.flatten(map(ModelingToolkit.vars, rhss))))
    params = setdiff(allvars, vars)

    u0map = default_u0(sys)
    # splatting to tighten the type
    u0 = [map(var->get(u0map, var, 1e-3), vars)...]
    # specialize on the scalar case
    isscalar = length(u0) == 1
    u0 = isscalar ? u0[1] : SVector(u0...)

    fname = gensym("fun")
    statearg = DestructuredArgs(vars)
    f = Func(
        [statearg,
         DestructuredArgs(params)
        ],
        [],
        isscalar ? rhss[1] : MakeArray(rhss, SVector)
    ) |> SymbolicUtils.Code.toexpr

    solver_call = LiteralExpr(quote
            $numerical_nlsolve(
                               $fname,
                               # initial guess
                               $u0,
                               # "captured variables"
                               ($(params...),)
                              )
        end)

    Dict(:solve_params => params,
         :solve_vars => vars,
         :u0 => u0,
         :function_def => fname ← @RuntimeGeneratedFunction(f),
         :output_type => LiteralExpr(:(typeof($u0))),
         :solver_call => solver_call)
end

function get_torn_eqs_vars(sys, parallelize)
    s = structure(sys)
    partitions = s.partitions
    vars = s.fullvars
    eqs = equations(sys)

    torn_eqs  = map(idxs-> eqs[idxs], map(x->x[3], partitions))
    torn_vars = map(idxs->vars[idxs], map(x->x[4], partitions))

    if parallelize
        dag = partitions_dag(s)
        spawned = Set{Int}()
        assigns = []
        while length(spawned) < size(dag, 2)
            next = findall(1:size(dag, 2)) do i
                !(i in spawned) && iszero(dag[:, i])
            end
            isempty(next) && error("Cannot parallelize")

            union!(spawned, next)

            dag[next, :] .= false

            nlsolves = gen_nlsolve.((sys,),
                                 torn_eqs[next],
                                 torn_vars[next],
                                 false)

            @show length(nlsolves)

            nt = Base.Threads.nthreads()
            batches = map(Iterators.partition(nlsolves, nt)) do batch
                params = reduce(union,
                                map(s->s[:solve_params],
                                    batch), init=[])
                solvercalls = map(s->s[:solver_call], batch)

                batchexpr = Let(map(s->s[:function_def],
                                    batch),
                    LiteralExpr(:(tuple($(solvercalls...)))))


                ot = LiteralExpr(
                    :(Tuple{$(map(s->s[:output_type],
                                  batch)...)}))
                ClosureExpr(
                    Func(params, [], batchexpr), params, ot)
            end

            varbatches = collect(Iterators.partition(map(DestructuredArgs, torn_vars[next]), nt))
            assign = DestructuredArgs(map(DestructuredArgs,
                                          varbatches)) ←
                SpawnFetch{MultithreadedForm}(batches, tuple)


            push!(assigns, assign)
        end
        return assigns
    end

    map(gen_nlsolve.((sys,), torn_eqs, torn_vars)) do s
        [s[:function_def],
         DestructuredArgs(s[:solve_vars]) ← s[:solver_call]]
    end |> Iterators.flatten |> collect
end

function build_torn_function(
        sys;
        expression=false,
        jacobian_sparsity=true,
        checkbounds=false,
        threaded=false,
        kw...
    )

    rhss = []
    for eq in equations(sys)
        isdiffeq(eq) && push!(rhss, eq.rhs)
    end

    s = structure(sys)

    out = Sym{Any}(gensym("out"))


    #=
    ps = vcat(parameters(sys),
              s.fullvars[diffvars_range(s)],
              s.fullvars[algvars_range(s)])

    Js = Symbolics.jacobian_sparsity(rhss, ps)
    idxs = findall(x->x>0, vec(sum(Js, dims=1)))

    closed_args = vcat(ps[idxs], [independent_variable(sys)])
    @show closed_args
    =#
    ps = parameters(sys)

    Js = Symbolics.jacobian_sparsity(rhss, ps)
    idxs = findall(x->x>0, vec(sum(Js, dims=1)))

    closed_args = vcat(s.fullvars[diffvars_range(s)],
                       s.fullvars[algvars_range(s)],
                       ps[idxs], [independent_variable(sys)])

    if threaded
        odefunbody = Symbolics.set_array(
            MultithreadedForm(Threads.nthreads()),
            closed_args,
            out,
            nothing, # outputidxs
            rhss,
            false, # checkbounds
            true   # skipzeros
        )
    else
        odefunbody = SetArray(
            !checkbounds,
            out,
            map(Symbolics.unflatten_long_ops, rhss)
        )
    end

    states = s.fullvars[diffvars_range(s)]
    syms = map(Symbol, states)

    expr = SymbolicUtils.Code.toexpr(
        Func(
             [
              out
              DestructuredArgs(states)
              DestructuredArgs(ps)
              independent_variable(sys)
             ],
             [],
             Let(
                 get_torn_eqs_vars(sys, threaded),
                 odefunbody
                )
            )
    )
    if isdefined(Main, :_debug_expr)
        Main._debug_expr[] = expr
    end
    if expression
        expr
    else
        observedfun = let sys = sys, dict = Dict()
            function generated_observed(obsvar, u, p, t)
                obs = get!(dict, value(obsvar)) do
                    build_observed_function(sys, obsvar)
                end
                obs(u, p, t)
            end
        end

        ODEFunction{true}(
                          @RuntimeGeneratedFunction(expr),
                          sparsity = jacobian_sparsity ? torn_system_jacobian_sparsity(sys) : nothing,
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
    subset = filter(x -> !isdisjoint(x[4], vars), partitions)
    isempty(subset) && return []
    vars′ = mapreduce(x->x[4], union, subset)
    if vars′ == vars
        return subset
    else
        return find_solve_sequence(partitions, vars′)
    end
end

function build_observed_function(
        sys, syms;
        expression=false,
        output_type=Array
    )

    if (isscalar = !(syms isa Vector))
        syms = [syms]
    end
    syms = value.(syms)
    syms_set = Set(syms)
    s = structure(sys)
    @unpack partitions, fullvars, graph = s
    diffvars = fullvars[diffvars_range(s)]
    algvars = fullvars[algvars_range(s)]

    required_algvars = Set(intersect(algvars, syms_set))
    obs = observed(sys)
    observed_idx = Dict(map(x->x.lhs, obs) .=> 1:length(obs))
    # FIXME: this is a rather rough estimate of dependencies.
    maxidx = 0
    for (i, s) in enumerate(syms)
        idx = get(observed_idx, s, nothing)
        idx === nothing && continue
        idx > maxidx && (maxidx = idx)
    end
    for idx in 1:maxidx
        vs = vars(obs[idx].rhs)
        union!(required_algvars, intersect(algvars, vs))
    end

    varidxs = findall(x->x in required_algvars, fullvars)
    subset = find_solve_sequence(partitions, varidxs)
    if !isempty(subset)
        eqs = equations(sys)

        torn_eqs  = map(idxs-> eqs[idxs[3]], subset)
        torn_vars = map(idxs->fullvars[idxs[4]], subset)

        solves = gen_nlsolve.((sys,), torn_eqs, torn_vars)
    else
        solves = []
    end

    solves = map(solves) do s
        [s[:function_def],
         DestructuredArgs(s[:solve_vars]) ← s[:solver_call]]
    end

    output = map(syms) do sym
        if sym in required_algvars
            sym
        else
            obs[observed_idx[sym]].rhs
        end
    end

    ex = Func(
        [
         DestructuredArgs(diffvars)
         DestructuredArgs(parameters(sys))
         independent_variable(sys)
        ],
        [],
        Let(
            [
             collect(Iterators.flatten(solves))
             map(eq -> eq.lhs←eq.rhs, obs[1:maxidx])
            ],
            isscalar ? output[1] : MakeArray(output, output_type)
           )
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
    dvs = fullvars[diffvars_range(s)]
    defaults = merge(default_p(sys), default_u0(sys))
    u0map′ = ModelingToolkit.lower_mapnames(u0map, independent_variable(sys))
    u0 = ModelingToolkit.varmap_to_vars(u0map′, dvs; defaults=defaults)

    ps = parameters(sys)
    if parammap isa DiffEqBase.NullParameters && isempty(default_p(sys))
        isempty(ps) || throw(ArgumentError("The model has non-empty parameters but no parameters are specified in the problem."))
        p = parammap
    else
        if parammap isa DiffEqBase.NullParameters
            pp = Pair[]
        else
            pp = ModelingToolkit.lower_mapnames(parammap)
        end
        p = ModelingToolkit.varmap_to_vars(pp, ps; defaults=defaults)
    end

    ODEProblem{iip}(build_torn_function(sys; kw...), u0, tspan, p; kw...)
end
