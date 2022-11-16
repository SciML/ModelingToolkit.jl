struct ClockInference
    ts::TearingState
    eq_domain::Vector{TimeDomain}
    var_domain::Vector{TimeDomain}
    inferred::BitSet
end

function ClockInference(ts::TearingState)
    @unpack fullvars, structure = ts
    @unpack graph = structure
    eq_domain = Vector{TimeDomain}(undef, nsrcs(graph))
    var_domain = Vector{TimeDomain}(undef, ndsts(graph))
    inferred = BitSet()
    for (i, v) in enumerate(fullvars)
        d = get_time_domain(v)
        if d isa Union{AbstractClock, Continuous}
            push!(inferred, i)
            dd = d
        else
            dd = Inferred()
        end
        var_domain[i] = dd
    end
    ClockInference(ts, eq_domain, var_domain, inferred)
end

function infer_clocks!(ci::ClockInference)
    @unpack ts, eq_domain, var_domain, inferred = ci
    @unpack fullvars = ts
    @unpack graph = ts.structure
    # TODO: add a graph type to do this lazily
    var_graph = SimpleGraph(ndsts(graph))
    for eq in 𝑠vertices(graph)
        vvs = 𝑠neighbors(graph, eq)
        if !isempty(vvs)
            fv, vs = Iterators.peel(vvs)
            for v in vs
                add_edge!(var_graph, fv, v)
            end
        end
    end
    cc = connected_components(var_graph)
    for c′ in cc
        c = BitSet(c′)
        idxs = intersect(c, inferred)
        isempty(idxs) && continue
        if !allequal(var_domain[i] for i in idxs)
            display(fullvars[c′])
            throw(ClockInferenceException("Clocks are not consistent in connected component $(fullvars[c′])"))
        end
        vd = var_domain[first(idxs)]
        for v in c′
            var_domain[v] = vd
        end
    end

    for v in 𝑑vertices(graph)
        vd = var_domain[v]
        eqs = 𝑑neighbors(graph, v)
        isempty(eqs) && continue
        #eq = first(eqs)
        for eq in eqs
            eq_domain[eq] = vd
        end
    end

    return ci
end

function resize_or_push!(v, val, idx)
    n = length(v)
    if idx > n
        for i in (n + 1):idx
            push!(v, Int[])
        end
        resize!(v, idx)
    end
    push!(v[idx], val)
end

function split_system(ci::ClockInference)
    @unpack ts, eq_domain, var_domain, inferred = ci
    @unpack fullvars = ts
    @unpack graph, var_to_diff = ts.structure
    continuous_id = Ref(0)
    clock_to_id = Dict{TimeDomain, Int}()
    id_to_clock = TimeDomain[]
    eq_to_cid = Vector{Int}(undef, nsrcs(graph))
    cid_to_eq = Vector{Int}[]
    var_to_cid = Vector{Int}(undef, ndsts(graph))
    cid_to_var = Vector{Int}[]
    cid_counter = Ref(0)
    for (i, d) in enumerate(eq_domain)
        cid = let cid_counter = cid_counter, id_to_clock = id_to_clock,
            continuous_id = continuous_id

            get!(clock_to_id, d) do
                cid = (cid_counter[] += 1)
                push!(id_to_clock, d)
                if d isa Continuous
                    continuous_id[] = cid
                end
                cid
            end
        end
        eq_to_cid[i] = cid
        resize_or_push!(cid_to_eq, i, cid)
    end
    continuous_id = continuous_id[]
    input_idxs = map(_ -> Int[], 1:cid_counter[])
    inputs = map(_ -> Any[], 1:cid_counter[])
    nvv = length(var_domain)
    for i in 1:nvv
        d = var_domain[i]
        cid = get(clock_to_id, d, 0)
        @assert cid!==0 "Internal error! Variable $(fullvars[i]) doesn't have a inferred time domain."
        var_to_cid[i] = cid
        v = fullvars[i]
        #TODO: remove Inferred*
        if istree(v) && (o = operation(v)) isa Operator &&
           input_timedomain(o) != output_timedomain(o)
            push!(input_idxs[cid], i)
            push!(inputs[cid], fullvars[i])
        end
        resize_or_push!(cid_to_var, i, cid)
    end

    eqs = equations(ts)
    tss = similar(cid_to_eq, TearingState)
    for (id, ieqs) in enumerate(cid_to_eq)
        vars = cid_to_var[id]
        ts_i = ts
        fadj = Vector{Int}[]
        eqs_i = Equation[]
        eq_to_diff = DiffGraph(length(ieqs))
        ne = 0
        for (j, eq_i) in enumerate(ieqs)
            vars = copy(graph.fadjlist[eq_i])
            ne += length(vars)
            push!(fadj, vars)
            push!(eqs_i, eqs[eq_i])
            eq_to_diff[j] = ts_i.structure.eq_to_diff[eq_i]
        end
        @set! ts_i.structure.graph = complete(BipartiteGraph(ne, fadj, ndsts(graph)))
        @set! ts_i.structure.only_discrete = id != continuous_id
        @set! ts_i.sys.eqs = eqs_i
        @set! ts_i.structure.eq_to_diff = eq_to_diff
        tss[id] = ts_i
    end
    return tss, inputs, continuous_id
end

function generate_discrete_affect(syss, inputs, continuous_id; checkbounds = true,
                                  eval_module = @__MODULE__, eval_expression = true)
    out = Sym{Any}(:out)
    appended_parameters = parameters(syss[continuous_id])
    param_to_idx = Dict{Any, Int}(reverse(en) for en in enumerate(appended_parameters))
    offset = length(appended_parameters)
    affect_funs = []
    svs = []
    for (i, (sys, input)) in enumerate(zip(syss, inputs))
        i == continuous_id && continue
        subs = get_substitutions(sys)
        assignments = map(s -> Assignment(s.lhs, s.rhs), subs.subs)
        let_body = SetArray(!checkbounds, out, rhss(equations(sys)))
        let_block = Let(assignments, let_body, false)
        needed_cont_to_disc_obs = map(v -> arguments(v)[1], input)
        # TODO: filter the needed ones
        needed_disc_to_cont_obs = map(v -> arguments(v)[1], inputs[continuous_id])
        append!(appended_parameters, input, states(sys))
        disc_to_cont_idxs = map(Base.Fix1(getindex, param_to_idx), inputs[continuous_id])
        cont_to_disc_obs = build_explicit_observed_function(syss[continuous_id],
                                                            needed_cont_to_disc_obs,
                                                            throw = false,
                                                            expression = true,
                                                            output_type = SVector)
        @set! sys.ps = appended_parameters
        disc_to_cont_obs = build_explicit_observed_function(sys, needed_disc_to_cont_obs,
                                                            throw = false,
                                                            expression = true,
                                                            output_type = SVector)
        ni = length(input)
        ns = length(states(sys))
        disc = Func([
                        out,
                        DestructuredArgs(states(sys)),
                        DestructuredArgs(appended_parameters),
                        get_iv(sys),
                    ], [],
                    let_block)
        cont_to_disc_idxs = (offset + 1):(offset += ni)
        input_offset = offset
        disc_range = (offset + 1):(offset += ns)
        save_tuple = Expr(:tuple)
        for i in 1:ns
            push!(save_tuple.args, :(p[$(input_offset + i)]))
        end
        affect! = :(function (integrator, saved_values)
                        @unpack u, p, t = integrator
                        c2d_obs = $cont_to_disc_obs
                        d2c_obs = $disc_to_cont_obs
                        c2d_view = view(p, $cont_to_disc_idxs)
                        d2c_view = view(p, $disc_to_cont_idxs)
                        disc_state = view(p, $disc_range)
                        disc = $disc
                        # Write continuous info to discrete
                        # Write discrete info to continuous
                        copyto!(c2d_view, c2d_obs(integrator.u, p, t))
                        copyto!(d2c_view, d2c_obs(disc_state, p, t))
                        push!(saved_values.t, t)
                        push!(saved_values.saveval, $save_tuple)
                        disc(disc_state, disc_state, p, t)
                    end)
        sv = SavedValues(Float64, NTuple{ns, Float64})
        push!(affect_funs, affect!)
        push!(svs, sv)
    end
    if eval_expression
        affects = map(affect_funs) do a
            @RuntimeGeneratedFunction(eval_module, toexpr(LiteralExpr(a)))
        end
    else
        affects = map(a -> toexpr(LiteralExpr(a)), affect_funs)
    end
    defaults = Dict{Any, Any}(v => 0.0 for v in Iterators.flatten(inputs))
    return affects, svs, appended_parameters, defaults
end
