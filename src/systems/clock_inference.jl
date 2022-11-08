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
        # TODO: just mark current and sample variables as inputs
    end
    return tss, inputs

    #id_to_clock, cid_to_eq, cid_to_var
end
