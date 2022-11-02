struct ClockInference
    ts::TearingState
    eq_domain::Vector{TimeDomain}
    var_domain::Vector{TimeDomain}
    # maybe BitSet
    inferred::Vector{Int}
end

function ClockInference(ts::TearingState)
    @unpack fullvars, structure = ts
    @unpack graph = structure
    eq_domain = Vector{TimeDomain}(undef, nsrcs(graph))
    var_domain = Vector{TimeDomain}(undef, ndsts(graph))
    inferred = Int[]
    for (i, v) in enumerate(fullvars)
        d = get_time_domain(v)
        if d === nothing
            dd = Inferred()
        else
            push!(inferred, i)
            dd = d
        end
        var_domain[i] = dd
    end
    ClockInference(ts, eq_domain, var_domain, inferred)
end

function infer_clocks!(ci::ClockInference)
    @unpack ts, eq_domain, var_domain, inferred = ci
    @unpack fullvars = ts
    @unpack graph = ts.structure
    # TODO: add a graph time to do this lazily
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
    inferred = BitSet(inferred)
    for c′ in cc
        c = BitSet(c′)
        idxs = intersect(c, inferred)
        isempty(idxs) && continue
        if !allequal(var_domain[i] for i in idxs)
            throw(ClockInferenceException("Clocks are not consistent in connected component $(fullvars[c])"))
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
        eq = first(eqs)
        eq_domain[eq] = vd
    end

    return ci
end
