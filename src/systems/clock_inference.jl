struct ClockInference{S}
    ts::S
    eq_domain::Vector{TimeDomain}
    var_domain::Vector{TimeDomain}
    inferred::BitSet
end

function ClockInference(ts::TransformationState)
    @unpack structure = ts
    @unpack graph = structure
    eq_domain = TimeDomain[Continuous() for _ in 1:nsrcs(graph)]
    var_domain = TimeDomain[Continuous() for _ in 1:ndsts(graph)]
    inferred = BitSet()
    for (i, v) in enumerate(get_fullvars(ts))
        d = get_time_domain(ts, v)
        if is_concrete_time_domain(d)
            push!(inferred, i)
            var_domain[i] = d
        end
    end
    ClockInference(ts, eq_domain, var_domain, inferred)
end

struct NotInferredTimeDomain end
function error_sample_time(eq)
    error("$eq\ncontains `SampleTime` but it is not an Inferred discrete equation.")
end
function substitute_sample_time(ci::ClockInference, ts::TearingState)
    @unpack eq_domain = ci
    eqs = copy(equations(ts))
    @assert length(eqs) == length(eq_domain)
    for i in eachindex(eqs)
        eq = eqs[i]
        domain = eq_domain[i]
        dt = sampletime(domain)
        neweq = substitute_sample_time(eq, dt)
        if neweq isa NotInferredTimeDomain
            error_sample_time(eq)
        end
        eqs[i] = neweq
    end
    @set! ts.sys.eqs = eqs
    @set! ci.ts = ts
end

function substitute_sample_time(eq::Equation, dt)
    substitute_sample_time(eq.lhs, dt) ~ substitute_sample_time(eq.rhs, dt)
end

function substitute_sample_time(ex, dt)
    iscall(ex) || return ex
    op = operation(ex)
    args = arguments(ex)
    if op == SampleTime
        dt === nothing && return NotInferredTimeDomain()
        return dt
    else
        new_args = similar(args)
        for (i, arg) in enumerate(args)
            ex_arg = substitute_sample_time(arg, dt)
            if ex_arg isa NotInferredTimeDomain
                return ex_arg
            end
            new_args[i] = ex_arg
        end
        maketerm(typeof(ex), op, new_args, metadata(ex))
    end
end

function infer_clocks!(ci::ClockInference)
    @unpack ts, eq_domain, var_domain, inferred = ci
    @unpack var_to_diff, graph = ts.structure
    fullvars = get_fullvars(ts)
    isempty(inferred) && return ci
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
    for v in vertices(var_to_diff)
        if (v′ = var_to_diff[v]) !== nothing
            add_edge!(var_graph, v, v′)
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
        for eq in eqs
            eq_domain[eq] = vd
        end
    end

    ci = substitute_sample_time(ci, ts)
    return ci
end

function resize_or_push!(v, val, idx)
    n = length(v)
    if idx > n
        for _ in (n + 1):idx
            push!(v, Int[])
        end
        resize!(v, idx)
    end
    push!(v[idx], val)
end

function is_time_domain_conversion(v)
    iscall(v) && (o = operation(v)) isa Operator &&
        input_timedomain(o) != output_timedomain(o)
end

function split_system(ci::ClockInference{S}) where {S}
    @unpack ts, eq_domain, var_domain, inferred = ci
    fullvars = get_fullvars(ts)
    @unpack graph = ts.structure
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
                if d == Continuous()
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
        if is_time_domain_conversion(v)
            push!(input_idxs[cid], i)
            push!(inputs[cid], fullvars[i])
        end
        resize_or_push!(cid_to_var, i, cid)
    end

    tss = similar(cid_to_eq, S)
    for (id, ieqs) in enumerate(cid_to_eq)
        ts_i = system_subset(ts, ieqs)
        if id != continuous_id
            ts_i = shift_discrete_system(ts_i)
            @set! ts_i.structure.only_discrete = true
        end
        tss[id] = ts_i
    end
    if continuous_id != 0
        tss[continuous_id], tss[end] = tss[end], tss[continuous_id]
        inputs[continuous_id], inputs[end] = inputs[end], inputs[continuous_id]
        id_to_clock[continuous_id], id_to_clock[end] = id_to_clock[end],
        id_to_clock[continuous_id]
        continuous_id = lastindex(tss)
    end
    return tss, inputs, continuous_id, id_to_clock
end
