@data ClockVertex begin
    Variable(Int)
    Equation(Int)
    InitEquation(Int)
    Clock(SciMLBase.AbstractClock)
end

struct ClockInference{S}
    """Tearing state."""
    ts::S
    """The time domain (discrete clock, continuous) of each equation."""
    eq_domain::Vector{TimeDomain}
    """The time domain (discrete clock, continuous) of each initialization equation."""
    init_eq_domain::Vector{TimeDomain}
    """The output time domain (discrete clock, continuous) of each variable."""
    var_domain::Vector{TimeDomain}
    inference_graph::HyperGraph{ClockVertex.Type}
    """The set of variables with concrete domains."""
    inferred::BitSet
end

function ClockInference(ts::TransformationState)
    @unpack structure = ts
    @unpack graph = structure
    eq_domain = TimeDomain[ContinuousClock() for _ in 1:nsrcs(graph)]
    init_eq_domain = TimeDomain[ContinuousClock()
                                for _ in 1:length(initialization_equations(ts.sys))]
    var_domain = TimeDomain[ContinuousClock() for _ in 1:ndsts(graph)]
    inferred = BitSet()
    for (i, v) in enumerate(get_fullvars(ts))
        d = get_time_domain(ts, v)
        if is_concrete_time_domain(d)
            push!(inferred, i)
            var_domain[i] = d
        end
    end
    inference_graph = HyperGraph{ClockVertex.Type}()
    for i in 1:nsrcs(graph)
        add_vertex!(inference_graph, ClockVertex.Equation(i))
    end
    for i in eachindex(initialization_equations(ts.sys))
        add_vertex!(inference_graph, ClockVertex.InitEquation(i))
    end
    for i in 1:ndsts(graph)
        varvert = ClockVertex.Variable(i)
        add_vertex!(inference_graph, varvert)
        v = ts.fullvars[i]
        d = get_time_domain(v)
        is_concrete_time_domain(d) || continue
        dvert = ClockVertex.Clock(d)
        add_vertex!(inference_graph, dvert)
        add_edge!(inference_graph, (varvert, dvert))
    end
    ClockInference(ts, eq_domain, init_eq_domain, var_domain, inference_graph, inferred)
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

"""
Update the equation-to-time domain mapping by inferring the time domain from the variables.
"""
function infer_clocks!(ci::ClockInference)
    @unpack ts, eq_domain, init_eq_domain, var_domain, inferred, inference_graph = ci
    @unpack var_to_diff, graph = ts.structure
    fullvars = get_fullvars(ts)
    isempty(inferred) && return ci

    var_to_idx = Dict(fullvars .=> eachindex(fullvars))

    # all shifted variables have the same clock as the unshifted variant
    for (i, v) in enumerate(fullvars)
        iscall(v) || continue
        operation(v) isa Shift || continue
        unshifted = only(arguments(v))
        add_edge!(inference_graph, (
            ClockVertex.Variable(i), ClockVertex.Variable(var_to_idx[unshifted])))
    end

    # preallocated buffers:
    # variables in each equation
    varsbuf = Set()
    # variables in each argument to an operator
    arg_varsbuf = Set()
    # hyperedge for each equation
    hyperedge = Set{ClockVertex.Type}()
    # hyperedge for each argument to an operator
    arg_hyperedge = Set{ClockVertex.Type}()
    # mapping from `i` in `InferredDiscrete(i)` to the vertices in that inferred partition
    relative_hyperedges = Dict{Int, Set{ClockVertex.Type}}()

    function infer_equation(ieq, eq, is_initialization_equation)
        empty!(varsbuf)
        empty!(hyperedge)
        # get variables in equation
        vars!(varsbuf, eq; op = Symbolics.Operator)
        # add the equation to the hyperedge
        eq_node = if is_initialization_equation
            ClockVertex.InitEquation(ieq)
        else
            ClockVertex.Equation(ieq)
        end
        push!(hyperedge, eq_node)
        for var in varsbuf
            idx = get(var_to_idx, var, nothing)
            # if this is just a single variable, add it to the hyperedge
            if idx isa Int
                push!(hyperedge, ClockVertex.Variable(idx))
                # we don't immediately `continue` here because this variable might be a
                # `Sample` or similar and we want the clock information from it if it is.
            end
            # now we only care about synchronous operators
            iscall(var) || continue
            op = operation(var)
            is_timevarying_operator(op) || continue

            # arguments and corresponding time domains
            args = arguments(var)
            tdomains = input_timedomain(op)
            if !(tdomains isa AbstractArray || tdomains isa Tuple)
                tdomains = [tdomains]
            end
            nargs = length(args)
            ndoms = length(tdomains)
            if nargs != ndoms
                throw(ArgumentError("""
                Operator $op applied to $nargs arguments $args but only returns $ndoms \
                domains $tdomains from `input_timedomain`.
                """))
            end

            # each relative clock mapping is only valid per operator application
            empty!(relative_hyperedges)
            for (arg, domain) in zip(args, tdomains)
                empty!(arg_varsbuf)
                empty!(arg_hyperedge)
                # get variables in argument
                vars!(arg_varsbuf, arg; op = Union{Differential, Shift})
                # get hyperedge for involved variables
                for v in arg_varsbuf
                    vidx = get(var_to_idx, v, nothing)
                    vidx === nothing && continue
                    push!(arg_hyperedge, ClockVertex.Variable(vidx))
                end

                Moshi.Match.@match domain begin
                    # If the time domain for this argument is a clock, then all variables in this edge have that clock.
                    x::SciMLBase.AbstractClock => begin
                        # add the clock to the edge
                        push!(arg_hyperedge, ClockVertex.Clock(x))
                        # add the edge to the graph
                        add_edge!(inference_graph, arg_hyperedge)
                    end
                    # We only know that this time domain is inferred. Treat it as a unique domain, all we know is that the
                    # involved variables have the same clock.
                    InferredClock.Inferred() => add_edge!(inference_graph, arg_hyperedge)
                    # All `InferredDiscrete` with the same `i` have the same clock (including output domain) so we don't
                    # add the edge, and instead add this to the `relative_hyperedges` mapping.
                    InferredClock.InferredDiscrete(i) => begin
                        relative_edge = get!(() -> Set{ClockVertex.Type}(), relative_hyperedges, i)
                        union!(relative_edge, arg_hyperedge)
                    end
                end
            end

            outdomain = output_timedomain(op)
            Moshi.Match.@match outdomain begin
                x::SciMLBase.AbstractClock => begin
                    push!(hyperedge, ClockVertex.Clock(x))
                end
                InferredClock.Inferred() => nothing
                InferredClock.InferredDiscrete(i) => begin
                    buffer = get(relative_hyperedges, i, nothing)
                    if buffer !== nothing
                        union!(hyperedge, buffer)
                        delete!(relative_hyperedges, i)
                    end
                end
            end

            for (_, relative_edge) in relative_hyperedges
                add_edge!(inference_graph, relative_edge)
            end
        end

        add_edge!(inference_graph, hyperedge)
    end
    for (ieq, eq) in enumerate(equations(ts))
        infer_equation(ieq, eq, false)
    end
    for (ieq, eq) in enumerate(initialization_equations(ts.sys))
        infer_equation(ieq, eq, true)
    end

    clock_partitions = connectionsets(inference_graph)
    for partition in clock_partitions
        clockidxs = findall(vert -> Moshi.Data.isa_variant(vert, ClockVertex.Clock), partition)
        if isempty(clockidxs)
            push!(partition, ClockVertex.Clock(ContinuousClock()))
            push!(clockidxs, length(partition))
        end
        if length(clockidxs) > 1
            vidxs = Int[vert.:1
                        for vert in partition
                        if Moshi.Data.isa_variant(vert, ClockVertex.Variable)]
            clks = [vert.:1 for vert in view(partition, clockidxs)]
            throw(ArgumentError("""
            Found clock partition with multiple associated clocks. Involved variables: \
            $(fullvars[vidxs]). Involved clocks: $(clks).
            """))
        end

        clock = partition[only(clockidxs)].:1
        for vert in partition
            Moshi.Match.@match vert begin
                ClockVertex.Variable(i) => (var_domain[i] = clock)
                ClockVertex.Equation(i) => (eq_domain[i] = clock)
                ClockVertex.InitEquation(i) => (init_eq_domain[i] = clock)
                ClockVertex.Clock(_) => nothing
            end
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
    iscall(v) || return false
    o = operation(v)
    o isa Operator || return false
    itd = input_timedomain(o)
    allequal(itd) || return true
    isempty(itd) && return true
    otd = output_timedomain(o)
    itd[1] == otd || return true
    return false
end

"""
For multi-clock systems, create a separate system for each clock in the system, along with associated equations. Return the updated tearing state, and the sets of clocked variables associated with each time domain.
"""
function split_system(ci::ClockInference{S}) where {S}
    @unpack ts, eq_domain, init_eq_domain, var_domain, inferred = ci
    fullvars = get_fullvars(ts)
    @unpack graph = ts.structure
    continuous_id = Ref(0)
    clock_to_id = Dict{TimeDomain, Int}()
    id_to_clock = TimeDomain[]
    eq_to_cid = Vector{Int}(undef, nsrcs(graph))
    cid_to_eq = Vector{Int}[]
    init_eq_to_cid = Vector{Int}(undef, length(initialization_equations(ts.sys)))
    cid_to_init_eq = Vector{Int}[]
    var_to_cid = Vector{Int}(undef, ndsts(graph))
    cid_to_var = Vector{Int}[]
    # cid_counter = number of clocks
    cid_counter = Ref(0)

    # populates clock_to_id and id_to_clock
    # checks if there is a continuous_id (for some reason? clock to id does this too)
    for (i, d) in enumerate(eq_domain)
        cid = let cid_counter = cid_counter, id_to_clock = id_to_clock,
            continuous_id = continuous_id

            # Fill the clock_to_id dict as you go,
            # ContinuousClock() => 1, ...
            get!(clock_to_id, d) do
                cid = (cid_counter[] += 1)
                push!(id_to_clock, d)
                if d == ContinuousClock()
                    continuous_id[] = cid
                end
                cid
            end
        end
        eq_to_cid[i] = cid
        resize_or_push!(cid_to_eq, i, cid)
    end
    # NOTE: This assumes that there is at least one equation for each clock
    for _ in 1:length(cid_to_eq)
        push!(cid_to_init_eq, Int[])
    end
    for (i, d) in enumerate(init_eq_domain)
        cid = clock_to_id[d]
        init_eq_to_cid[i] = cid
        push!(cid_to_init_eq[cid], i)
    end
    continuous_id = continuous_id[]
    # for each clock partition what are the input (indexes/vars)
    input_idxs = map(_ -> Int[], 1:cid_counter[])
    inputs = map(_ -> Any[], 1:cid_counter[])
    # var_domain corresponds to fullvars/all variables in the system
    nvv = length(var_domain)
    # put variables into the right clock partition
    # keep track of inputs to each partition
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

    # breaks the system up into a continuous and 0 or more discrete systems
    tss = similar(cid_to_eq, S)
    for (id, (ieqs, iieqs, ivars)) in enumerate(zip(cid_to_eq, cid_to_init_eq, cid_to_var))
        ts_i = system_subset(ts, ieqs, iieqs, ivars)
        if id != continuous_id
            ts_i = shift_discrete_system(ts_i)
            @set! ts_i.structure.only_discrete = true
        end
        tss[id] = ts_i
    end
    # put the continuous system at the back
    if continuous_id != 0
        tss[continuous_id], tss[end] = tss[end], tss[continuous_id]
        inputs[continuous_id], inputs[end] = inputs[end], inputs[continuous_id]
        id_to_clock[continuous_id],
        id_to_clock[end] = id_to_clock[end],
        id_to_clock[continuous_id]
        continuous_id = lastindex(tss)
    end
    return tss, inputs, continuous_id, id_to_clock
end
