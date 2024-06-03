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
        d = get_time_domain(v)
        if is_concrete_time_domain(d)
            push!(inferred, i)
            var_domain[i] = d
        end
    end
    ClockInference(ts, eq_domain, var_domain, inferred)
end

struct NotInferedTimeDomain end
function error_sample_time(eq)
    error("$eq\ncontains `SampleTime` but it is not an infered discrete equation.")
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
        if neweq isa NotInferedTimeDomain
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
    istree(ex) || return ex
    op = operation(ex)
    args = arguments(ex)
    if op == SampleTime
        dt === nothing && return NotInferedTimeDomain()
        return dt
    else
        new_args = similar(args)
        for (i, arg) in enumerate(args)
            ex_arg = substitute_sample_time(arg, dt)
            if ex_arg isa NotInferedTimeDomain
                return ex_arg
            end
            new_args[i] = ex_arg
        end
        similarterm(ex, op, new_args, symtype(ex); metadata = metadata(ex))
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
    istree(v) && (o = operation(v)) isa Operator &&
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
    return tss, inputs, continuous_id, id_to_clock
end

function generate_discrete_affect(
        osys::AbstractODESystem, syss, inputs, continuous_id, id_to_clock;
        checkbounds = true,
        eval_module = @__MODULE__, eval_expression = true)
    @static if VERSION < v"1.7"
        error("The `generate_discrete_affect` function requires at least Julia 1.7")
    end
    use_index_cache = has_index_cache(osys) && get_index_cache(osys) !== nothing
    out = Sym{Any}(:out)
    appended_parameters = full_parameters(syss[continuous_id])
    offset = length(appended_parameters)
    param_to_idx = if use_index_cache
        Dict{Any, ParameterIndex}(p => parameter_index(osys, p)
        for p in appended_parameters)
    else
        Dict{Any, Int}(reverse(en) for en in enumerate(appended_parameters))
    end
    affect_funs = []
    init_funs = []
    svs = []
    clocks = TimeDomain[]
    for (i, (sys, input)) in enumerate(zip(syss, inputs))
        i == continuous_id && continue
        push!(clocks, id_to_clock[i])
        subs = get_substitutions(sys)
        assignments = map(s -> Assignment(s.lhs, s.rhs), subs.subs)
        let_body = SetArray(!checkbounds, out, rhss(equations(sys)))
        let_block = Let(assignments, let_body, false)
        needed_cont_to_disc_obs = map(v -> arguments(v)[1], input)
        # TODO: filter the needed ones
        fullvars = Set{Any}(eq.lhs for eq in observed(sys))
        for s in unknowns(sys)
            push!(fullvars, s)
        end
        needed_disc_to_cont_obs = []
        if use_index_cache
            disc_to_cont_idxs = ParameterIndex[]
        else
            disc_to_cont_idxs = Int[]
        end
        for v in inputs[continuous_id]
            _v = arguments(v)[1]
            if _v in fullvars
                push!(needed_disc_to_cont_obs, _v)
                push!(disc_to_cont_idxs, param_to_idx[v])
                continue
            end

            # If the held quantity is calculated through observed
            # it will be shifted forward by 1
            _v = Shift(get_iv(sys), 1)(_v)
            if _v in fullvars
                push!(needed_disc_to_cont_obs, _v)
                push!(disc_to_cont_idxs, param_to_idx[v])
                continue
            end
        end
        append!(appended_parameters, input)
        cont_to_disc_obs = build_explicit_observed_function(
            use_index_cache ? osys : syss[continuous_id],
            needed_cont_to_disc_obs,
            throw = false,
            expression = true,
            output_type = SVector)
        disc_to_cont_obs = build_explicit_observed_function(sys, needed_disc_to_cont_obs,
            throw = false,
            expression = true,
            output_type = SVector,
            op = Shift,
            ps = reorder_parameters(osys, appended_parameters))
        ni = length(input)
        ns = length(unknowns(sys))
        disc = Func(
            [
                out,
                DestructuredArgs(unknowns(osys)),
                DestructuredArgs.(reorder_parameters(osys, full_parameters(osys)))...,
                get_iv(sys)
            ],
            [],
            let_block) |> toexpr
        if use_index_cache
            cont_to_disc_idxs = [parameter_index(osys, sym) for sym in input]
            disc_range = [parameter_index(osys, sym) for sym in unknowns(sys)]
        else
            cont_to_disc_idxs = (offset + 1):(offset += ni)
            input_offset = offset
            disc_range = (offset + 1):(offset += ns)
        end
        save_vec = Expr(:ref, :Float64)
        if use_index_cache
            for unk in unknowns(sys)
                idx = parameter_index(osys, unk)
                push!(save_vec.args, :($(parameter_values)(p, $idx)))
            end
        else
            for i in 1:ns
                push!(save_vec.args, :(p[$(input_offset + i)]))
            end
        end
        empty_disc = isempty(disc_range)
        disc_init = if use_index_cache
            :(function (u, p, t)
                c2d_obs = $cont_to_disc_obs
                d2c_obs = $disc_to_cont_obs
                result = c2d_obs(u, p..., t)
                for (val, i) in zip(result, $cont_to_disc_idxs)
                    $(_set_parameter_unchecked!)(p, val, i; update_dependent = false)
                end

                disc_state = Tuple($(parameter_values)(p, i) for i in $disc_range)
                result = d2c_obs(disc_state, p..., t)
                for (val, i) in zip(result, $disc_to_cont_idxs)
                    # prevent multiple updates to dependents
                    _set_parameter_unchecked!(p, val, i; update_dependent = false)
                end
                discretes, repack, _ = $(SciMLStructures.canonicalize)(
                    $(SciMLStructures.Discrete()), p)
                repack(discretes) # to force recalculation of dependents
            end)
        else
            :(function (u, p, t)
                c2d_obs = $cont_to_disc_obs
                d2c_obs = $disc_to_cont_obs
                c2d_view = view(p, $cont_to_disc_idxs)
                d2c_view = view(p, $disc_to_cont_idxs)
                disc_unknowns = view(p, $disc_range)
                copyto!(c2d_view, c2d_obs(u, p, t))
                copyto!(d2c_view, d2c_obs(disc_unknowns, p, t))
            end)
        end

        # @show disc_to_cont_idxs
        # @show cont_to_disc_idxs
        # @show disc_range
        affect! = :(function (integrator, saved_values)
            @unpack u, p, t = integrator
            c2d_obs = $cont_to_disc_obs
            d2c_obs = $disc_to_cont_obs
            $(
                if use_index_cache
                quote
                    disc_unknowns = [$(parameter_values)(p, i) for i in $disc_range]
                    cache = copy(disc_unknowns) # Cache needed for atomic state update
                end
            else
                quote
                    c2d_view = view(p, $cont_to_disc_idxs)
                    d2c_view = view(p, $disc_to_cont_idxs)
                    disc_unknowns = view(p, $disc_range)
                    cache = copy(disc_unknowns)
                end
            end
            )
            # TODO: find a way to do this without allocating
            disc = $disc

            # Write continuous into to discrete: handles `Sample`
            # Write discrete into to continuous
            # Update discrete unknowns

            # At a tick, c2d must come first
            # state update comes in the middle
            # d2c comes last
            # @show t
            # @show "incoming", p
            $(
                if use_index_cache
                quote
                    result = c2d_obs(integrator.u, p..., t)
                    for (val, i) in zip(result, $cont_to_disc_idxs)
                        $(_set_parameter_unchecked!)(p, val, i; update_dependent = false)
                    end
                    if !$empty_disc
                        # NOTE: the first and third arguments to `disc` MAY NOT be aliased
                        disc(cache, integrator.u, p..., t) # Cache needed for atomic state update
                    end
                    copyto!(disc_unknowns, cache)
                    for (val, i) in zip(cache, $disc_range)
                        $(_set_parameter_unchecked!)(p, val, i; update_dependent = false)
                    end
                    result = d2c_obs(disc_unknowns, p..., t)
                    for (val, i) in zip(result, $disc_to_cont_idxs)
                        $(_set_parameter_unchecked!)(p, val, i; update_dependent = false)
                    end
                end
            else
                quote
                    copyto!(c2d_view, c2d_obs(integrator.u, p, t))
                    $empty_disc || disc(cache, disc_unknowns, p, t) # Cache needed for atomic state update
                    copyto!(disc_unknowns, cache)
                    copyto!(d2c_view, d2c_obs(disc_unknowns, p, t))
                end
            end
            )

            push!(saved_values.t, t)
            push!(saved_values.saveval, $save_vec)

            # @show "after d2c", p
            $(
                if use_index_cache
                quote
                    discretes, repack, _ = $(SciMLStructures.canonicalize)(
                        $(SciMLStructures.Discrete()), p)
                    repack(discretes)
                end
            end
            )
        end)
        sv = SavedValues(Float64, Vector{Float64})
        push!(affect_funs, affect!)
        push!(init_funs, disc_init)
        push!(svs, sv)
    end
    if eval_expression
        affects = map(affect_funs) do a
            drop_expr(@RuntimeGeneratedFunction(eval_module, toexpr(LiteralExpr(a))))
        end
        inits = map(init_funs) do a
            drop_expr(@RuntimeGeneratedFunction(eval_module, toexpr(LiteralExpr(a))))
        end
    else
        affects = map(a -> toexpr(LiteralExpr(a)), affect_funs)
        inits = map(a -> toexpr(LiteralExpr(a)), init_funs)
    end
    defaults = Dict{Any, Any}(v => 0.0 for v in Iterators.flatten(inputs))
    return affects, inits, clocks, svs, appended_parameters, defaults
end
