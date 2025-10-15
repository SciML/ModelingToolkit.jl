const REPEATED_SIMPLIFICATION_MESSAGE = "Structural simplification cannot be applied to a completed system. Double simplification is not allowed."
struct RepeatedStructuralSimplificationError <: Exception end
function Base.showerror(io::IO, e::RepeatedStructuralSimplificationError)
    print(io, REPEATED_SIMPLIFICATION_MESSAGE)
end
""""""
function mtkcompile(
        sys::System; additional_passes = (), simplify = false, split = true,
        allow_symbolic = false, allow_parameter = true, conservative = false, fully_determined = true,
        inputs = SymbolicT[], outputs = SymbolicT[],
        disturbance_inputs = SymbolicT[],
        kwargs...)
    isscheduled(sys) && throw(RepeatedStructuralSimplificationError())
    inputs = unwrap_vars(inputs)
    outputs = unwrap_vars(outputs)
    disturbance_inputs = unwrap_vars(disturbance_inputs)
    newsys = __mtkcompile(sys; simplify,
        allow_symbolic, allow_parameter, conservative, fully_determined,
        inputs, outputs, disturbance_inputs, additional_passes,
        kwargs...)
    for pass in additional_passes
        newsys = pass(newsys)
    end
    @set! newsys.parent = complete(sys; split = false, flatten = false)
    newsys = complete(newsys; split)
    return newsys
end
function __mtkcompile(sys::AbstractSystem; simplify = false,
        inputs::Vector{SymbolicT} = SymbolicT[], outputs::Vector{SymbolicT} = SymbolicT[],
        disturbance_inputs::Vector{SymbolicT} = SymbolicT[],
        sort_eqs = true,
        kwargs...)
    if has_noise_eqs(sys) && get_noise_eqs(sys) !== nothing
        sys = noise_to_brownians(sys; names = :αₘₜₖ)
    end
    if !isempty(jumps(sys))
        return sys
    end
    if isempty(equations(sys)) && !is_time_dependent(sys) && !_iszero(cost(sys))
        return simplify_optimization_system(sys; kwargs..., sort_eqs, simplify)::System
    end
    sys, statemachines = extract_top_level_statemachines(sys)
    sys = expand_connections(sys)
    state = TearingState(sys)
    append!(state.statemachines, statemachines)
    @unpack structure, fullvars = state
    @unpack graph, var_to_diff, var_types = structure
    brown_vars = Int[]
    new_idxs = zeros(Int, length(var_types))
    idx = 0
    for (i, vt) in enumerate(var_types)
        if vt === BROWNIAN
            push!(brown_vars, i)
        else
            new_idxs[i] = (idx += 1)
        end
    end
    if isempty(brown_vars)
        return mtkcompile!(
            state; simplify, inputs, outputs, disturbance_inputs, kwargs...)
    else
        Is = Int[]
        Js = Int[]
        vals = Num[]
        make_eqs_zero_equals!(state)
        new_eqs = copy(equations(state))
        dvar2eq = Dict{Any, Int}()
        for (v, dv) in enumerate(var_to_diff)
            dv === nothing && continue
            deqs = 𝑑neighbors(graph, dv)
            if length(deqs) != 1
                error("$(eqs[deqs]) is not handled.")
            end
            dvar2eq[fullvars[dv]] = only(deqs)
        end
        for (j, bj) in enumerate(brown_vars), i in 𝑑neighbors(graph, bj)
            push!(Is, i)
            push!(Js, j)
            eq = new_eqs[i]
            brown = fullvars[bj]
            (coeff, residual, islinear) = Symbolics.linear_expansion(eq, brown)
            islinear || error("$brown isn't linear in $eq")
            new_eqs[i] = 0 ~ residual
            push!(vals, coeff)
        end
        g = Matrix(sparse(Is, Js, vals))
        sys = state.sys
        @set! sys.eqs = new_eqs
        @set! sys.unknowns = [v
                              for (i, v) in enumerate(fullvars)
                              if !iszero(new_idxs[i]) &&
                                 invview(var_to_diff)[i] === nothing]
        ode_sys = mtkcompile(
            sys; simplify, inputs, outputs, disturbance_inputs, kwargs...)
        eqs = equations(ode_sys)
        sorted_g_rows = zeros(Num, length(eqs), size(g, 2))
        for (i, eq) in enumerate(eqs)
            dvar = eq.lhs
            _iszero(dvar) && break
            g_row = get(dvar2eq, dvar, 0)
            iszero(g_row) && error("$dvar isn't handled.")
            g_row > size(g, 1) && continue
            @views copyto!(sorted_g_rows[i, :], g[g_row, :])
        end
        if sorted_g_rows isa AbstractMatrix && size(sorted_g_rows, 2) == 1
            noise_eqs = reshape(sorted_g_rows[:, 1], (:, 1))
            is_scalar_noise = true
        elseif __num_isdiag_noise(sorted_g_rows)
            noise_eqs = __get_num_diag_noise(sorted_g_rows)
            is_scalar_noise = false
        else
            noise_eqs = sorted_g_rows
            is_scalar_noise = false
        end
        noise_eqs = substitute_observed(ode_sys, noise_eqs)
        ssys = System(Vector{Equation}(full_equations(ode_sys)),
            get_iv(ode_sys), unknowns(ode_sys), parameters(ode_sys); noise_eqs,
            name = nameof(ode_sys), observed = observed(ode_sys), defaults = defaults(sys),
            assertions = assertions(sys),
            guesses = guesses(sys), initialization_eqs = initialization_equations(sys),
            continuous_events = continuous_events(sys),
            discrete_events = discrete_events(sys))
        @set! ssys.parameter_dependencies = get_parameter_dependencies(sys)
        return ssys
    end
end
function simplify_optimization_system(sys::System; split = true, kwargs...)
    sys = flatten(sys)
    cons = constraints(sys)
    econs = Equation[]
    icons = Inequality[]
    for e in cons
        if e isa Equation
            push!(econs, e)
        elseif e isa Inequality
            push!(icons, e)
        end
    end
    irreducible_subs = Dict{SymbolicT, SymbolicT}()
    dvs = SymbolicT[]
    for var in unknowns(sys)
        sh = SU.shape(var)::SU.ShapeVecT
        if isempty(sh)
            push!(dvs, var)
        else
            append!(dvs, vec(collect(var)::Array{SymbolicT})::Vector{SymbolicT})
        end
    end
    for i in eachindex(dvs)
        var = dvs[i]
        if hasbounds(var)
            irreducible_subs[var] = irrvar = setirreducible(var, true)::SymbolicT
            dvs[i] = irrvar
        end
    end
    subst = SU.Substituter{false}(irreducible_subs, SU.default_substitute_filter)
    for i in eachindex(econs)
        econs[i] = subst(econs[i])
    end
    nlsys = System(econs, dvs, parameters(sys); name = :___tmp_nlsystem)
    snlsys = mtkcompile(nlsys; kwargs..., fully_determined = false)::System
    obs = observed(snlsys)
    seqs = equations(snlsys)
    trueobs, _ = unhack_observed(obs, seqs)
    subs = Dict{SymbolicT, SymbolicT}()
    for eq in trueobs
        subs[eq.lhs] = eq.rhs
    end
    cons_simplified = Union{Equation, Inequality}[]
    for eq in seqs
        push!(cons_simplified, fixpoint_sub(eq, subs))
    end
    for eq in icons
        push!(cons_simplified, fixpoint_sub(eq, subs))
    end
    setdiff!(dvs, keys(subs))
    newsts = dvs
    @set! sys.constraints = cons_simplified
    newobs = copy(observed(sys))
    append!(newobs, obs)
    @set! sys.observed = newobs
    newcosts = copy(get_costs(sys))
    for i in eachindex(newcosts)
        newcosts[i] = fixpoint_sub(newcosts[i], subs)
    end
    @set! sys.costs = newcosts
    @set! sys.unknowns = newsts
    return sys
end
function __num_isdiag_noise(mat)
    for i in axes(mat, 1)
        nnz = 0
        for j in axes(mat, 2)
            if !isequal(mat[i, j], 0)
                nnz += 1
            end
        end
        if nnz > 1
            return (false)
        end
    end
    true
end
function __get_num_diag_noise(mat)
    map(axes(mat, 1)) do i
        for j in axes(mat, 2)
            mij = mat[i, j]
            if !isequal(mij, 0)
                return mij
            end
        end
        0
    end
end
""""""
function map_variables_to_equations(sys::AbstractSystem; rename_dummy_derivatives = true)
    if !has_tearing_state(sys)
        throw(ArgumentError("$(typeof(sys)) is not supported."))
    end
    ts = get_tearing_state(sys)
    if ts === nothing
        throw(ArgumentError("`map_variables_to_equations` requires a simplified system. Call `mtkcompile` on the system before calling this function."))
    end
    dummy_sub = Dict()
    if rename_dummy_derivatives && has_schedule(sys) && (sc = get_schedule(sys)) !== nothing
        dummy_sub = Dict(v => k for (k, v) in sc.dummy_sub if isequal(default_toterm(k), v))
    end
    mapping = Dict{Union{Num, BasicSymbolic}, Equation}()
    eqs = equations(sys)
    for eq in eqs
        isdifferential(eq.lhs) || continue
        var = arguments(eq.lhs)[1]
        var = get(dummy_sub, var, var)
        mapping[var] = eq
    end
    graph = ts.structure.graph
    algvars = BitSet(findall(
        Base.Fix1(StructuralTransformations.isalgvar, ts.structure), 1:ndsts(graph)))
    algeqs = BitSet(findall(1:nsrcs(graph)) do eq
        all(!Base.Fix1(isdervar, ts.structure), 𝑠neighbors(graph, eq))
    end)
    alge_var_eq_matching = complete(maximal_matching(graph, in(algeqs), in(algvars)))
    for (i, eq) in enumerate(alge_var_eq_matching)
        eq isa Unassigned && continue
        mapping[get(dummy_sub, ts.fullvars[i], ts.fullvars[i])] = eqs[eq]
    end
    for eq in observed(sys)
        mapping[get(dummy_sub, eq.lhs, eq.lhs)] = eq
    end
    return mapping
end
""""""
discrete_compile_pass(p) = false
