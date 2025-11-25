"""
    $(TYPEDSIGNATURES)

Generate the initialization system for `sys`. The initialization system is a system of
nonlinear equations that solve for the full set of initial conditions of `sys` given
specified constraints.

The initialization system can be of two types: time-dependent and time-independent.
Time-dependent initialization systems solve for the initial values of unknowns as well as
the values of solvable parameters of the system. Time-independent initialization systems
only solve for solvable parameters of the system.

# Keyword arguments

- `time_dependent_init`: Whether to create an initialization system for a time-dependent
  system. A time-dependent initialization requires a time-dependent `sys`, but a time-
  independent initialization can be created regardless.
- `op`: The operating point of user-specified initial conditions of variables in `sys`.
- `initialization_eqs`: Additional initialization equations to use apart from those in
  `initialization_equations(sys)`.
- `guesses`: Additional guesses to use apart from those in `guesses(sys)`.
- `default_dd_guess`: Default guess for dummy derivative variables in time-dependent
  initialization.
- `algebraic_only`: If `false`, does not use initialization equations (provided via the
  keyword or part of the system) to construct initialization.
- `check_defguess`: Whether to error when a variable does not have a default or guess
  despite ModelingToolkitBase expecting it to.
- `name`: The name of the initialization system.

All other keyword arguments are forwarded to the [`System`](@ref) constructor.
"""
function generate_initializesystem(
        sys::AbstractSystem; time_dependent_init = is_time_dependent(sys), kwargs...)
    if time_dependent_init
        generate_initializesystem_timevarying(sys; kwargs...)
    else
        generate_initializesystem_timeindependent(sys; kwargs...)
    end
end

"""
$(TYPEDSIGNATURES)

Generate `System` of nonlinear equations which initializes a problem from specified initial conditions of a time-dependent `AbstractSystem`.
"""
function generate_initializesystem_timevarying(sys::AbstractSystem;
        op = SymmapT(),
        initialization_eqs = Equation[],
        guesses = SymmapT(),
        default_dd_guess = Bool(0),
        fast_path = false,
        algebraic_only = false,
        check_units = true, check_defguess = false,
        name = nameof(sys), kwargs...)
    eqs = equations(sys)
    trueobs, eqs = unhack_observed(observed(sys), eqs)
    # remove any observed equations that directly or indirectly contain
    # delayed unknowns
    isempty(trueobs) || filter_delay_equations_variables!(sys, trueobs)

    # Firstly, all variables and observables are initialization unknowns
    init_vars_set = AtomicArraySet{OrderedDict{SymbolicT, Nothing}}()
    add_trivial_initsys_vars!(init_vars_set, unknowns(sys), trueobs)

    eqs_ics = Equation[]

    inps = copy(get_inputs(sys))
    ps = parameters(sys; initial_parameters = true)
    init_ps = AtomicArraySet{OrderedDict{SymbolicT, Nothing}}()

    for v in get_all_discretes_fast(sys)
        push!(is_variable_floatingpoint(v) ? init_vars_set : init_ps, v)
    end
    for v in inps
        Moshi.Match.@match v begin
            BSImpl.Term(; f, args) => begin
                if f === getindex
                    push!(init_ps, args[1])
                elseif f isa SymbolicT
                    push!(init_ps, v)
                else
                    error("Unexpected input $v.")
                end
            end
            # Intentionally no fallback case. All inputs are originally variables.
        end
    end
    push!(init_ps, get_iv(sys)::SymbolicT)
    initsys_sort_system_parameters!(init_vars_set, init_ps, ps)

    guesses = as_atomic_dict_with_defaults(Dict{SymbolicT, SymbolicT}(guesses), COMMON_NOTHING)
    left_merge!(guesses, ModelingToolkitBase.guesses(sys))

    # Anything with a binding of `missing` is solvable.
    binds = bindings(sys)
    newbinds = SymmapT()
    # All bound parameters are solvable. The corresponding equation comes from the binding
    for v in bound_parameters(sys)
         push!(is_variable_floatingpoint(v) ? init_vars_set : init_ps, v)
    end
    initsys_sort_system_bindings!(init_vars_set, init_ps, eqs_ics, binds, newbinds, guesses)

    derivative_rules = DerivativeDict()
    dd_guess_sym = BSImpl.Const{VartypeT}(default_dd_guess)
    banned_derivatives = Set{SymbolicT}()
    if has_schedule(sys) && (schedule = get_schedule(sys); schedule isa Schedule)
        for (k, v) in schedule.dummy_sub
            if !has_possibly_indexed_key(guesses, k)
                write_possibly_indexed_array!(guesses, k, dd_guess_sym, COMMON_NOTHING)
            end
            ttk = default_toterm(k)
            # For DDEs, the derivatives can have delayed terms
            if _has_delays(sys, v, banned_derivatives)
                push!(banned_derivatives, ttk)
                continue
            end
            push_as_atomic_array!(init_vars_set, ttk)
            isequal(ttk, v) || push!(eqs_ics, ttk ~ v)
            derivative_rules[k] = ttk
        end
        merge!(derivative_rules, as_atomic_dict_with_defaults(Dict{SymbolicT, SymbolicT}(derivative_rules), COMMON_NOTHING))
    end
    for eq in eqs
        if _has_delays(sys, eq.rhs, banned_derivatives)
            isdiffeq(eq) && push!(banned_derivatives, default_toterm(eq.lhs))
            continue
        end
        if isdiffeq(eq)
            get!(derivative_rules, eq.lhs) do
                ttk = default_toterm(eq.lhs)
                push_as_atomic_array!(init_vars_set, ttk)
                push!(eqs_ics, ttk ~ eq.rhs)
                ttk
            end
        else
            push!(eqs_ics, eq)
        end
    end
    D = Differential(get_iv(sys))
    for eq in trueobs
        # Observed derivatives aren't added the same way as dummy_sub/diffeqs because
        # doing so would require all observed equations to be symbolically differentiable.
        get!(derivative_rules, D(eq.lhs), D(eq.rhs))
        # Add as guesses
        if !has_possibly_indexed_key(guesses, eq.lhs)
            write_possibly_indexed_array!(guesses, eq.lhs, eq.rhs, COMMON_NOTHING)
        end
    end
    op::SymmapT = if fast_path
        op
    else
        build_operating_point(sys, op)
    end
    timevaring_initsys_process_op!(init_vars_set, init_ps, eqs_ics, op, derivative_rules, guesses)

    # process explicitly provided initialization equations
    if !algebraic_only
        initialization_eqs = [get_initialization_eqs(sys); initialization_eqs]
        for eq in initialization_eqs
            eq = fixpoint_sub(eq, derivative_rules; maxiters = get_maxiters(derivative_rules)) # expand dummy derivatives
            push!(eqs_ics, eq)
        end
    end
    # TODO
    # 8) use observed equations for guesses of observed variables if not provided
    # guessed = Set(keys(defs)) # x(t), D(x(t)), ...
    # guessed = union(guessed, Set(default_toterm.(guessed))) # x(t), D(x(t)), xˍt(t), ...
    # for eq in trueobs
    #     if !(eq.lhs in guessed)
    #         defs[eq.lhs] = eq.rhs
    #         #push!(guessed, eq.lhs) # should not encounter eq.lhs twice, so don't need to track it
    #     end
    # end

    append!(eqs_ics, trueobs)

    vars = collect(init_vars_set)
    pars = collect(init_ps)
    System(Vector{Equation}(eqs_ics),
        vars,
        pars;
        bindings = newbinds,
        initial_conditions = guesses,
        checks = check_units,
        name,
        is_initializesystem = true,
        discover_from_metadata = false,
        kwargs...)
end

get_maxiters(subrules::AbstractDict) = max(3, min(1000, length(subrules)))

"""
$(TYPEDSIGNATURES)

Generate `System` of nonlinear equations which initializes a problem from specified initial conditions of a time-independent `AbstractSystem`.
"""
function generate_initializesystem_timeindependent(sys::AbstractSystem;
        op = Dict(),
        initialization_eqs = [],
        guesses = Dict(),
        algebraic_only = false,
        check_units = true, check_defguess = false,
        fast_path = false,
        name = nameof(sys), kwargs...)
    eqs = equations(sys)
    trueobs, eqs = unhack_observed(observed(sys), eqs)
    # remove any observed equations that directly or indirectly contain
    # delayed unknowns
    isempty(trueobs) || filter_delay_equations_variables!(sys, trueobs)

    og_dvs = as_atomic_array_set(unknowns(sys))
    union!(og_dvs, as_atomic_array_set(observables(sys)))

    init_vars_set = AtomicArraySet{OrderedDict{SymbolicT, Nothing}}()

    eqs_ics = Equation[]

    ps = parameters(sys; initial_parameters = true)
    init_ps = AtomicArraySet{OrderedDict{SymbolicT, Nothing}}()
    initsys_sort_system_parameters!(init_vars_set, init_ps, ps)

    guesses = SymmapT(guesses)
    left_merge!(guesses, ModelingToolkitBase.guesses(sys))

    # Anything with a binding of `missing` is solvable.
    binds = bindings(sys)
    newbinds = SymmapT()
    # All bound parameters are solvable. The corresponding equation comes from the binding
    for v in bound_parameters(sys)
         push!(is_variable_floatingpoint(v) ? init_vars_set : init_ps, v)
    end
    # Anything with a binding of `missing` is solvable.
    for (k, v) in binds
        if v === COMMON_MISSING
            push!(init_vars_set, k)
            delete!(init_ps, k)
            continue
        end
        k in og_dvs && continue
        if is_variable_floatingpoint(k)
            push!(eqs_ics, k ~ v)
            get!(guesses, k, v)
        else
            newbinds[k] = v
        end
    end

    op::SymmapT = if fast_path
        op
    else
        build_operating_point(sys, op)
    end

    valid_initial_parameters = AtomicArraySet{OrderedDict{SymbolicT, Nothing}}()
    for (k, v) in op
        if is_variable(sys, k) || has_observed_with_lhs(sys, k) ||
            Moshi.Match.@match k begin
                BSImpl.Term(; f, args) && if f isa Differential end => is_variable(sys, args[1])
                _ => false
            end

            isconst(v) && push!(valid_initial_parameters, Initial(k))
            continue
        end

        if v === COMMON_MISSING
            push!(init_vars_set, k)
            delete!(init_ps, k)
            continue
        end

        # No need to process any non-solvables
        if k in init_ps
            continue
        end

        if isconst(v)
            push!(eqs_ics, k ~ Initial(k))
            op[Initial(k)] = v
        else
            push!(eqs_ics, k ~ v)
        end
    end

    for k in valid_initial_parameters
        op[k] = Moshi.Data.variant_getfield(k, BSImpl.Term, :args)[1]
    end

    # process explicitly provided initialization equations
    if !algebraic_only
        initialization_eqs = [get_initialization_eqs(sys); initialization_eqs]
    end

    # only include initialization equations where all the involved `Initial`
    # parameters are valid.
    vs = Set{SymbolicT}()
    allpars = as_atomic_array_set(parameters(sys; initial_parameters = true))
    union!(allpars, bound_parameters(sys))
    initialization_eqs = filter(initialization_eqs) do eq
        empty!(vs)
        SU.search_variables!(vs, eq; is_atomic = OperatorIsAtomic{Initial}())
        # error if non-parameters are present in the initialization equations
        non_params = filter(!Base.Fix1(contains_possibly_indexed_element, allpars), vs)
        if !isempty(non_params)
            throw(UnknownsInTimeIndependentInitializationError(eq, non_params))
        end
        filter!(x -> iscall(x) && isinitial(x), vs)
        return issubset(vs, valid_initial_parameters)
        invalid_initials = setdiff(vs, valid_initial_parameters)
        return isempty(invalid_initials)
    end

    append!(eqs_ics, initialization_eqs)

    vars = collect(init_vars_set)
    pars = collect(init_ps)
    System(Vector{Equation}(eqs_ics),
        vars,
        pars;
        initial_conditions = guesses,
        checks = check_units,
        name,
        is_initializesystem = true,
        discover_from_metadata = false,
        kwargs...)
end

function add_trivial_initsys_vars!(init_vars_set::AtomicArraySet{OrderedDict{SymbolicT, Nothing}}, dvs::Vector{SymbolicT}, trueobs::Vector{Equation})
    for v in dvs
        push!(init_vars_set, split_indexed_var(v)[1])
    end
    for eq in trueobs
        push!(init_vars_set, split_indexed_var(eq.lhs)[1])
    end
end

function initsys_sort_system_parameters!(init_vars_set::AtomicArraySet{OrderedDict{SymbolicT, Nothing}},
                                         init_ps::AtomicArraySet{OrderedDict{SymbolicT, Nothing}},
                                         ps::Vector{SymbolicT})
    for v in ps
        arr, _ = split_indexed_var(v)
        arr in init_vars_set && continue
        push!(init_ps, arr)
    end
end

function initsys_sort_system_bindings!(init_vars_set::AtomicArraySet{OrderedDict{SymbolicT, Nothing}},
                                       init_ps::AtomicArraySet{OrderedDict{SymbolicT, Nothing}},
                                       eqs_ics::Vector{Equation}, binds::ROSymmapT,
                                       newbinds::SymmapT, guesses::SymmapT)
    # Anything with a binding of `missing` is solvable.
    for (k, v) in binds
        if v === COMMON_MISSING
            push!(init_vars_set, k)
            delete!(init_ps, k)
            continue
        end
        if is_variable_floatingpoint(k)
            push!(eqs_ics, k ~ v)
            get!(guesses, k, v)
        else
            newbinds[k] = v
        end
    end
end

function timevaring_initsys_process_op!(init_vars_set::AtomicArraySet{OrderedDict{SymbolicT, Nothing}},
                                        init_ps::AtomicArraySet{OrderedDict{SymbolicT, Nothing}},
                                        eqs_ics::Vector{Equation}, op::SymmapT,
                                        derivative_rules::DerivativeDict, guesses::SymmapT)
    for (k, v) in op
        # Late binding `missing` also makes the key solvable
        if v === COMMON_MISSING
            push!(init_vars_set, k)
            delete!(init_ps, k)
            continue
        end
        # No need to process any non-solvables
        if k in init_ps
            continue
        end

        # At this point, not only is `k` solvable but it should also have
        # `Initial(k)` defined if required.
        ik = Initial(k)
        @assert ik in init_ps
        subk = fixpoint_sub(k, derivative_rules; maxiters = get_maxiters(derivative_rules))
        # FIXME: DAEs can have initial conditions that require reducing the system
        # to index zero. If `isdifferential(y)`, an initial condition was given for the
        # derivative of an algebraic variable, so ignore it. Otherwise, the initialization
        # system gets a `D(y) ~ ...` equation and errors. This is the same behavior as v9.
        if isdifferential(subk)
            continue
        end
        shk = SU.shape(k)
        if SU.isconst(v)
            # The operating point already has `Initial(x) => x`. This same operating point
            # will be passed to the `NonlinearProblem` constructor, and guesses will not take
            # priority over it. So instead of adding `Initial(x) => v` as a guess, add `x => v`.
            op[ik] = k
            left_merge!(guesses, AtomicArrayDict(k => v))
            if !SU.is_array_shape(shk)
                push!(eqs_ics, subk ~ ik)
                continue
            end
            for i in SU.stable_eachindex(k)
                v[i] === COMMON_NOTHING && continue
                push!(eqs_ics, subk[i] ~ ik[i])
            end
            continue
        end
        if Symbolics.isarraysymbolic(k)
            for idx in SU.stable_eachindex(k)
                vv = v[idx]
                # `as_atomic_dict_with_defaults` is used to build `op`, which in
                # `build_operating_point` will put `COMMON_NOTHING` for missing
                # entries. Ignore them.
                vv === COMMON_NOTHING && continue
                kk = k[idx]
                subkk = subk[idx]
                ikk = Initial(kk)
                if SU.isconst(vv)
                    push!(eqs_ics, subkk ~ ikk)
                else
                    write_possibly_indexed_array!(guesses, kk, vv, COMMON_FALSE)
                    vv = fixpoint_sub(vv, derivative_rules; maxiters = get_maxiters(derivative_rules))
                    push!(eqs_ics, subkk ~ vv)
                end
            end
        else
            v = fixpoint_sub(v, derivative_rules; maxiters = get_maxiters(derivative_rules))
            isequal(subk, v) || push!(eqs_ics, subk ~ v)
        end
    end
end

"""
    $(TYPEDSIGNATURES)

Given `sys` and a list of observed equations `trueobs`, remove all the equations that
directly or indirectly contain a delayed unknown of `sys`.
"""
function filter_delay_equations_variables!(sys::AbstractSystem, trueobs::Vector{Equation})
    is_time_dependent(sys) || return trueobs
    banned_vars = Set{SymbolicT}()
    idxs_to_remove = Int[]
    for (i, eq) in enumerate(trueobs)
        _has_delays(sys, eq.rhs, banned_vars) || continue
        push!(idxs_to_remove, i)
        push!(banned_vars, eq.lhs)
    end
    return deleteat!(trueobs, idxs_to_remove)
end

"""
    $(TYPEDSIGNATURES)

Check if the expression `ex` contains a delayed unknown of `sys` or a term in
`banned`.
"""
function _has_delays(sys::AbstractSystem, ex, banned)
    ex = unwrap(ex)
    ex in banned && return true
    if symbolic_type(ex) == NotSymbolic()
        if is_array_of_symbolics(ex)
            return any(x -> _has_delays(sys, x, banned), ex)
        end
        return false
    end
    iscall(ex) || return false
    op = operation(ex)
    args = arguments(ex)
    if iscalledparameter(ex)
        return any(x -> _has_delays(sys, x, banned), args)
    end
    if issym(op) && length(args) == 1 && is_variable(sys, op(get_iv(sys))) &&
       iscall(args[1]) && get_iv(sys) in SU.search_variables(args[1])
        return true
    end
    return any(x -> _has_delays(sys, x, banned), args)
end

function SciMLBase.remake_initialization_data(
        sys::AbstractSystem, odefn, u0, t0, p, newu0, newp)
    if u0 === missing && p === missing
        return odefn.initialization_data
    end

    oldinitdata = odefn.initialization_data

    # We _always_ build initialization now. So if we didn't build  it before, don't do
    # it now
    oldinitdata === nothing && return nothing
    meta = oldinitdata.metadata
    meta isa InitializationMetadata || return oldinitdata

    if !(eltype(u0) <: Pair) && !(eltype(p) <: Pair)
        oldinitprob = oldinitdata.initializeprob
        oldinitprob === nothing && return nothing

        reconstruct_fn = meta.oop_reconstruct_u0_p
        # the history function doesn't matter because `reconstruct_fn` is only going to
        # update the values of parameters, which aren't time dependent. The reason it
        # is called is because `Initial` parameters are calculated from the corresponding
        # state values.
        history_fn = is_time_dependent(sys) && !is_markovian(sys) ? Returns(newu0) : nothing
        new_initu0,
        new_initp = reconstruct_fn(
            ProblemState(; u = newu0, p = newp, t = t0, h = history_fn), oldinitprob)
        if oldinitprob.f.resid_prototype === nothing
            newf = oldinitprob.f
        else
            newf = remake(oldinitprob.f;
                resid_prototype = calculate_resid_prototype(
                    length(oldinitprob.f.resid_prototype), new_initu0, new_initp))
        end
        initprob = remake(oldinitprob; f = newf, u0 = new_initu0, p = new_initp)
        return @set oldinitdata.initializeprob = initprob
    end

    dvs = unknowns(sys)
    ps = parameters(sys)
    if eltype(u0) <: Pair
        if u0 isa Array
            u0 = Dict(u0)
        end
        if keytype(u0) === Any || keytype(u0) <: Symbol
            u0 = anydict(u0)
            symbols_to_symbolics!(sys, u0)
        end
    else
        u0 = to_varmap(u0, dvs)
        symbols_to_symbolics!(sys, u0)
    end
    u0map = as_atomic_dict_with_defaults(Dict{SymbolicT, SymbolicT}(u0), COMMON_NOTHING)
    if eltype(p) <: Pair
        if p isa Array
            p = Dict(p)
        end
        if keytype(p) === Any || keytype(p) <: Symbol
            p = anydict(p)
            symbols_to_symbolics!(sys, p)
        end
    else
        p = to_varmap(p, ps)
        symbols_to_symbolics!(sys, p)
    end
    pmap = as_atomic_dict_with_defaults(Dict{SymbolicT, SymbolicT}(p), COMMON_NOTHING)
    op = merge!(u0map, pmap)
    guesses = SymmapT()
    use_scc = true
    initialization_eqs = Equation[]

    left_merge!(op, meta.op)
    filter!(Base.Fix2(!==, COMMON_NOTHING) ∘ last, op)
    merge!(guesses, meta.guesses)
    use_scc = meta.use_scc
    initialization_eqs = meta.additional_initialization_eqs
    time_dependent_init = meta.time_dependent_init

    if t0 === nothing && is_time_dependent(sys)
        t0 = 0.0
    end

    floatT = float_type_from_varmap(op)
    u0_constructor = get_u0_constructor(identity, typeof(newu0), floatT, false)
    p_constructor = get_p_constructor(identity, typeof(newu0), floatT)
    kws = maybe_build_initialization_problem(
        sys, SciMLBase.isinplace(odefn), op, t0, guesses;
        time_dependent_init, use_scc, initialization_eqs, floatT, fast_path = true,
        u0_constructor, p_constructor, allow_incomplete = true, check_units = false)

    odefn = remake(odefn; kws...)
    return SciMLBase.remake_initialization_data(sys, odefn, newu0, t0, newp, newu0, newp)
end

promote_type_with_nothing(::Type{T}, ::Nothing) where {T} = T
promote_type_with_nothing(::Type{T}, ::StaticVector{0}) where {T} = T
function promote_type_with_nothing(::Type{T}, ::AbstractArray{T2}) where {T, T2}
    promote_type(T, T2)
end
function promote_type_with_nothing(::Type{T}, p::MTKParameters) where {T}
    promote_type_with_nothing(promote_type_with_nothing(T, p.tunable), p.initials)
end

promote_with_nothing(::Type, ::Nothing) = nothing
promote_with_nothing(::Type, x::StaticVector{0}) = x
promote_with_nothing(::Type{T}, x::AbstractArray{T}) where {T} = x
function promote_with_nothing(::Type{T}, x::AbstractArray{T2}) where {T, T2}
    if ArrayInterface.ismutable(x)
        y = similar(x, T)
        copyto!(y, x)
        return y
    else
        yT = similar_type(x, T)
        return yT(x)
    end
end
function promote_with_nothing(::Type{T}, p::MTKParameters) where {T}
    tunables = promote_with_nothing(T, p.tunable)
    p = SciMLStructures.replace(SciMLStructures.Tunable(), p, tunables)
    initials = promote_with_nothing(T, p.initials)
    p = SciMLStructures.replace(SciMLStructures.Initials(), p, initials)
    return p
end

function promote_u0_p(u0, p, t0)
    T = Union{}
    T = promote_type_with_nothing(T, u0)
    T = promote_type_with_nothing(T, p)

    u0 = promote_with_nothing(T, u0)
    p = promote_with_nothing(T, p)
    return u0, p
end

function SciMLBase.late_binding_update_u0_p(
        prob, sys::AbstractSystem, u0, p, t0, newu0, newp)
    supports_initialization(sys) || return newu0, newp
    prob isa IntervalNonlinearProblem && return newu0, newp
    prob isa LinearProblem && return newu0, newp

    initdata = prob.f.initialization_data
    meta = initdata === nothing ? nothing : initdata.metadata

    newu0, newp = promote_u0_p(newu0, newp, t0)

    # non-symbolic u0 updates initials...
    if eltype(u0) <: Pair
        syms = []
        vals = []
        allsyms = all_symbols(sys)
        for (k, v) in u0
            v === nothing && continue
            (symbolic_type(v) == NotSymbolic() && !is_array_of_symbolics(v)) || continue
            if k isa Symbol
                k2 = symbol_to_symbolic(sys, k; allsyms)
                # if it is returned as-is, there is no match so skip it
                k2 === k && continue
                k = k2
            end
            is_parameter(sys, Initial(k)) || continue
            push!(syms, Initial(k))
            push!(vals, v)
        end
        newp = setp_oop(sys, syms)(newp, vals)
    else
        allsyms = nothing
        # if `p` is not provided or is symbolic
        p === missing || eltype(p) <: Pair || return newu0, newp
        (newu0 === nothing || isempty(newu0)) && return newu0, newp
        initdata === nothing && return newu0, newp
        meta = initdata.metadata
        meta isa InitializationMetadata || return newu0, newp
        newp = p === missing ? copy(newp) : newp

        if length(newu0) != length(prob.u0)
            throw(ArgumentError("Expected `newu0` to be of same length as unknowns ($(length(prob.u0))). Got $(typeof(newu0)) of length $(length(newu0))"))
        end
        newp = meta.set_initial_unknowns!(newp, newu0)
    end

    if eltype(p) <: Pair
        syms = []
        vals = []
        if allsyms === nothing
            allsyms = all_symbols(sys)
        end
        for (k, v) in p
            v === nothing && continue
            (symbolic_type(v) == NotSymbolic() && !is_array_of_symbolics(v)) || continue
            if k isa Symbol
                k2 = symbol_to_symbolic(sys, k; allsyms)
                # if it is returned as-is, there is no match so skip it
                k2 === k && continue
                k = k2
            end
            is_parameter(sys, Initial(k)) || continue
            push!(syms, Initial(k))
            push!(vals, v)
        end
        newp = setp_oop(sys, syms)(newp, vals)
    end

    return newu0, newp
end

function DiffEqBase.get_updated_symbolic_problem(
        sys::AbstractSystem, prob; u0 = state_values(prob),
        p = parameter_values(prob), kw...)
    supports_initialization(sys) || return prob
    initdata = prob.f.initialization_data
    initdata isa SciMLBase.OverrideInitData || return prob
    meta = initdata.metadata
    meta isa InitializationMetadata || return prob
    meta.get_updated_u0 === nothing && return prob

    u0 === nothing && return remake(prob; p)

    t0 = is_time_dependent(prob) ? current_time(prob) : nothing

    if p isa MTKParameters
        buffer = p.initials
    else
        buffer = p
    end

    u0 = DiffEqBase.promote_u0(u0, buffer, t0)

    if ArrayInterface.ismutable(u0)
        T = typeof(u0)
    else
        T = StaticArrays.similar_type(u0)
    end

    return remake(prob; u0 = T(meta.get_updated_u0(prob, initdata.initializeprob)), p)
end

"""
    $(TYPEDSIGNATURES)

Check if the given system is an initialization system.
"""
function is_initializesystem(sys::AbstractSystem)
    has_is_initializesystem(sys) && get_is_initializesystem(sys)
end

"""
Counteracts the CSE/array variable hacks in `symbolics_tearing.jl` so it works with
initialization.
"""
function unhack_observed(obseqs, eqs)
    return obseqs, eqs
end

function UnknownsInTimeIndependentInitializationError(eq, non_params)
    ArgumentError("""
    Initialization equations for time-independent systems can only contain parameters. \
    Found $non_params in $eq. If the equations refer to the initial guess for unknowns, \
    use the `Initial` operator.
    """)
end
