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
  despite ModelingToolkit expecting it to.
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
        op = Dict(),
        initialization_eqs = [],
        guesses = Dict(),
        default_dd_guess = Bool(0),
        algebraic_only = false,
        check_units = true, check_defguess = false,
        name = nameof(sys), kwargs...)
    eqs = equations(sys)
    if !(eqs isa Vector{Equation})
        eqs = Equation[x for x in eqs if x isa Equation]
    end
    trueobs, eqs = unhack_observed(observed(sys), eqs)
    # remove any observed equations that directly or indirectly contain
    # delayed unknowns
    isempty(trueobs) || filter_delay_equations_variables!(sys, trueobs)
    vars = unique([unknowns(sys); getfield.(trueobs, :lhs)])
    vars_set = Set(vars) # for efficient in-lookup
    arrvars = Set()
    for var in vars
        if iscall(var) && operation(var) === getindex
            push!(arrvars, first(arguments(var)))
        end
    end

    eqs_ics = Equation[]
    defs = copy(defaults(sys)) # copy so we don't modify sys.defaults
    additional_guesses = anydict(guesses)
    guesses = merge(get_guesses(sys), additional_guesses)
    idxs_diff = isdiffeq.(eqs)

    # PREPROCESSING
    op = anydict(op)
    if isempty(op)
        op = copy(defs)
    end
    scalarize_vars_in_varmap!(op, arrvars)
    u0map = anydict()
    pmap = anydict()
    build_operating_point!(sys, op, u0map, pmap, Dict(), unknowns(sys),
        parameters(sys; initial_parameters = true))
    for (k, v) in op
        if has_parameter_dependency_with_lhs(sys, k) && is_variable_floatingpoint(k)
            pmap[k] = v
        end
    end
    initsys_preprocessing!(u0map, defs)

    # 1) Use algebraic equations of system as initialization constraints
    idxs_alge = .!idxs_diff
    append!(eqs_ics, eqs[idxs_alge]) # start equation list with algebraic equations

    eqs_diff = eqs[idxs_diff]
    D = Differential(get_iv(sys))
    diffmap = merge(
        Dict(eq.lhs => eq.rhs for eq in eqs_diff),
        Dict(D(eq.lhs) => D(eq.rhs) for eq in trueobs)
    )

    if has_schedule(sys) && (schedule = get_schedule(sys); !isnothing(schedule))
        # 2) process dummy derivatives and u0map into initialization system
        # prepare map for dummy derivative substitution
        for x in filter(x -> !isnothing(x[1]), schedule.dummy_sub)
            # set dummy derivatives to default_dd_guess unless specified
            push!(defs, x[1] => get(guesses, x[1], default_dd_guess))
        end
        function process_u0map_with_dummysubs(y, x)
            y = get(schedule.dummy_sub, y, y)
            y = fixpoint_sub(y, diffmap)
            # FIXME: DAEs provide initial conditions that require reducing the system
            # to index zero. If `isdifferential(y)`, an initial condition was given for an
            # algebraic variable, so ignore it. Otherwise, the initialization system
            # gets a `D(y) ~ ...` equation and errors. This is the same behavior as v9.
            if isdifferential(y)
                return
            end
            # If we have `D(x) ~ x` and provide [D(x) => x, x => 1.0] to `u0map`, then
            # without this condition `defs` would get `x => x` instead of retaining
            # `x => 1.0`.
            isequal(y, x) && return
            if y ∈ vars_set
                # variables specified in u0 overrides defaults
                push!(defs, y => x)
            elseif y isa Symbolics.Arr
                # TODO: don't scalarize arrays
                merge!(defs, Dict(scalarize(y .=> x)))
            elseif y isa Symbolics.BasicSymbolic
                # y is a derivative expression expanded; add it to the initialization equations
                push!(eqs_ics, y ~ x)
            else
                error("Initialization expression $y is currently not supported. If its a higher order derivative expression, then only the dummy derivative expressions are supported.")
            end
        end
        for (y, x) in u0map
            if Symbolics.isarraysymbolic(y)
                process_u0map_with_dummysubs.(collect(y), collect(x))
            else
                process_u0map_with_dummysubs(y, x)
            end
        end
    else
        # TODO: Check if this is still necessary
        # 2) System doesn't have a schedule, so dummy derivatives don't exist/aren't handled (SDESystem)
        for (k, v) in u0map
            defs[k] = v
        end
    end

    # 3) process other variables
    for var in vars
        if var ∈ keys(op)
            push!(eqs_ics, var ~ op[var])
        elseif var ∈ keys(guesses)
            push!(defs, var => guesses[var])
        elseif check_defguess
            error("Invalid setup: variable $(var) has no default value or initial guess")
        end
    end

    # 4) process explicitly provided initialization equations
    if !algebraic_only
        initialization_eqs = [get_initialization_eqs(sys); initialization_eqs]
        for eq in initialization_eqs
            eq = fixpoint_sub(eq, diffmap) # expand dummy derivatives
            push!(eqs_ics, eq)
        end
    end

    # 5) process parameters as initialization unknowns
    solved_params = setup_parameter_initialization!(
        sys, pmap, defs, guesses, eqs_ics; check_defguess)

    # 6) parameter dependencies become equations, their LHS become unknowns
    # non-numeric dependent parameters stay as parameter dependencies
    new_parameter_deps = solve_parameter_dependencies!(
        sys, solved_params, eqs_ics, defs, guesses)

    # 7) handle values provided for dependent parameters similar to values for observed variables
    handle_dependent_parameter_constraints!(sys, pmap, eqs_ics)

    # parameters do not include ones that became initialization unknowns
    pars = Vector{SymbolicParam}(filter(
        !in(solved_params), parameters(sys; initial_parameters = true)))
    push!(pars, get_iv(sys))

    # 8) use observed equations for guesses of observed variables if not provided
    guessed = Set(keys(defs)) # x(t), D(x(t)), ...
    guessed = union(guessed, Set(default_toterm.(guessed))) # x(t), D(x(t)), xˍt(t), ...
    for eq in trueobs
        if !(eq.lhs in guessed)
            defs[eq.lhs] = eq.rhs
            #push!(guessed, eq.lhs) # should not encounter eq.lhs twice, so don't need to track it
        end
    end
    append!(eqs_ics, trueobs)

    vars = [vars; collect(solved_params)]

    initials = Dict(k => v for (k, v) in pmap if isinitial(k))
    merge!(defs, initials)
    isys = System(Vector{Equation}(eqs_ics),
        vars,
        pars;
        defaults = defs,
        checks = check_units,
        name,
        is_initializesystem = true,
        kwargs...)
    @set isys.parameter_dependencies = new_parameter_deps
end

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
        name = nameof(sys), kwargs...)
    eqs = equations(sys)
    trueobs, eqs = unhack_observed(observed(sys), eqs)
    vars = unique([unknowns(sys); getfield.(trueobs, :lhs)])

    eqs_ics = Equation[]
    defs = copy(defaults(sys)) # copy so we don't modify sys.defaults
    additional_guesses = anydict(guesses)
    guesses = merge(get_guesses(sys), additional_guesses)

    # PREPROCESSING
    op = anydict(op)
    u0map = anydict()
    pmap = anydict()
    build_operating_point!(sys, op, u0map, pmap, Dict(), unknowns(sys),
        parameters(sys; initial_parameters = true))
    for (k, v) in op
        if has_parameter_dependency_with_lhs(sys, k) && is_variable_floatingpoint(k)
            pmap[k] = v
        end
    end
    initsys_preprocessing!(u0map, defs)

    # Calculate valid `Initial` parameters. These are unknowns for
    # which constant initial values were provided. By this point,
    # they have been separated into `x => Initial(x)` in `u0map`
    # and `Initial(x) => val` in `pmap`.
    valid_initial_parameters = Set{BasicSymbolic}()
    for (k, v) in u0map
        isequal(Initial(k), v) || continue
        push!(valid_initial_parameters, v)
    end

    # get the initialization equations
    if !algebraic_only
        initialization_eqs = [get_initialization_eqs(sys); initialization_eqs]
    end

    # only include initialization equations where all the involved `Initial`
    # parameters are valid.
    vs = Set()
    initialization_eqs = filter(initialization_eqs) do eq
        empty!(vs)
        vars!(vs, eq; op = Initial)
        allpars = full_parameters(sys)
        for p in allpars
            if symbolic_type(p) == ArraySymbolic() &&
               Symbolics.shape(p) != Symbolics.Unknown()
                append!(allpars, Symbolics.scalarize(p))
            end
        end
        allpars = Set(allpars)
        non_params = filter(!in(allpars), vs)
        # error if non-parameters are present in the initialization equations
        if !isempty(non_params)
            throw(UnknownsInTimeIndependentInitializationError(eq, non_params))
        end
        filter!(x -> iscall(x) && isinitial(x), vs)
        invalid_initials = setdiff(vs, valid_initial_parameters)
        return isempty(invalid_initials)
    end

    append!(eqs_ics, initialization_eqs)

    # process parameters as initialization unknowns
    solved_params = setup_parameter_initialization!(
        sys, pmap, defs, guesses, eqs_ics; check_defguess)

    # parameter dependencies become equations, their LHS become unknowns
    # non-numeric dependent parameters stay as parameter dependencies
    new_parameter_deps = solve_parameter_dependencies!(
        sys, solved_params, eqs_ics, defs, guesses)

    # handle values provided for dependent parameters similar to values for observed variables
    handle_dependent_parameter_constraints!(sys, pmap, eqs_ics)

    # parameters do not include ones that became initialization unknowns
    pars = Vector{SymbolicParam}(filter(
        !in(solved_params), parameters(sys; initial_parameters = true)))
    vars = collect(solved_params)

    initials = Dict(k => v for (k, v) in pmap if isinitial(k))
    merge!(defs, initials)
    isys = System(Vector{Equation}(eqs_ics),
        vars,
        pars;
        defaults = defs,
        checks = check_units,
        name,
        is_initializesystem = true,
        kwargs...)
    @set isys.parameter_dependencies = new_parameter_deps
end

"""
    $(TYPEDSIGNATURES)

Preprocessing step for initialization. Currently removes key `k` from `defs` and `u0map`
if `k => nothing` is present in `u0map`.
"""
function initsys_preprocessing!(u0map::AbstractDict, defs::AbstractDict)
    for (k, v) in u0map
        v === nothing || continue
        delete!(defs, k)
    end
    filter_missing_values!(u0map)
end

"""
    $(TYPEDSIGNATURES)

Update `defs` and `eqs_ics` appropriately for parameter initialization. Return a dictionary
mapping solvable parameters to their `tovar` variants.
"""
function setup_parameter_initialization!(
        sys::AbstractSystem, pmap::AbstractDict, defs::AbstractDict,
        guesses::AbstractDict, eqs_ics::Vector{Equation}; check_defguess = false)
    solved_params = Set()
    for p in parameters(sys)
        if is_parameter_solvable(p, pmap, defs, guesses)
            # If either of them are `missing` the parameter is an unknown
            # But if the parameter is passed a value, use that as an additional
            # equation in the system
            _val1 = get_possibly_array_fallback_singletons(pmap, p)
            _val2 = get_possibly_array_fallback_singletons(defs, p)
            _val3 = get_possibly_array_fallback_singletons(guesses, p)
            varp = tovar(p)
            push!(solved_params, p)
            # Has a default of `missing`, and (either an equation using the value passed to `ODEProblem` or a guess)
            if _val2 === missing
                if _val1 !== nothing && _val1 !== missing
                    push!(eqs_ics, varp ~ _val1)
                    push!(defs, varp => _val1)
                elseif _val3 !== nothing
                    # assuming an equation exists (either via algebraic equations or initialization_eqs)
                    push!(defs, varp => _val3)
                elseif check_defguess
                    error("Invalid setup: parameter $(p) has no default value, initial value, or guess")
                end
                # `missing` passed to `ODEProblem`, and (either an equation using default or a guess)
            elseif _val1 === missing
                if _val2 !== nothing && _val2 !== missing
                    push!(eqs_ics, varp ~ _val2)
                    push!(defs, varp => _val2)
                elseif _val3 !== nothing
                    push!(defs, varp => _val3)
                elseif check_defguess
                    error("Invalid setup: parameter $(p) has no default value, initial value, or guess")
                end
                # given a symbolic value to ODEProblem
            elseif symbolic_type(_val1) != NotSymbolic() || is_array_of_symbolics(_val1)
                push!(eqs_ics, varp ~ _val1)
                push!(defs, varp => _val3)
                # No value passed to `ODEProblem`, but a default and a guess are present
                # _val2 !== missing is implied by it falling this far in the elseif chain
            elseif _val1 === nothing && _val2 !== nothing
                push!(eqs_ics, varp ~ _val2)
                push!(defs, varp => _val3)
            else
                # _val1 !== missing and _val1 !== nothing, so a value was provided to ODEProblem
                # This would mean `is_parameter_solvable` returned `false`, so we never end up
                # here
                error("This should never be reached")
            end
        end
    end

    return solved_params
end

"""
    $(TYPEDSIGNATURES)

Add appropriate parameter dependencies as initialization equations. Return the new list of
parameter dependencies for the initialization system.
"""
function solve_parameter_dependencies!(sys::AbstractSystem, solved_params::AbstractSet,
        eqs_ics::Vector{Equation}, defs::AbstractDict, guesses::AbstractDict)
    new_parameter_deps = Equation[]
    for eq in parameter_dependencies(sys)
        if !is_variable_floatingpoint(eq.lhs)
            push!(new_parameter_deps, eq)
            continue
        end
        varp = tovar(eq.lhs)
        push!(solved_params, eq.lhs)
        push!(eqs_ics, eq)
        guessval = get(guesses, eq.lhs, eq.rhs)
        push!(defs, varp => guessval)
    end

    return new_parameter_deps
end

"""
    $(TYPEDSIGNATURES)

Turn values provided for parameter dependencies into initialization equations.
"""
function handle_dependent_parameter_constraints!(sys::AbstractSystem, pmap::AbstractDict,
        eqs_ics::Vector{Equation})
    for (k, v) in merge(defaults(sys), pmap)
        if is_variable_floatingpoint(k) && has_parameter_dependency_with_lhs(sys, k)
            push!(eqs_ics, k ~ v)
        end
    end

    return nothing
end

"""
    $(TYPEDSIGNATURES)

Get a new symbolic variable of the same type and size as `sym`, which is a parameter.
"""
function get_initial_value_parameter(sym)
    sym = default_toterm(unwrap(sym))
    name = hasname(sym) ? getname(sym) : Symbol(sym)
    if iscall(sym) && operation(sym) === getindex
        name = Symbol(name, :_, join(arguments(sym)[2:end], "_"))
    end
    name = Symbol(name, :ₘₜₖ_₀)
    newvar = unwrap(similar_variable(sym, name; use_gensym = false))
    return toparam(newvar)
end

"""
    $(TYPEDSIGNATURES)

Given `sys` and a list of observed equations `trueobs`, remove all the equations that
directly or indirectly contain a delayed unknown of `sys`.
"""
function filter_delay_equations_variables!(sys::AbstractSystem, trueobs::Vector{Equation})
    is_time_dependent(sys) || return trueobs
    banned_vars = Set()
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
       iscall(args[1]) && get_iv(sys) in vars(args[1])
        return true
    end
    return any(x -> _has_delays(sys, x, banned), args)
end

function get_possibly_array_fallback_singletons(varmap, p)
    if haskey(varmap, p)
        return varmap[p]
    end
    if symbolic_type(p) == ArraySymbolic()
        is_sized_array_symbolic(p) || return nothing
        scal = collect(p)
        if all(x -> haskey(varmap, x), scal)
            res = [varmap[x] for x in scal]
            if any(x -> x === nothing, res)
                return nothing
            elseif any(x -> x === missing, res)
                return missing
            end
            return res
        end
    elseif iscall(p) && operation(p) == getindex
        arrp = arguments(p)[1]
        val = get_possibly_array_fallback_singletons(varmap, arrp)
        if val === nothing
            return nothing
        elseif val === missing
            return missing
        else
            return val
        end
    end
    return nothing
end

function is_parameter_solvable(p, pmap, defs, guesses)
    p = unwrap(p)
    is_variable_floatingpoint(p) || return false
    _val1 = pmap isa AbstractDict ? get_possibly_array_fallback_singletons(pmap, p) :
            nothing
    _val2 = get_possibly_array_fallback_singletons(defs, p)
    _val3 = get_possibly_array_fallback_singletons(guesses, p)
    # either (missing is a default or was passed to the ODEProblem) or (nothing was passed to
    # the ODEProblem and it has a default and a guess)
    return ((_val1 === missing || _val2 === missing) ||
            (symbolic_type(_val1) != NotSymbolic() || is_array_of_symbolics(_val1) ||
             _val1 === nothing && _val2 !== nothing)) && _val3 !== nothing
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

    if !(eltype(u0) <: Pair) && !(eltype(p) <: Pair)
        oldinitdata === nothing && return nothing

        oldinitprob = oldinitdata.initializeprob
        oldinitprob === nothing && return nothing

        meta = oldinitdata.metadata
        meta isa InitializationMetadata || return oldinitdata

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
    u0map = to_varmap(u0, dvs)
    symbols_to_symbolics!(sys, u0map)
    add_toterms!(u0map)
    pmap = to_varmap(p, ps)
    symbols_to_symbolics!(sys, pmap)
    guesses = Dict()
    defs = defaults(sys)
    use_scc = true
    initialization_eqs = Equation[]
    op = anydict()

    if oldinitdata !== nothing && oldinitdata.metadata isa InitializationMetadata
        meta = oldinitdata.metadata
        op = copy(meta.op)
        merge!(guesses, meta.guesses)
        use_scc = meta.use_scc
        initialization_eqs = meta.additional_initialization_eqs
        time_dependent_init = meta.time_dependent_init
    else
        # there is no initializeprob, so the original problem construction
        # had no solvable parameters and had the differential variables
        # specified in `u0map`.
        if u0 === missing
            # the user didn't pass `u0` to `remake`, so they want to retain
            # existing values. Fill the differential variables in `u0map`,
            # initialization will either be elided or solve for the algebraic
            # variables
            diff_idxs = isdiffeq.(equations(sys))
            for i in eachindex(dvs)
                diff_idxs[i] || continue
                u0map[dvs[i]] = newu0[i]
            end
        end
        # ensure all unknowns have guesses in case they weren't given one
        # and become solvable
        for i in eachindex(dvs)
            haskey(guesses, dvs[i]) && continue
            guesses[dvs[i]] = newu0[i]
        end
        if p === missing
            # the user didn't pass `p` to `remake`, so they want to retain
            # existing values. Fill all parameters in `pmap` so that none of
            # them are solvable.
            for p in ps
                pmap[p] = getp(sys, p)(newp)
            end
        end
        # all non-solvable parameters need values regardless
        for p in ps
            haskey(pmap, p) && continue
            is_parameter_solvable(p, pmap, defs, guesses) && continue
            pmap[p] = getp(sys, p)(newp)
        end
    end
    if t0 === nothing && is_time_dependent(sys)
        t0 = 0.0
    end
    merge!(op, u0map, pmap)
    filter_missing_values!(op)

    u0map = anydict()
    pmap = anydict()
    missing_unknowns,
    missing_pars = build_operating_point!(sys, op,
        u0map, pmap, defs, dvs, ps)
    floatT = float_type_from_varmap(op)
    u0_constructor = p_constructor = identity
    if newu0 isa StaticArray
        u0_constructor = vals -> SymbolicUtils.Code.create_array(
            typeof(newu0), floatT, Val(1), Val(length(vals)), vals...)
    end
    if newp isa StaticArray || newp isa MTKParameters && newp.initials isa StaticArray
        p_constructor = vals -> SymbolicUtils.Code.create_array(
            typeof(newp.initials), floatT, Val(1), Val(length(vals)), vals...)
    end
    kws = maybe_build_initialization_problem(
        sys, SciMLBase.isinplace(odefn), op, t0, defs, guesses,
        missing_unknowns; time_dependent_init, use_scc, initialization_eqs, floatT,
        u0_constructor, p_constructor, allow_incomplete = true, check_units = false)

    odefn = remake(odefn; kws...)
    return SciMLBase.remake_initialization_data(sys, odefn, newu0, t0, newp, newu0, newp)
end

promote_type_with_nothing(::Type{T}, ::Nothing) where {T} = T
promote_type_with_nothing(::Type{T}, ::SizedVector{0}) where {T} = T
function promote_type_with_nothing(::Type{T}, ::AbstractArray{T2}) where {T, T2}
    promote_type(T, T2)
end
function promote_type_with_nothing(::Type{T}, p::MTKParameters) where {T}
    promote_type_with_nothing(promote_type_with_nothing(T, p.tunable), p.initials)
end

promote_with_nothing(::Type, ::Nothing) = nothing
promote_with_nothing(::Type, x::SizedVector{0}) = x
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
function unhack_observed(obseqs::Vector{Equation}, eqs::Vector{Equation})
    rm_idxs = Int[]
    for (i, eq) in enumerate(obseqs)
        iscall(eq.rhs) || continue
        if operation(eq.rhs) == StructuralTransformations.change_origin
            push!(rm_idxs, i)
            continue
        end
    end

    obseqs = obseqs[setdiff(eachindex(obseqs), rm_idxs)]
    return obseqs, eqs
end

function UnknownsInTimeIndependentInitializationError(eq, non_params)
    ArgumentError("""
    Initialization equations for time-independent systems can only contain parameters. \
    Found $non_params in $eq. If the equations refer to the initial guess for unknowns, \
    use the `Initial` operator.
    """)
end
