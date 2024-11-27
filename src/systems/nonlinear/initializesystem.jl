"""
$(TYPEDSIGNATURES)

Generate `NonlinearSystem` which initializes an ODE problem from specified initial conditions of an `ODESystem`.
"""
function generate_initializesystem(sys::ODESystem;
        u0map = Dict(),
        pmap = Dict(),
        initialization_eqs = [],
        guesses = Dict(),
        default_dd_guess = 0.0,
        algebraic_only = false,
        check_units = true, check_defguess = false,
        implicit_dae = false,
        name = nameof(sys), kwargs...)
    trueobs, eqs = unhack_observed(observed(sys), equations(sys))
    vars = unique([unknowns(sys); getfield.(trueobs, :lhs)])

    if implicit_dae
        pre_simplification_sys = sys
        while get_parent(pre_simplification_sys) !== nothing
            pre_simplification_sys = get_parent(pre_simplification_sys)
        end
        schedule = get_schedule(sys)
        if schedule === nothing
            throw(ArgumentError("The system must be structurally simplified to create an initialization system for an implicit DAE."))
        end
        old_eqs = equations(pre_simplification_sys)
        inv_dummy_sub = Dict()
        for (k, v) in schedule.dummy_sub
            if isequal(default_toterm(k), v)
                inv_dummy_sub[v] = k
            end
        end
        new_eqs = Symbolics.fast_substitute.([trueobs; eqs], (inv_dummy_sub,))
        filter!(eq -> !isequal(eq.lhs, eq.rhs), new_eqs)
        new_sys = ODESystem(new_eqs, get_iv(sys); name = nameof(sys))
        new_sys = dummy_derivative(new_sys; to_index_zero = true, array_hack = false, cse_hack = false)
        trueobs = observed(new_sys)
        eqs = equations(new_sys)
        vars = unique([unknowns(new_sys); getfield.(trueobs, :lhs)])
    end
    vars_set = Set(vars) # for efficient in-lookup

    idxs_diff = isdiffeq.(eqs)
    idxs_alge = .!idxs_diff

    # prepare map for dummy derivative substitution
    eqs_diff = eqs[idxs_diff]
    D = Differential(get_iv(sys))
    diffmap = merge(
        Dict(eq.lhs => eq.rhs for eq in eqs_diff),
        Dict(D(eq.lhs) => D(eq.rhs) for eq in trueobs)
    )

    # 1) process dummy derivatives and u0map into initialization system
    eqs_ics = eqs[idxs_alge] # start equation list with algebraic equations
    defs = copy(defaults(sys)) # copy so we don't modify sys.defaults
    guesses = merge(get_guesses(sys), todict(guesses))
    schedule = getfield(sys, :schedule)
    if !isnothing(schedule)
        for x in filter(x -> !isnothing(x[1]), schedule.dummy_sub)
            # set dummy derivatives to default_dd_guess unless specified
            push!(defs, x[1] => get(guesses, x[1], default_dd_guess))
        end
        function process_u0map_with_dummysubs(y, x)
            y = get(schedule.dummy_sub, y, y)
            y = fixpoint_sub(y, diffmap)
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
    end

    # 2) process other variables
    for var in vars
        if var ∈ keys(defs)
            push!(eqs_ics, var ~ defs[var])
        elseif var ∈ keys(guesses)
            push!(defs, var => guesses[var])
        elseif check_defguess
            error("Invalid setup: variable $(var) has no default value or initial guess")
        end
    end

    # 3) process explicitly provided initialization equations
    if !algebraic_only
        initialization_eqs = [get_initialization_eqs(sys); initialization_eqs]
        for eq in initialization_eqs
            eq = fixpoint_sub(eq, diffmap) # expand dummy derivatives
            push!(eqs_ics, eq)
        end
    end

    # 4) process parameters as initialization unknowns
    paramsubs = Dict()
    if pmap isa SciMLBase.NullParameters
        pmap = Dict()
    end
    pmap = todict(pmap)
    for p in parameters(sys)
        if is_parameter_solvable(p, pmap, defs, guesses)
            # If either of them are `missing` the parameter is an unknown
            # But if the parameter is passed a value, use that as an additional
            # equation in the system
            _val1 = get(pmap, p, nothing)
            _val2 = get(defs, p, nothing)
            _val3 = get(guesses, p, nothing)
            varp = tovar(p)
            paramsubs[p] = varp
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
            elseif symbolic_type(_val1) != NotSymbolic()
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

    # 5) parameter dependencies become equations, their LHS become unknowns
    # non-numeric dependent parameters stay as parameter dependencies
    new_parameter_deps = Equation[]
    for eq in parameter_dependencies(sys)
        if !is_variable_floatingpoint(eq.lhs)
            push!(new_parameter_deps, eq)
            continue
        end
        varp = tovar(eq.lhs)
        paramsubs[eq.lhs] = varp
        push!(eqs_ics, eq)
        guessval = get(guesses, eq.lhs, eq.rhs)
        push!(defs, varp => guessval)
    end

    # 6) handle values provided for dependent parameters similar to values for observed variables
    for (k, v) in merge(defaults(sys), pmap)
        if is_variable_floatingpoint(k) && has_parameter_dependency_with_lhs(sys, k)
            push!(eqs_ics, paramsubs[k] ~ v)
        end
    end

    # parameters do not include ones that became initialization unknowns
    pars = vcat(
        [get_iv(sys)], # include independent variable as pseudo-parameter
        [p for p in parameters(sys) if !haskey(paramsubs, p)]
    )

    # 7) use observed equations for guesses of observed variables if not provided
    for eq in trueobs
        haskey(defs, eq.lhs) && continue
        any(x -> isequal(default_toterm(x), eq.lhs), keys(defs)) && continue

        defs[eq.lhs] = eq.rhs
    end

    eqs_ics = Symbolics.substitute.([eqs_ics; trueobs], (paramsubs,))
    vars = [vars; collect(values(paramsubs))]
    for k in keys(defs)
        defs[k] = substitute(defs[k], paramsubs)
    end
    meta = InitializationSystemMetadata(Dict{Any, Any}(u0map), Dict{Any, Any}(pmap))
    return NonlinearSystem(eqs_ics,
        vars,
        pars;
        defaults = defs,
        checks = check_units,
        parameter_dependencies = new_parameter_deps,
        name,
        metadata = meta,
        kwargs...)
end

struct InitializationSystemMetadata
    u0map::Dict{Any, Any}
    pmap::Dict{Any, Any}
end

function is_parameter_solvable(p, pmap, defs, guesses)
    p = unwrap(p)
    is_variable_floatingpoint(p) || return false
    _val1 = pmap isa AbstractDict ? get(pmap, p, nothing) : nothing
    _val2 = get(defs, p, nothing)
    _val3 = get(guesses, p, nothing)
    # either (missing is a default or was passed to the ODEProblem) or (nothing was passed to
    # the ODEProblem and it has a default and a guess)
    return ((_val1 === missing || _val2 === missing) ||
            (symbolic_type(_val1) != NotSymbolic() ||
             _val1 === nothing && _val2 !== nothing)) && _val3 !== nothing
end

function SciMLBase.remake_initializeprob(sys::ODESystem, odefn, u0, t0, p)
    if u0 === missing && p === missing
        return odefn.initializeprob, odefn.update_initializeprob!, odefn.initializeprobmap,
        odefn.initializeprobpmap
    end
    if !(eltype(u0) <: Pair) && !(eltype(p) <: Pair)
        oldinitprob = odefn.initializeprob
        if oldinitprob === nothing || !SciMLBase.has_sys(oldinitprob.f) ||
           !(oldinitprob.f.sys isa NonlinearSystem)
            return oldinitprob, odefn.update_initializeprob!, odefn.initializeprobmap,
            odefn.initializeprobpmap
        end
        pidxs = ParameterIndex[]
        pvals = []
        u0idxs = Int[]
        u0vals = []
        for sym in variable_symbols(oldinitprob)
            if is_variable(sys, sym) || has_observed_with_lhs(sys, sym)
                u0 !== missing || continue
                idx = variable_index(oldinitprob, sym)
                push!(u0idxs, idx)
                push!(u0vals, eltype(u0)(state_values(oldinitprob, idx)))
            else
                p !== missing || continue
                idx = variable_index(oldinitprob, sym)
                push!(u0idxs, idx)
                push!(u0vals, typeof(getp(sys, sym)(p))(state_values(oldinitprob, idx)))
            end
        end
        if p !== missing
            for sym in parameter_symbols(oldinitprob)
                push!(pidxs, parameter_index(oldinitprob, sym))
                if isequal(sym, get_iv(sys))
                    push!(pvals, t0)
                else
                    push!(pvals, getp(sys, sym)(p))
                end
            end
        end
        if isempty(u0idxs)
            newu0 = state_values(oldinitprob)
        else
            newu0 = remake_buffer(
                oldinitprob.f.sys, state_values(oldinitprob), u0idxs, u0vals)
        end
        if isempty(pidxs)
            newp = parameter_values(oldinitprob)
        else
            newp = remake_buffer(
                oldinitprob.f.sys, parameter_values(oldinitprob), pidxs, pvals)
        end
        initprob = remake(oldinitprob; u0 = newu0, p = newp)
        return initprob, odefn.update_initializeprob!, odefn.initializeprobmap,
        odefn.initializeprobpmap
    end
    if u0 === missing || isempty(u0)
        u0 = Dict()
    elseif !(eltype(u0) <: Pair)
        u0 = Dict(unknowns(sys) .=> u0)
    end
    if p === missing
        p = Dict()
    end
    if t0 === nothing
        t0 = 0.0
    end
    u0 = todict(u0)
    defs = defaults(sys)
    varmap = merge(defs, u0)
    for k in collect(keys(varmap))
        if varmap[k] === nothing
            delete!(varmap, k)
        end
    end
    varmap = canonicalize_varmap(varmap)
    missingvars = setdiff(unknowns(sys), collect(keys(varmap)))
    setobserved = filter(keys(varmap)) do var
        has_observed_with_lhs(sys, var) || has_observed_with_lhs(sys, default_toterm(var))
    end
    p = todict(p)
    guesses = ModelingToolkit.guesses(sys)
    solvablepars = [par
                    for par in parameters(sys)
                    if is_parameter_solvable(par, p, defs, guesses)]
    pvarmap = merge(defs, p)
    setparobserved = filter(keys(pvarmap)) do var
        has_parameter_dependency_with_lhs(sys, var)
    end
    if (((!isempty(missingvars) || !isempty(solvablepars) ||
          !isempty(setobserved) || !isempty(setparobserved)) &&
         ModelingToolkit.get_tearing_state(sys) !== nothing) ||
        !isempty(initialization_equations(sys)))
        if SciMLBase.has_initializeprob(odefn)
            oldsys = odefn.initializeprob.f.sys
            meta = get_metadata(oldsys)
            if meta isa InitializationSystemMetadata
                u0 = merge(meta.u0map, u0)
                p = merge(meta.pmap, p)
            end
        end
        for k in collect(keys(u0))
            if u0[k] === nothing
                delete!(u0, k)
            end
        end
        for k in collect(keys(p))
            if p[k] === nothing
                delete!(p, k)
            end
        end

        initprob = InitializationProblem(sys, t0, u0, p)
        initprobmap = getu(initprob, unknowns(sys))
        punknowns = [p for p in all_variable_symbols(initprob) if is_parameter(sys, p)]
        getpunknowns = getu(initprob, punknowns)
        setpunknowns = setp(sys, punknowns)
        initprobpmap = GetUpdatedMTKParameters(getpunknowns, setpunknowns)
        reqd_syms = parameter_symbols(initprob)
        update_initializeprob! = UpdateInitializeprob(
            getu(sys, reqd_syms), setu(initprob, reqd_syms))
        return initprob, update_initializeprob!, initprobmap, initprobpmap
    else
        return nothing, nothing, nothing, nothing
    end
end

"""
Counteracts the CSE/array variable hacks in `symbolics_tearing.jl` so it works with
initialization.
"""
function unhack_observed(obseqs::Vector{Equation}, eqs::Vector{Equation})
    subs = Dict()
    tempvars = Set()
    rm_idxs = Int[]
    for (i, eq) in enumerate(obseqs)
        iscall(eq.rhs) || continue
        if operation(eq.rhs) == StructuralTransformations.change_origin
            push!(rm_idxs, i)
            continue
        end
        if operation(eq.rhs) == StructuralTransformations.getindex_wrapper
            var, idxs = arguments(eq.rhs)
            subs[eq.rhs] = var[idxs...]
            push!(tempvars, var)
        end
    end

    for (i, eq) in enumerate(eqs)
        iscall(eq.rhs) || continue
        if operation(eq.rhs) == StructuralTransformations.getindex_wrapper
            var, idxs = arguments(eq.rhs)
            subs[eq.rhs] = var[idxs...]
            push!(tempvars, var)
        end
    end

    for (i, eq) in enumerate(obseqs)
        if eq.lhs in tempvars
            subs[eq.lhs] = eq.rhs
            push!(rm_idxs, i)
        end
    end

    obseqs = obseqs[setdiff(eachindex(obseqs), rm_idxs)]
    obseqs = map(obseqs) do eq
        fixpoint_sub(eq.lhs, subs) ~ fixpoint_sub(eq.rhs, subs)
    end
    eqs = map(eqs) do eq
        fixpoint_sub(eq.lhs, subs) ~ fixpoint_sub(eq.rhs, subs)
    end
    return obseqs, eqs
end
