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
        name = nameof(sys), kwargs...)
    vars = unique([unknowns(sys); getfield.((observed(sys)), :lhs)])
    vars_set = Set(vars) # for efficient in-lookup

    eqs = equations(sys)
    idxs_diff = isdiffeq.(eqs)
    idxs_alge = .!idxs_diff

    # prepare map for dummy derivative substitution
    eqs_diff = eqs[idxs_diff]
    D = Differential(get_iv(sys))
    diffmap = merge(
        Dict(eq.lhs => eq.rhs for eq in eqs_diff),
        Dict(D(eq.lhs) => D(eq.rhs) for eq in observed(sys))
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
        for (y, x) in u0map
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
    for eq in parameter_dependencies(sys)
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

    eqs_ics = Symbolics.substitute.([eqs_ics; observed(sys)], (paramsubs,))
    vars = [vars; collect(values(paramsubs))]
    for k in keys(defs)
        defs[k] = substitute(defs[k], paramsubs)
    end
    return NonlinearSystem(eqs_ics,
        vars,
        pars;
        defaults = defs,
        checks = check_units,
        name,
        kwargs...)
end

function is_parameter_solvable(p, pmap, defs, guesses)
    _val1 = pmap isa AbstractDict ? get(pmap, p, nothing) : nothing
    _val2 = get(defs, p, nothing)
    _val3 = get(guesses, p, nothing)
    # either (missing is a default or was passed to the ODEProblem) or (nothing was passed to
    # the ODEProblem and it has a default and a guess)
    return ((_val1 === missing || _val2 === missing) ||
           (_val1 === nothing && _val2 !== nothing)) && _val3 !== nothing
end
