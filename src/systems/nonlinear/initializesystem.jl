"""
$(TYPEDSIGNATURES)

Generate `NonlinearSystem` which initializes an ODE problem from specified initial conditions of an `ODESystem`.
"""
function generate_initializesystem(sys::ODESystem;
        u0map = Dict(),
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
        if !isnothing(u0map)
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

    pars = [parameters(sys); get_iv(sys)] # include independent variable as pseudo-parameter
    eqs_ics = [eqs_ics; observed(sys)]
    return NonlinearSystem(
        eqs_ics, vars, pars;
        defaults = defs, parameter_dependencies = parameter_dependencies(sys),
        checks = check_units,
        name, kwargs...
    )
end
