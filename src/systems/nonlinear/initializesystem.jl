"""
$(TYPEDSIGNATURES)

Generate `NonlinearSystem` which initializes an ODE problem from specified initial conditions of an `ODESystem`.
"""
function generate_initializesystem(sys::ODESystem;
        u0map = Dict(),
        name = nameof(sys),
        guesses = Dict(), check_defguess = false,
        default_dd_value = 0.0,
        kwargs...)
    sts, eqs = unknowns(sys), equations(sys)
    idxs_diff = isdiffeq.(eqs)
    idxs_alge = .!idxs_diff
    num_alge = sum(idxs_alge)

    # Start the equations list with algebraic equations
    eqs_ics = eqs[idxs_alge]
    u0 = Vector{Pair}(undef, 0)

    full_states = [sts; getfield.((observed(sys)), :lhs)]
    set_full_states = Set(full_states)
    guesses = todict(guesses)
    schedule = getfield(sys, :schedule)

    if schedule !== nothing
        guessmap = [x[2] => get(guesses, x[1], default_dd_value)
                    for x in schedule.dummy_sub]
        dd_guess = Dict(filter(x -> !isnothing(x[1]), guessmap))
        if u0map === nothing || isempty(u0map)
            filtered_u0 = u0map
        else
            # TODO: Don't scalarize arrays
            filtered_u0 = map(u0map) do x
                y = get(schedule.dummy_sub, x[1], x[1])
                y isa Symbolics.Arr ? collect(x[1]) .=> x[2] : x[1] => x[2]
            end
            filtered_u0 = todict(reduce(vcat, filtered_u0))
        end
    else
        dd_guess = Dict()
        filtered_u0 = u0map
    end

    defs = merge(defaults(sys), filtered_u0)
    guesses = merge(get_guesses(sys), todict(guesses), dd_guess)

    for st in full_states
        if st ∈ keys(defs)
            def = defs[st]

            if def isa Equation
                st ∉ keys(guesses) && check_defguess &&
                    error("Invalid setup: unknown $(st) has an initial condition equation with no guess.")
                push!(eqs_ics, def)
                push!(u0, st => guesses[st])
            else
                push!(eqs_ics, st ~ def)
                push!(u0, st => def)
            end
        elseif st ∈ keys(guesses)
            push!(u0, st => guesses[st])
        elseif check_defguess
            error("Invalid setup: unknown $(st) has no default value or initial guess")
        end
    end

    pars = [parameters(sys); get_iv(sys)]
    nleqs = [eqs_ics; get_initialization_eqs(sys); observed(sys)]

    sys_nl = NonlinearSystem(nleqs,
        full_states,
        pars;
        defaults = merge(ModelingToolkit.defaults(sys), todict(u0), dd_guess),
        name,
        kwargs...)

    return sys_nl
end
