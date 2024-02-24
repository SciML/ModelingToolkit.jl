"""
$(TYPEDSIGNATURES)

Generate `NonlinearSystem` which initializes an ODE problem from specified initial conditions of an `ODESystem`.
"""
function generate_initializesystem(sys::ODESystem; name = nameof(sys), guesses = Dict(), check_defguess = false, kwargs...)
    sts, eqs = unknowns(sys), equations(sys)
    idxs_diff = isdiffeq.(eqs)
    idxs_alge = .!idxs_diff
    num_alge = sum(idxs_alge)

    # Start the equations list with algebraic equations
    eqs_ics = eqs[idxs_alge]
    u0 = Vector{Pair}(undef, 0)
    defs = defaults(sys)

    full_states = [sts; getfield.((observed(sys)), :lhs)]

    # Refactor to ODESystem construction
    # should be ModelingToolkit.guesses(sys)
     
    guesses = merge(get_guesses(sys), todict(guesses))

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

    pars = parameters(sys)
    nleqs = [eqs_ics; observed(sys)]

    sys_nl = NonlinearSystem(nleqs,
        full_states,
        pars;
        defaults = merge(ModelingToolkit.defaults(sys), todict(u0)),
        name,
        kwargs...)

    return sys_nl
end
