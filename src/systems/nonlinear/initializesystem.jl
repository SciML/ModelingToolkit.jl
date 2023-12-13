"""
$(TYPEDSIGNATURES)

Generate `NonlinearSystem` which initializes an ODE problem from specified initial conditions of an `ODESystem`.
"""
function initializesystem(sys::ODESystem; name = nameof(sys), kwargs...)
    if has_parent(sys) && (parent = get_parent(sys); parent !== nothing)
        sys = parent
    end
    sts, eqs = states(sys), equations(sys)

    idxs_diff = isdiffeq.(eqs)
    idxs_alge = .!idxs_diff

    # Algebraic equations and initial guesses are unchanged
    eqs_ics = similar(eqs)
    u0 = Vector{Any}(undef, length(sts))

    eqs_ics[idxs_alge] .= eqs[idxs_alge]
    u0[idxs_alge] .= getmetadata.(unwrap.(sts[idxs_alge]),
        Symbolics.VariableDefaultValue,
        nothing)

    for idx in findall(idxs_diff)
        st = sts[idx]
        if !hasdefault(st)
            error("Invalid setup: unknown $(st) has no default value or equation.")
        end

        def = getdefault(st)
        if def isa Equation
            if !hasguess(st)
                error("Invalid setup: unknown $(st) has an initial condition equation with no guess.")
            end
            guess = getguess(st)
            eqs_ics[idx] = def

            u0[idx] = guess
        else
            eqs_ics[idx] = st ~ def

            u0[idx] = def
        end
    end

    pars = parameters(sys)
    sys_nl = NonlinearSystem(eqs_ics,
        sts,
        pars;
        defaults = Dict(sts .=> u0),
        name,
        kwargs...)

    return sys_nl
end
