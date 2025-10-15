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
    mask = trues(length(obseqs))
    for (i, eq) in enumerate(obseqs)
        mask[i] = !iscall(eq.rhs) || operation(eq.rhs) !== StructuralTransformations.change_origin
    end

    obseqs = obseqs[mask]
    return obseqs, eqs
end

function UnknownsInTimeIndependentInitializationError(eq, non_params)
    ArgumentError("""
    Initialization equations for time-independent systems can only contain parameters. \
    Found $non_params in $eq. If the equations refer to the initial guess for unknowns, \
    use the `Initial` operator.
    """)
end
