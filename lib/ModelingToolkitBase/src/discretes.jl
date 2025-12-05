function todiscrete_validate(s::SymbolicT)
    if !iscall(s)
        error("""
        `@discretes` cannot create time-independent variables. Encountered $s. Use \
        `@parameters` for this purpose.
        """)
    end
    toparam(s)
end
function todiscrete_validate(s::Union{Num, Symbolics.Arr, Symbolics.CallAndWrap})
    typeof(s)(todiscrete_validate(unwrap(s)))
end

"""
$(SIGNATURES)

Define one or more discrete variables, for use in events of continuous systems. All
symbolics declare with this macro must be dependent variables.

See also [`@independent_variables`](@ref), [`@variables`](@ref) and [`@constants`](@ref).
"""
macro discretes(xs...)
    Symbolics.parse_vars(:discretes,
        Real,
        xs,
        todiscrete_validate)
end

