"""
Test whether `x` is a constant-type Sym.
"""
function isconstant(x)
    x = unwrap(x)
    return x isa SymbolicT && !getmetadata(x, VariableTunable, true)
end

"""
    toconstant(s)

Maps the parameter to a constant. The parameter must have a default.
"""
function toconstant(s)
    s = toparam(s)
    return setmetadata(s, VariableTunable, false)
end

toconstant(s::Union{Num, Symbolics.Arr}) = wrap(toconstant(value(s)))

"""
$(SIGNATURES)

Define one or more constants.

See also [`@independent_variables`](@ref), [`@parameters`](@ref) and [`@variables`](@ref).
"""
macro constants(xs...)
    return Symbolics.parse_vars(
        :constants,
        Real,
        xs,
        toconstant
    )
end
