"""
Test whether `x` is a constant-type Sym.
"""
function isconstant(x)
    x = unwrap(x)
    x isa Symbolic && !getmetadata(x, VariableTunable, true)
end

"""
    toconstant(s)

Maps the parameter to a constant. The parameter must have a default.
"""
function toconstant(s)
    s = toparam(s)
    setmetadata(s, VariableTunable, false)
end

toconstant(s::Union{Num, Symbolics.Arr}) = wrap(toconstant(value(s)))

"""
$(SIGNATURES)

Define one or more constants.

See also [`@independent_variables`](@ref), [`@parameters`](@ref) and [`@variables`](@ref).
"""
macro constants(xs...)
    Symbolics._parse_vars(:constants,
        Real,
        xs,
        toconstant) |> esc
end
