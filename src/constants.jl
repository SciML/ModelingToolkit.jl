""""""
function isconstant(x)
    x = unwrap(x)
    x isa SymbolicT && !getmetadata(x, VariableTunable, true)
end
""""""
function toconstant(s)
    s = toparam(s)
    setmetadata(s, VariableTunable, false)
end
toconstant(s::Union{Num, Symbolics.Arr}) = wrap(toconstant(value(s)))
""""""
macro constants(xs...)
    Symbolics.parse_vars(:constants,
        Real,
        xs,
        toconstant)
end
