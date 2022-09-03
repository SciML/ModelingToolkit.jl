import SymbolicUtils: symtype, term, hasmetadata, issym
struct MTKConstantCtx end

function isconstant(x)
    x = unwrap(x)
    x isa Symbolic &&  getmetadata(x, MTKConstantCtx, false)
end
isconstant(x::Num) = isconstant(unwrap(x))
"""
    toconst(s::Sym)

Maps the parameter to a constant. The parameter must have a default.
"""
function toconstant(s)
    if s isa Symbolics.Arr
        Symbolics.wrap(toconstant(Symbolics.unwrap(s)))
    elseif s isa AbstractArray
        map(toconstant, s)
    else
        hasmetadata(s, Symbolics.VariableDefaultValue) || throw(ArgumentError("Constant `$(s)` must be assigned a default value."))
        setmetadata(s, MTKConstantCtx, true)
    end
end
toconstant(s::Num) = wrap(toconstant(value(s)))

"""
$(SIGNATURES)

Define one or more known variables.
"""
macro constants(xs...)
    Symbolics._parse_vars(:constants,
                          Real,
                          xs,
                          toconstant) |> esc
end
