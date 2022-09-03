import SymbolicUtils: symtype, term, hasmetadata, issym
struct MTKConstantCtx end

function isconstant(x)
    x = unwrap(x)
    x isa Symbolic &&  getmetadata(x, MTKConstantCtx, false)
end

"""
    toconst(s::Sym)

Maps the parameter to a constant. The parameter must have a default.
"""
function toconst(s)
    if s isa Symbolics.Arr
        Symbolics.wrap(toconst(Symbolics.unwrap(s)))
    elseif s isa AbstractArray
        map(toconst, s)
    else
        assert(hasmetadata(s,VariableDefaultValue))
        setmetadata(s, MTKConstCtx, true)
    end
end
toconst(s::Num) = wrap(toconst(value(s)))

"""
$(SIGNATURES)

Define one or more known variables.
"""
macro constants(xs...)
    Symbolics._parse_vars(:constants,
                          Real,
                          xs,
                          toconst) |> esc
end
