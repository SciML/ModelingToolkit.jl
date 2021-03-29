import SymbolicUtils: symtype, term, hasmetadata
struct MTKParameterCtx end

isparameter(x::Num) = isparameter(value(x))
isparameter(x::Symbolic) = hasmetadata(x, MTKParameterCtx) && getmetadata(x, MTKParameterCtx)
isparameter(x) = false

"""
    toparam(s::Sym)

Maps the variable to a paramter.
"""
toparam(s::Sym) = setmetadata(s, MTKParameterCtx, true)

"""
    tovar(s::Sym)

Maps the variable to a state.
"""
tovar(s::Sym) = setmetadata(s, MTKParameterCtx, false)

"""
$(SIGNATURES)

Define one or more known variables.
"""
macro parameters(xs...)
    Symbolics._parse_vars(:parameters,
                          Real,
                          xs,
                          s->setmetadata(s, MTKParameterCtx, true)) |> esc
end
