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
toparam(s::Num) = Num(toparam(value(s)))

"""
    tovar(s::Sym)

Maps the variable to a state.
"""
tovar(s::Sym) = setmetadata(s, MTKParameterCtx, false)
tovar(s::Num) = Num(tovar(value(s)))

"""
$(SIGNATURES)

Define one or more known variables.
"""
macro parameters(xs...)
    Symbolics._parse_vars(:parameters,
                          Real,
                          xs,
                          x -> x isa Array ? toparam.(x) : toparam(x)
                         ) |> esc
end
