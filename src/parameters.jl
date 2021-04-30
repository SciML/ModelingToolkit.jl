import SymbolicUtils: symtype, term, hasmetadata
struct MTKParameterCtx end

isparameter(x::Num) = isparameter(value(x))
isparameter(x::Symbolic) = getmetadata(x, MTKParameterCtx, false)
isparameter(x) = false

"""
    toparam(s::Sym)

Maps the variable to a paramter.
"""
toparam(s::Symbolic) = setmetadata(s, MTKParameterCtx, true)
toparam(s::Num) = Num(toparam(value(s)))

"""
    tovar(s::Sym)

Maps the variable to a state.
"""
tovar(s::Symbolic) = setmetadata(s, MTKParameterCtx, false)
tovar(s::Num) = Num(tovar(value(s)))

"""
$(SIGNATURES)

Define one or more known variables.
"""
macro parameters(xs...)
    Symbolics._parse_vars(:parameters,
                          Real,
                          xs,
                          x -> x isa AbstractArray ?
                              Symbolics.getindex_posthook((a, b...)->toparam(a), x) : toparam(x)
                         ) |> esc
end
