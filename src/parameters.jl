import SymbolicUtils: symtype, term, hasmetadata
struct MTKParameterCtx end

isparameter(x::Num) = isparameter(value(x))
isparameter(x::Symbolic) = getmetadata(x, MTKParameterCtx, false)
isparameter(x) = false

"""
    toparam(s::Sym)

Maps the variable to a paramter.
"""
function toparam(s)
    if s isa Symbolics.Arr
        Symbolics.wrap(toparam(Symbolics.unwrap(s)))
    elseif s isa AbstractArray
        map(toparam, s)
    elseif symtype(s) <: AbstractArray
        Symbolics.recurse_and_apply(toparam, s)
    else
        setmetadata(s, MTKParameterCtx, true)
    end
end
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
                          toparam,
                         ) |> esc
end
