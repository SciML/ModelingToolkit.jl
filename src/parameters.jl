import SymbolicUtils: symtype, term, hasmetadata, issym
@enum VariableType VARIABLE PARAMETER BROWNIAN
struct MTKVariableTypeCtx end
getvariabletype(x, def = VARIABLE) = safe_getmetadata(MTKVariableTypeCtx, unwrap(x), def)::Union{typeof(def), VariableType}
isparameter(x::Union{Num, Symbolics.Arr, Symbolics.CallAndWrap}) = isparameter(unwrap(x))
function isparameter(x::SymbolicT)
    varT = getvariabletype(x, nothing)
    return varT === PARAMETER
end
isparameter(x) = false
function toparam(s)
    if s isa Symbolics.Arr
        Symbolics.wrap(toparam(Symbolics.unwrap(s)))
    elseif s isa AbstractArray
        map(toparam, s)
    else
        setmetadata(s, MTKVariableTypeCtx, PARAMETER)
    end
end
toparam(s::Num) = wrap(toparam(value(s)))
tovar(s::SymbolicT) = setmetadata(s, MTKVariableTypeCtx, VARIABLE)
tovar(s::Union{Num, Symbolics.Arr}) = wrap(tovar(unwrap(s)))
macro parameters(xs...)
    Symbolics.parse_vars(:parameters,
        Real,
        xs,
        toparam)
end
