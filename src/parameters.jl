import SymbolicUtils: symtype, term, hasmetadata, issym
@enum VariableType VARIABLE PARAMETER BROWNIAN
struct MTKVariableTypeCtx end

getvariabletype(x, def = VARIABLE) = getmetadata(unwrap(x), MTKVariableTypeCtx, def)

function isparameter(x)
    x = unwrap(x)

    if x isa Symbolic && (varT = getvariabletype(x, nothing)) !== nothing
        return varT === PARAMETER
        #TODO: Delete this branch
    elseif x isa Symbolic && Symbolics.getparent(x, false) !== false
        p = Symbolics.getparent(x)
        isparameter(p) ||
            (hasmetadata(p, Symbolics.VariableSource) &&
             getmetadata(p, Symbolics.VariableSource)[1] == :parameters)
    elseif istree(x) && operation(x) isa Symbolic
        varT === PARAMETER || isparameter(operation(x))
    elseif istree(x) && operation(x) == (getindex)
        isparameter(arguments(x)[1])
    elseif x isa Symbolic
        varT === PARAMETER
    else
        false
    end
end

"""
    toparam(s)

Maps the variable to a parameter.
"""
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

"""
    tovar(s)

Maps the variable to a state.
"""
tovar(s::Symbolic) = setmetadata(s, MTKVariableTypeCtx, VARIABLE)
tovar(s::Num) = Num(tovar(value(s)))

"""
$(SIGNATURES)

Define one or more known parameters.
"""
macro parameters(xs...)
    Symbolics._parse_vars(:parameters,
                          Real,
                          xs,
                          toparam) |> esc
end
