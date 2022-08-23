import SymbolicUtils: symtype, term, hasmetadata, issym
struct MTKParameterCtx end

function isparameter(x)
    x = unwrap(x)

    #TODO: Delete this branch
    if x isa Symbolic && Symbolics.getparent(x, false) !== false
        p = Symbolics.getparent(x)
        isparameter(p) ||
            (hasmetadata(p, Symbolics.VariableSource) &&
             getmetadata(p, Symbolics.VariableSource)[1] == :parameters)
    elseif istree(x) && operation(x) isa Symbolic
        getmetadata(x, MTKParameterCtx, false) ||
            isparameter(operation(x))
    elseif istree(x) && operation(x) == (getindex)
        isparameter(arguments(x)[1])
    elseif x isa Symbolic
        getmetadata(x, MTKParameterCtx, false)
    else
        false
    end
end

"""
    toparam(s::Sym)

Maps the variable to a paramter.
"""
function toparam(s)
    if s isa Symbolics.Arr
        Symbolics.wrap(toparam(Symbolics.unwrap(s)))
    elseif s isa AbstractArray
        map(toparam, s)
    else
        setmetadata(s, MTKParameterCtx, true)
    end
end
toparam(s::Num) = wrap(toparam(value(s)))

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
                          toparam) |> esc
end
