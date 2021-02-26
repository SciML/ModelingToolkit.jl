# `_toexpr` is only used for latexify
function _toexpr(O; canonicalize=true)
    if canonicalize
        canonical, O = canonicalexpr(O)
        canonical && return O
    else
        !istree(O) && return O
    end

    op = operation(O)
    args = arguments(O)
    if op isa Differential
        ex = _toexpr(args[1]; canonicalize=canonicalize)
        wrt = _toexpr(op.x; canonicalize=canonicalize)
        return :(_derivative($ex, $wrt))
    elseif op isa Sym
        isempty(args) && return nameof(op)
        return Expr(:call, _toexpr(op; canonicalize=canonicalize), _toexpr(args; canonicalize=canonicalize)...)
    end
    return Expr(:call, op, _toexpr(args; canonicalize=canonicalize)...)
end
_toexpr(s::Sym; kw...) = nameof(s)

function canonicalexpr(O)
    !istree(O) && return true, O
    op = operation(O)
    args = arguments(O)
    if op === (^)
        if length(args) == 2 && args[2] isa Number && args[2] < 0
            ex = _toexpr(args[1])
            if args[2] == -1
                expr = Expr(:call, inv, ex)
            else
                expr = Expr(:call, ^, Expr(:call, inv, ex), -args[2])
            end
            return true, expr
        end
    end
    return false, O
end

for fun in [:toexpr, :_toexpr]
    @eval begin
        function $fun(eq::Equation; kw...)
            Expr(:(=), $fun(eq.lhs; kw...), $fun(eq.rhs; kw...))
        end

        $fun(eqs::AbstractArray; kw...) = map(eq->$fun(eq; kw...), eqs)
        $fun(x::Integer; kw...) = x
        $fun(x::AbstractFloat; kw...) = x
    end
end
_toexpr(x::Num; kw...) = _toexpr(value(x); kw...)
