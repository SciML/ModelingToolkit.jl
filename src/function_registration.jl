# Register functions and handle literals
Base.convert(::Type{Expression}, n::Number) = Constant(n)
macro register(fun)
    esc(:($fun(x::Expression) = Operation($fun, Expression[x])))
end
macro register2(fun)
    case0 = esc(:($fun(x::Expression, y::Expression) = Operation($fun, Expression[x, y])))
    case1 = esc(:($fun(x::Expression, y::Number) = Operation($fun, Expression[x, y])))
    case2 = esc(:($fun(x::Number, y::Expression) = Operation($fun, Expression[x, y])))
    quote
        $case0
        $case1
        $case2
    end
end

# Binary operators and functions
for fun = (:+, :-, :*, :/, :\, :%, :^, :<, :>, :(==), :!, :&, :|, :div, :rem,
           :atan2, :max, :min)
    basefun = Expr(:., Base, QuoteNode(fun))
    @eval @register2 $basefun
end


# Unary operators and functions
for fun = (:-, :sin,:sind,:sinh,:asin,:asind,:asinh,:cos,:cosd,:cosh,:acos,
           :acosd,:acosh,:cot,:cotd,:coth,:acot,:acotd,:acoth,:csc,:cscd,:csch,
           :acsc,:acscd,:acsch,:tan,:tand,:tanh,:atan,:atand,:atanh,:sec,
           :secd,:sech,:asec,:asecd,:asech,:hypot,:sqrt,:cbrt,:exp,:exp2,:expm1,
           :log,:log10,:log1p,:log2,:abs,:abs2)
    basefun = Expr(:., Base, QuoteNode(fun))
    @eval @register $basefun
end

# ifelse
Base.ifelse(cond::Expression, t::Expression, f::Expression) =
                                       Operation(ifelse, Expression[cond, t, f])
Base.ifelse(cond::Expression, t::Number, f::Expression) =
                                       Operation(ifelse, Expression[cond, t, f])
Base.ifelse(cond::Expression, t::Expression, f::Number) =
                                       Operation(ifelse, Expression[cond, t, f])
Base.ifelse(cond::Expression, t::Number, f::Number) =
                                       Operation(ifelse, Expression[cond, t, f])
