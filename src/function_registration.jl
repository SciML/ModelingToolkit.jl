import MacroTools: splitdef, combinedef
import IterTools: product

# Literals treated as constants
Base.convert(::Type{Expression}, n::Number) = Constant(n)
# Register functions and handle literals
macro register(sig)
    splitsig = splitdef(:($sig = nothing))
    name = splitsig[:name]
    args = splitsig[:args]
    typargs = typed_args(args)
    defs = :()
    for typarg in typargs
        splitsig[:args] = typarg
        splitsig[:body] = :(Operation($name, Expression[$(args...)]))
        defs = :($defs; $(combinedef(splitsig)))
    end
    esc(defs)
end
# Create all valid combinations of Expression,Number for function signature
function typed_args(args)
    cases = product(Tuple((:Expression, :Number) for i = 1:length(args))...)
    typargs = Vector{Any}[]
    for case in cases
        typarg = [Expr(:(::), args[i], case[i]) for i = 1:length(args)]
        push!(typargs, typarg)
    end
    typargs[1:end-1] # Last case is all Numbers but that should be specialized
end

# Binary operators and functions
for fun = (:+, :-, :*, :/, :\, :%, :^, :<, :>, :(==), :!, :&, :|, :div, :atan2,
           :max, :min)
    basefun = Expr(:., Base, QuoteNode(fun))
    sig = :($basefun(x,y))
    @eval @register $sig
end


# Unary operators and functions
for fun = (:-, :sin,:sind,:sinh,:asin,:asind,:asinh,:cos,:cosd,:cosh,:acos,
           :acosd,:acosh,:cot,:cotd,:coth,:acot,:acotd,:acoth,:csc,:cscd,:csch,
           :acsc,:acscd,:acsch,:tan,:tand,:tanh,:atan,:atand,:atanh,:sec,
           :secd,:sech,:asec,:asecd,:asech,:hypot,:sqrt,:cbrt,:exp,:exp2,:expm1,
           :log,:log10,:log1p,:log2,:abs,:abs2)
    basefun = Expr(:., Base, QuoteNode(fun))
    sig = :($basefun(x))
    @eval @register $sig
end

# ifelse
@register Base.ifelse(cond,t,f)
