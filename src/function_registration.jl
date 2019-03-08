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
    cases = Iterators.product(Tuple((:Expression, :Number) for i = 1:length(args))...)
    typargs = Vector{Any}[]
    for case in cases
        typarg = [Expr(:(::), args[i], case[i]) for i = 1:length(args)]
        push!(typargs, typarg)
    end
    typargs[1:end-1] # Last case is all Numbers but that should be specialized
end

# Binary & unary operators and functions
import DiffRules, SpecialFunctions, NaNMath
for (M, f, arity) in DiffRules.diffrules()
    fun = :($M.$f)
    sig = arity == 1 ? :($fun(x)) :
          arity == 2 && :($fun(x,y))
    @eval @register $sig
end

for fun ∈ [:!]
    basefun = Expr(:., Base, QuoteNode(fun))
    sig = :($basefun(x))
    @eval @register $sig
end

for fun ∈ [:<, :>, :(==), :&, :|, :div]
    basefun = Expr(:., Base, QuoteNode(fun))
    sig = :($basefun(x,y))
    @eval @register $sig
end

# ifelse
#@register Base.ifelse(cond,t,f)

# special cases
Base.:^(x::Expression,y::T) where T <: Integer = Operation(Base.:^, Expression[x, y])
