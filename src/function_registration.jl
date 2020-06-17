# Register functions and handle literals
"""
$(SIGNATURES)

Registers a function call as a primitive for the `Operation` graph of the
ModelingToolkit IR. Example:

```julia
@register f(x,y)
```

registers `f` as a possible two-argument function.
"""
macro register(sig)
    splitsig = splitdef(:($sig = nothing))
    name = splitsig[:name]

    # Extract the module and function name from the signature
    if name isa Symbol
        mod = __module__  # Calling module
        funcname = name
    else
        mod = name.args[1]
        funcname = name.args[2].value
    end

    args = splitsig[:args]
    typargs = typed_args(args)
    defs = :()
    for typarg in typargs
        splitsig[:args] = typarg
        if mod == (@__MODULE__)  # If the calling module is ModelingToolkit itself...
            splitsig[:body] = :(Operation($name, Expression[$(args...)]))
        else
            # Register the function's associated model so we can inject it in later.
            splitsig[:body] = quote
                get!(ModelingToolkit.registered_external_functions, Symbol($("$funcname")), $mod)
                Operation($name, Expression[$(args...)])
            end
        end
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
import DiffRules
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
Base.:^(x::Expression,y::T) where T <: Rational = Operation(Base.:^, Expression[x, y])

@register Base.conj(x)
@register Base.getindex(x,i)
@register Base.binomial(n,k)
@register Base.copysign(x,y)

Base.getindex(x::Operation,i::Int64) = Operation(getindex,[x,i])
Base.one(::Operation) = 1

# Ensure that Operations that get @registered from outside the ModelingToolkit
# module can work without having to bring in the associated function into the
# ModelingToolkit namespace. We basically store information about functions
# registered at runtime in a ModelingToolkit variable,
# `registered_external_functions`. It's not pretty, but we are limited by the
# way GeneralizedGenerated builds a function (adding "ModelingToolkit" to every
# function call).
# ---
const registered_external_functions = Dict{Symbol,Module}()
function inject_registered_module_functions(expr)
    MacroTools.postwalk(expr) do x
        # Find all function calls in the expression and extract the function
        # name and calling module.
        MacroTools.@capture(x, f_module_.f_name_(xs__))
        if isnothing(f_module)
            MacroTools.@capture(x, f_name_(xs__))
        end

        if !isnothing(f_name)
            # Set the calling module to the module that registered it.
            mod = get(registered_external_functions, f_name, f_module)
            if !isnothing(mod)
                x.args[1] = :(getproperty($mod, $(Meta.quot(f_name))))
            end
        end

        return x
    end
end