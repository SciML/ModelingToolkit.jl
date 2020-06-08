# Register functions and handle literals
"""
$(SIGNATURES)

Registers a function call as a primitive for the `Operation` graph of the
ModelingToolkit IR. Example:

```julia
@register f(x,y)
```

registers `f` as a possible two-argument function.

NOTE: If registering outside of ModelingToolkit (i.e. in your own package),
this should be done at runtime (e.g. in a package `__init__()`, or inside a
method that is called at runtime and not during precompile) to ensure that
any generated functions will use your registered method. See
`inject_registered_module_functions()`.
"""
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
    if (@__MODULE__) != ModelingToolkit
        # Register external module registrations so we can rewrite any generated
        # functions through JuliaVariables.solve().
        get!(ModelingToolkit.registered_external_functions, name, @__MODULE__)
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

Base.getindex(x::Operation,i::Int64) = Operation(getindex,[x,i])
Base.one(::Operation) = 1

# ---
# Ensure that Operations that get @registered from outside the ModelingToolkit module can work without
# having to bring in the associated function into the ModelingToolkit namespace.
# We basically store information about functions registered at runtime in a ModelingToolkit variable,
# `registered_external_functions`. It's not pretty, but we are limited by the way GeneralizedGenerated
# builds a function (adding "ModelingToolkit" to every function call).
# ---
import JuliaVariables
const registered_external_functions = Dict{Symbol,Module}()
function inject_registered_module_functions(expr)
    MacroTools.postwalk(expr) do x
        MacroTools.@capture(x, f_(xs__))  # We need to find all function calls in the expression.
        # If the function call has been converted to a JuliaVariables.Var and matches
        # one of the functions we've registered...
        if !isnothing(f) && x.args[1] isa JuliaVariables.Var && x.args[1].name in keys(registered_external_functions)
            # Rewrite it from a Var to a regular function call.
            x.args[1] = getproperty(registered_external_functions[x.args[1].name], x.args[1].name)
        end
        return x  # Make sure we rebuild the expression as is.
    end
end

# TODO: Overwriting this function works, but is quite ugly. Is there a nicer way to inject the module names?
function GeneralizedGenerated.mk_function(mod::Module, ex)
    ex = macroexpand(mod, ex)
    ex = GeneralizedGenerated.simplify_ex(ex)
    ex = GeneralizedGenerated.solve(ex)

    # We need to modify the expression built by the JuliaVariables.solve(ex)
    # method, before GeneralizedGenerated.closure_conv(mod, ex) converts it to a
    # RuntimeFn (as done in GeneralizedGenerated.mk_function(mod, ex)).
    ex = inject_registered_module_functions(ex)

    fn = GeneralizedGenerated.closure_conv(mod, ex)
    if !(fn isa GeneralizedGenerated.RuntimeFn)
        error("Expect an unnamed function expression. ")
    end
    fn
end