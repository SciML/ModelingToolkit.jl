using ModelingToolkit

@IVar t 
@register k(x,t)
@register f(x)  #Limit Function
@register g(x)  #Limit Function
@Deriv D'~t

lebnitztransform = ModelingToolkit.Derivative(::typeof(k(x,t)),args,::Type{Val{t}})+k(x,g(x))*ModelingToolkit.Derivative(::typeof(g(x),args,::Type{Val{x}})-k(x,f(x))*ModelingToolkit.Derivative(::typeof(f(x),args,::Type{Val{x}})

