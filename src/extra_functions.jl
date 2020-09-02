function ifelse end
ifelse(args...) = Core.ifelse(args...)
@register Base.conj(x)
@register Base.getindex(x,i)
@register Base.binomial(n,k)
@register Base.copysign(x,y)

@register Base.signbit(x)
ModelingToolkit.derivative(::typeof(signbit), args::NTuple{1,Any}, ::Val{1}) = 0

@register Base.abs(x)
ModelingToolkit.derivative(::typeof(abs), args::NTuple{1,Any}, ::Val{1}) = ModelingToolkit.ifelse(signbit(args[1]),-one(args[1]),one(args[1]))

@register ModelingToolkit.ifelse(x,y,z)
ModelingToolkit.derivative(::typeof(ModelingToolkit.ifelse), args::NTuple{3,Any}, ::Val{1}) = 0
ModelingToolkit.derivative(::typeof(ModelingToolkit.ifelse), args::NTuple{3,Any}, ::Val{2}) = ModelingToolkit.ifelse(args[1],1,0)
ModelingToolkit.derivative(::typeof(ModelingToolkit.ifelse), args::NTuple{3,Any}, ::Val{3}) = ModelingToolkit.ifelse(args[1],0,1)
