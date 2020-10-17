@register Base.getindex(x,i::Integer)
@register Base.getindex(x,i)
@register Base.binomial(n,k)

@register Base.signbit(x)
ModelingToolkit.derivative(::typeof(signbit), args::NTuple{1,Any}, ::Val{1}) = 0

ModelingToolkit.derivative(::typeof(abs), args::NTuple{1,Any}, ::Val{1}) = IfElse.ifelse(signbit(args[1]),-one(args[1]),one(args[1]))

@register IfElse.ifelse(x,y,z::Any)
ModelingToolkit.derivative(::typeof(IfElse.ifelse), args::NTuple{3,Any}, ::Val{1}) = 0
ModelingToolkit.derivative(::typeof(IfElse.ifelse), args::NTuple{3,Any}, ::Val{2}) = IfElse.ifelse(args[1],1,0)
ModelingToolkit.derivative(::typeof(IfElse.ifelse), args::NTuple{3,Any}, ::Val{3}) = IfElse.ifelse(args[1],0,1)
