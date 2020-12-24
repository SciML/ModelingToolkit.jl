@register Base.getindex(x,i::Integer)
@register Base.getindex(x,i)
@register Base.binomial(n,k)

@register Base.signbit(x)
ModelingToolkit.derivative(::typeof(signbit), args::NTuple{1,Any}, ::Val{1}) = 0
ModelingToolkit.derivative(::typeof(abs), args::NTuple{1,Any}, ::Val{1}) = IfElse.ifelse(signbit(args[1]),-one(args[1]),one(args[1]))
function ModelingToolkit.derivative(::typeof(min), args::NTuple{2,Any}, ::Val{1})
    x, y = args
    IfElse.ifelse(x < y, one(x), zero(x))
end
function ModelingToolkit.derivative(::typeof(min), args::NTuple{2,Any}, ::Val{2})
    x, y = args
    IfElse.ifelse(x < y, zero(y), one(y))
end
function ModelingToolkit.derivative(::typeof(max), args::NTuple{2,Any}, ::Val{1})
    x, y = args
    IfElse.ifelse(x > y, one(x), zero(x))
end
function ModelingToolkit.derivative(::typeof(max), args::NTuple{2,Any}, ::Val{2})
    x, y = args
    IfElse.ifelse(x > y, zero(y), one(y))
end

@register IfElse.ifelse(x,y,z::Any)
ModelingToolkit.derivative(::typeof(IfElse.ifelse), args::NTuple{3,Any}, ::Val{1}) = 0
ModelingToolkit.derivative(::typeof(IfElse.ifelse), args::NTuple{3,Any}, ::Val{2}) = IfElse.ifelse(args[1],1,0)
ModelingToolkit.derivative(::typeof(IfElse.ifelse), args::NTuple{3,Any}, ::Val{3}) = IfElse.ifelse(args[1],0,1)

ModelingToolkit.@register Base.rand(x)
ModelingToolkit.@register Base.randn(x)

ModelingToolkit.@register Distributions.pdf(dist,x)
ModelingToolkit.@register Distributions.logpdf(dist,x)
ModelingToolkit.@register Distributions.cdf(dist,x)
ModelingToolkit.@register Distributions.logcdf(dist,x)
ModelingToolkit.@register Distributions.quantile(dist,x)

ModelingToolkit.@register Distributions.Uniform(mu,sigma)
ModelingToolkit.@register Distributions.Normal(mu,sigma)
