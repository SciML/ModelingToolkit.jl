@register Base.getindex(x,i::Integer)
@register Base.getindex(x,i)
@register Base.binomial(n,k)

@register Base.signbit(x)
ModelingToolkit.derivatives[signbit] =  (args,i) -> 0

ModelingToolkit.derivatives[abs] = (args, i) -> IfElse.ifelse(signbit(args[1]),-one(args[1]),one(args[1]))

@register IfElse.ifelse(x,y,z::Any)
ModelingToolkit.derivatives[IfElse.ifelse] = Dict(1=>args -> 0,
                                                  2=>args -> IfElse.ifelse(args[1],1,0),
                                                  3=>args -> IfElse.ifelse(args[1],0,1))


ModelingToolkit.@register Base.rand(x)
ModelingToolkit.@register Base.randn(x)

ModelingToolkit.@register Distributions.pdf(dist,x)
ModelingToolkit.@register Distributions.logpdf(dist,x)
ModelingToolkit.@register Distributions.cdf(dist,x)
ModelingToolkit.@register Distributions.logcdf(dist,x)
ModelingToolkit.@register Distributions.quantile(dist,x)

ModelingToolkit.@register Distributions.Uniform(mu,sigma)
ModelingToolkit.@register Distributions.Normal(mu,sigma)
