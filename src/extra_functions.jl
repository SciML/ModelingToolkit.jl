@register Base.getindex(x,i::Integer) false
@register Base.getindex(x,i) false
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

IfElse.ifelse(x::Num,y,z) = Num(Term{Real}(IfElse.ifelse, [value(x), value(y), value(z)]))
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

ModelingToolkit.@register Distributions.Uniform(mu,sigma) false
ModelingToolkit.@register Distributions.Normal(mu,sigma) false

@register ∈(x::Num, y::AbstractArray)
@register ∪(x, y)
@register ∩(x, y)
@register ∨(x, y)
@register ∧(x, y)
@register ⊆(x, y)

function LinearAlgebra.det(A::AbstractMatrix{<:Num}; laplace=true)
    if laplace
        n = LinearAlgebra.checksquare(A)
        if n == 1 
            return A[1, 1]
        elseif n == 2
            return A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]
        else
            temp = 0
            # Laplace expansion along the first column
            M′ = A[:, 2:end]
            for i in axes(A, 1)
                M = M′[(1:n) .!= i, :]
                d′ = A[i, 1] * det(M)
                if iseven(i)
                    temp = iszero(temp) ? d′ : temp - d′
                else
                    temp = iszero(temp) ? d′ : temp + d′
                end
            end
        end
        return temp
    else
        if istriu(A) || istril(A)
            return det(UpperTriangular(A))
        end
        return det(lu(A; check = false))
    end
end
