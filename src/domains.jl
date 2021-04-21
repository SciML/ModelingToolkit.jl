using DomainSets

abstract type AbstractDomain{T,N} <: Domain{T} end

struct VarDomainPairing
  variables
  domain::Domain
end
Base.:∈(variable::ModelingToolkit.Num,domain::Domain) = VarDomainPairing(value(variable),domain)
Base.:∈(variables::NTuple{N,ModelingToolkit.Num},domain::Domain) where N = VarDomainPairing(value.(variables),domain)

## Specific Domains

struct IntervalDomain{T} <: AbstractDomain{T,1}
  lower::T
  upper::T
end


struct ProductDomain{D,T,N} <: AbstractDomain{T,N}
  domains::D
end
⊗(args::AbstractDomain{T}...) where T = ProductDomain{typeof(args),T,length(args)}(args)

struct CircleDomain <: AbstractDomain{Float64,2}
  polar::Bool
  CircleDomain(polar=false) = new(polar)
end
