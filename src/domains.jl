Base.ndims(::AbstractDomain) = error("Domain dimension is undefined")

struct DiscreteDomain <: AbstractDomain
    length::Int
end

struct Interval{T} <: AbstractDomain
  start::T
  stop::T
end

struct ProductDomain{T} <: AbstractDomain
    domains::T
end

Base.:*(x::AbstractDomain...) = ProductDomain((x...))

export DiscreteDomain, Interval, ProductDomain
