using ModelingToolkit
using Base.Test

d = DiscreteDomain(5) * Interval(0.0,5.0)
typeof(d) <: ProductDomain
typeof(d.domains[1]) <: DiscreteDomain
typeof(d.domains[2]) <: Interval
