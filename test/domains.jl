using ModelingToolkit

@parameters t x
domains = [t ∈ IntervalDomain(0.0,1.0),
           x ∈ IntervalDomain(0.0,1.0)]

@parameters z
z ∈ (IntervalDomain(0.0,1.0) ⊗ IntervalDomain(0.0,1.0))

@parameters y
(x,y) ∈ CircleDomain()
@parameters r θ
(r,θ) ∈ CircleDomain(true)
