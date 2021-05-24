using ModelingToolkit
using DomainSets

domain = Interval(0, 1)
@test infimum(domain) == 0
@test supremum(domain) == 1

domain = -1.0..2.0
@test infimum(domain) == -1.0
@test supremum(domain) == 2.0

domain = Ball()
@test radius(ball) == 1.0
@test center(ball) == [0.0,0.0,0.0]

domain = Ball(2.5, [1,2,3])
@test radius(ball) == 2.5
@test center(ball) == [1,2,3]

@parameters t x
domains = [t ∈ Interval(0.0,1.0),
           x ∈ Interval(0.0,1.0)]

@parameters z
z ∈ (Interval(0.0,1.0) ⊗ Interval(0.0,1.0))

@parameters y
(x,y) ∈ CircleDomain()
@parameters r θ
(r,θ) ∈ CircleDomain(true)

(x,y) ∈ UnitDisk()
(r,θ) ∈ UnitDisk()

(x,y,z) ∈ UnitBall()
@parameters ϕ
(r,θ,ϕ) ∈ UnitBall()