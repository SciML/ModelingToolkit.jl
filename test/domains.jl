using ModelingToolkit, DomainSets, Test

domain = Interval(0, 1)
@test infimum(domain) == 0
@test supremum(domain) == 1

domain = -1.0..2.0
@test infimum(domain) == -1.0
@test supremum(domain) == 2.0

ball = Ball()
@test radius(ball) == 1.0
@test center(ball) == [0.0,0.0,0.0]

ball = Ball(2.5, [1,2])
@test radius(ball) == 2.5
@test center(ball) == [1,2]

@parameters t x
domains = [t ∈ Interval(0.0,1.0),
           x ∈ Interval(0.0,1.0)]

@parameters y, z
(x,y) ∈ Ball(2.0, [0,0])
(x,y,z) ∈ Ball(1.5, [1,2,3])