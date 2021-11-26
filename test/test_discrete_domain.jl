using ModelingToolkit
using ModelingToolkit: hashold, hassample, hasshift

@variables t x

## Shift

D1 = Shift(t; dt=0.01)
D2 = Shift(t; dt=0.01)

@test D1 == D2
@test Base.isequal(D1, D2)
@test Base.hash(D1) == Base.hash(D2)

@test D1(x) isa Num

@test isequal((D1^2)(x), D1(D1(x)))

@test hasshift(D1(x) ~ x)
@test hasshift(x ~ D1(x))
@test !hasshift(t ~ x)
@test hasshift((D1^2)(x) ~ x)
@test hasshift(D1(x) ~ D2(x))


## Sample

D1 = Sample(t; dt=0.01)
D2 = Sample(t; dt=0.01)

@test D1 == D2
@test Base.isequal(D1, D2)
@test Base.hash(D1) == Base.hash(D2)

@test D1(x) isa Num

@test hassample(D1(x) ~ x)
@test hassample(x ~ D1(x))
@test !hassample(t ~ x)
@test hassample(D1(x) ~ D2(x))


## Hold

D1 = Hold()

@test D1(x) isa Num

@test hashold(D1(x) ~ x)
@test hashold(x ~ D1(x))
@test !hashold(t ~ x)
@test hashold(D1(x) ~ D2(x))

