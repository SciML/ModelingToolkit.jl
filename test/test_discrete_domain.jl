using ModelingToolkit
using ModelingToolkit: hashold, hassample, hasshift, value

@variables t x

## Shift

D1 = Shift(t)
D2 = Shift(t)

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

@test D1.steps == 1
@test (D1^2).steps == 2
@test (D1^(-2)).steps == -2


## Sample

D1 = Sample(t, 0.01)
D2 = Sample(t, 0.01)

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


## Test ShiftIndex

@variables t t2 x(t) y(t) y2(t2) y3(t, t2)
k = ShiftIndex(t, 0.1)


xk = x(k)
@test xk isa Num
@test ModelingToolkit.has_time_domain(xk)
@test ModelingToolkit.get_time_domain(xk) == Clock(t, 0.1)

k1 = k+1
@test k1.steps == 1
@test k1.clock === k.clock
@test sampletime(k1) === sampletime(k)


xk1 = x(k+1)
@test xk1 isa Num
shift = value(xk1)
@test shift isa Term
@test value(xk1).f isa Shift

# more complicated expressions
dt = 0.1
z = Shift(t)
x2 = x^2
@test isequal(x2(k+1), z(x2))

@test_throws ErrorException y2(k)

@test_throws ErrorException (x + y)(k)
@test_throws ErrorException y3(k)

# test that the sample operator is present when continuous variables are indexed
d = Clock(t, 1)
@variables t xd(t) [timedomain=d]
k = ShiftIndex(d, 0)
eq = x(k+1) ~ x(k) + x(k-1)
@test !ModelingToolkit.hassample(eq)
@test  ModelingToolkit.hasshift(eq)
