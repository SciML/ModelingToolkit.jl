using Symbolics, ModelingToolkit
using Symbolics: hasderiv, hasdiff
using ModelingToolkit: hashold, hassample, hasshift
using Test

@variables t x


op = Difference(t; dt=0.01)
@test sampletime(op) == 0.01


@test has_continuous_domain(op(x) ~ x)
@test has_continuous_domain(x ~ op(x))
@test !has_continuous_domain(t ~ x)
@test has_continuous_domain((op^2)(x) ~ x)

@test !has_discrete_domain(op(x) ~ x)
@test !has_discrete_domain(x ~ op(x))
@test !has_discrete_domain(t ~ x)
@test !has_discrete_domain((op^2)(x) ~ x)

@test !transitions_timedomain(op(x) ~ x)
@test !transitions_timedomain(x ~ op(x))


@test is_continuous_domain(op(x) ~ x)
@test !is_discrete_domain(op(x) ~ x)



## Shift

op = Shift(t; dt=0.01)
@test sampletime(op) == 0.01


@test !has_continuous_domain(op(x) ~ x)
@test !has_continuous_domain(x ~ op(x))
@test !has_continuous_domain(t ~ x)

@test has_discrete_domain(op(x) ~ x)
@test has_discrete_domain(x ~ op(x))
@test !has_discrete_domain(t ~ x)

@test !transitions_timedomain(op(x) ~ x)
@test !transitions_timedomain(x ~ op(x))

@test !is_continuous_domain(op(x) ~ x)
@test is_discrete_domain(op(x) ~ x)


## Sample

op = Sample(t; dt=0.01)
@test sampletime(op) == 0.01

@test has_continuous_domain(op(x) ~ x)
@test has_continuous_domain(x ~ op(x))

@test has_discrete_domain(op(x) ~ x)
@test has_discrete_domain(x ~ op(x))

@test transitions_timedomain(op(x) ~ x)
@test transitions_timedomain(x ~ op(x))

@test !is_continuous_domain(op(x) ~ x)
@test !is_discrete_domain(op(x) ~ x)


## Hold

op = Hold()
@test sampletime(op) === nothing

@test has_continuous_domain(op(x) ~ x)
@test has_continuous_domain(x ~ op(x))

@test has_discrete_domain(op(x) ~ x)
@test has_discrete_domain(x ~ op(x))

@test transitions_timedomain(op(x) ~ x)
@test transitions_timedomain(x ~ op(x))

@test !is_continuous_domain(op(x) ~ x)
@test !is_discrete_domain(op(x) ~ x)