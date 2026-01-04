using ModelingToolkitBase, OrdinaryDiffEq
using ModelingToolkitBase: D_nounits as D, t_nounits as t
using Test

# Standard case
@testset "Observed equation" begin
    @variables x(t) [guess = 2]
    @variables y(t)
    eqs = [D(x) ~ 1]
    initialization_eqs = [1 ~ exp(1 + x)]

    @named sys = System(eqs, t; initialization_eqs, observed = [y ~ x])
    sys = complete(sys)
    tspan = (0.0, 0.2)
    prob = ODEProblem(sys, [], tspan)

    @test prob.f.initializeprob[y] == 2.0
    @test prob.f.initializeprob[x] == 2.0
    sol = solve(prob.f.initializeprob; show_trace = Val(true))
end

@testset "Guess via parameter" begin
    @parameters a = -1.0
    @variables x(t) [guess = a]

    eqs = [D(x) ~ a]

    initialization_eqs = [1 ~ exp(1 + x)]

    @named sys = System(eqs, t; initialization_eqs)
    sys = complete(mtkcompile(sys))

    tspan = (0.0, 0.2)
    prob = ODEProblem(sys, [], tspan)

    @test prob.f.initializeprob[x] == -1.0
    sol = solve(prob.f.initializeprob; show_trace = Val(true))
end

if @isdefined(ModelingToolkit)
    @testset "Guess via observed and parameter" begin
        @parameters a = -1.0
        @variables x(t)
        @variables y(t) [guess = a]

        eqs = [D(x) ~ a]

        initialization_eqs = [1 ~ exp(1 + x)]

        @named sys = System(eqs, t; initialization_eqs, observed = [y ~ x])
        sys = complete(sys)

        tspan = (0.0, 0.2)
        prob = ODEProblem(sys, [], tspan)

        @test prob.f.initializeprob[x] == -1.0
        sol = solve(prob.f.initializeprob; show_trace = Val(true))
    end
end

# Test parameters + defaults
# https://github.com/SciML/ModelingToolkit.jl/issues/2774

@parameters x0
@variables x(t)
@variables y(t) = x
@mtkcompile sys = System([x ~ x0, D(y) ~ x], t)
prob = ODEProblem(sys, [x0 => 1.0], (0.0, 1.0))
@test prob[x] == 1.0
@test prob[y] == 1.0

@parameters x0
@variables x(t)
@variables y(t) = x0
@mtkcompile sys = System([x ~ x0, D(y) ~ x], t)
prob = ODEProblem(sys, [x0 => 1.0], (0.0, 1.0))
@test prob[x] == 1.0
@test prob[y] == 1.0

@parameters x0
@variables x(t)
@variables y(t) = x0
@mtkcompile sys = System([x ~ y, D(y) ~ x], t)
prob = ODEProblem(sys, [x0 => 1.0], (0.0, 1.0))
@test prob[x] == 1.0
@test prob[y] == 1.0

@parameters x0
@variables x(t) = x0
@variables y(t) = x
@mtkcompile sys = System([x ~ y, D(y) ~ x], t)
prob = ODEProblem(sys, [x0 => 1.0], (0.0, 1.0))
@test prob[x] == 1.0
@test prob[y] == 1.0
