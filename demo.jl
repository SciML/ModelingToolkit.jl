macro mtkcompile(ex...)
    quote
        @mtkbuild $(ex...)
    end
end

function mtkcompile(args...; kwargs...)
    structural_simplify(args...; kwargs...)
end

#################################

using ModelingToolkit, OrdinaryDiffEq, StochasticDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

## ODEs

@parameters g
@variables x(t) y(t) λ(t)
eqs = [D(D(x)) ~ λ * x
       D(D(y)) ~ λ * y - g
       x^2 + y^2 ~ 1]
@mtkcompile pend = System(eqs, t)
prob = ODEProblem(pend, [x => -1, y => 0], (0.0, 10.0), [g => 1], guesses = [λ => 1])

sol = solve(prob, FBDF())

## SDEs and unified `System`

@variables x(t) y(t) z(t)
@parameters σ ρ β
@brownian a

eqs = [
    D(x) ~ σ * (y - x) + 0.1x * a,
    D(y) ~ x * (ρ - z) - y + 0.1y * a,
    D(z) ~ x * y - β * z + 0.1z * a
]

@mtkcompile sys1 = System(eqs, t)

eqs = [
    D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z
]

noiseeqs = [0.1*x;
            0.1*y;
            0.1*z;;]

@mtkcompile sys2 = SDESystem(eqs, noiseeqs, t)

u0 = [
    x => 1.0,
    y => 0.0,
    z => 0.0]

p = [σ => 28.0,
    ρ => 10.0,
    β => 8 / 3]

sdeprob = SDEProblem(sys1, u0, (0.0, 10.0), p)
sdesol = solve(sdeprob, ImplicitEM())

odeprob = ODEProblem(sys1, u0, (0.0, 10.0), p) # error!
odeprob = ODEProblem(sys1, u0, (0.0, 10.0), p; check_compatibility = false)

@variables x y z
@parameters σ ρ β

# Define a nonlinear system
eqs = [0 ~ σ * (y - x),
    y ~ x * (ρ - z),
    β * z ~ x * y]
@mtkcompile sys = System(eqs)

## ImplicitDiscrete Affects

@parameters g
@variables x(t) y(t) λ(t)
eqs = [D(D(x)) ~ λ * x
       D(D(y)) ~ λ * y - g
       x^2 + y^2 ~ 1]
c_evt = [t ~ 5.0] => [x ~ Pre(x) + 0.1]
@mtkcompile pend = System(eqs, t, continuous_events = c_evt)
prob = ODEProblem(pend, [x => -1, y => 0], (0.0, 10.0), [g => 1], guesses = [λ => 1])

sol = solve(prob, FBDF())

## `@named` and `ParentScope`

function SysA(; name, var1)
    @variables x(t)
    return System([D(x) ~ var1], t; name)
end
function SysB(; name, var1)
    @variables x(t)
    @named subsys = SysA(; var1)
    return System([D(x) ~ x], t; systems = [subsys], name)
end
function SysC(; name)
    @variables x(t)
    @named subsys = SysB(; var1 = x)
    return System([D(x) ~ x], t; systems = [subsys], name)
end
@mtkcompile sys = SysC()
