using ModelingToolkit
using OrdinaryDiffEq
using Test

@variables t
D = Differential(t)

function oscillator_ce(k = 1.0; name)
    sts = @variables x(t)=1.0 v(t)=0.0 F(t)
    ps = @parameters k=k Θ=0.5
    eqs = [D(x) ~ v, D(v) ~ -k * x + F]
    ev = [x ~ Θ] => [x ~ 1.0, v ~ 0.0]
    ODESystem(eqs, t, sts, ps, continuous_events = [ev]; name)
end

@named oscce = oscillator_ce()
eqs = [oscce.F ~ 0]
@named eqs_sys = ODESystem(eqs, t)
@named oneosc_ce = compose(eqs_sys, oscce)
oneosc_ce_simpl = structural_simplify(oneosc_ce)

prob = ODEProblem(oneosc_ce_simpl, [], (0.0, 2.0), [])
sol = solve(prob, Tsit5(), saveat = 0.1)

@test typeof(oneosc_ce_simpl) == ODESystem
@test sol[1, 6] < 1.0 # test whether x(t) decreases over time
@test sol[1, 18] > 0.5 # test whether event happened
