using Test
using ModelingToolkit, OrdinaryDiffEq

include("../examples/rc_model.jl")

sys = structural_simplify(rc_model)
@test !isempty(ModelingToolkit.defaults(sys))
u0 = [
      capacitor.v => 0.0
      capacitor.p.i => 0.0
      resistor.v => 0.0
     ]
prob = ODEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas4())

@test sol[resistor.p.i] == sol[capacitor.p.i]
@test sol[resistor.n.i] == -sol[capacitor.p.i]
@test sol[capacitor.n.i] == -sol[capacitor.p.i]
@test iszero(sol[ground.g.i])
@test iszero(sol[ground.g.v])
@test sol[resistor.v] == sol[source.p.v] - sol[capacitor.p.v]

u0 = [
      capacitor.v => 0.0
     ]
prob = ODAEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Tsit5())

@test sol[resistor.p.i] == sol[capacitor.p.i]
@test sol[resistor.n.i] == -sol[capacitor.p.i]
@test sol[capacitor.n.i] == -sol[capacitor.p.i]
@test iszero(sol[ground.g.i])
@test iszero(sol[ground.g.v])
@test sol[resistor.v] == sol[source.p.v] - sol[capacitor.p.v]
#using Plots
#plot(sol)

include("../examples/serial_inductor.jl")
sys = structural_simplify(ll_model)
u0 = [
      inductor1.i => 0.0
      inductor2.i => 0.0
      inductor2.v => 0.0
     ]
prob = ODEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Rodas4())

prob = ODAEProblem(sys, u0, (0, 10.0))
sol = solve(prob, Tsit5())
