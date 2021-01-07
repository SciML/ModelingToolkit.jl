using ModelingToolkit, OrdinaryDiffEq, Test

@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

eqs = [D(D(x)) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

sys = ODESystem(eqs)
sys = ode_order_lowering(sys)

u0 = [D(x) => 2.0,
      x => 1.0,
      y => 0.0,
      z => 0.0]

p  = [σ => 28.0,
      ρ => 10.0,
      β => 8/3]

tspan = (0.0,100.0)
prob = ODEProblem(sys,u0,tspan,p,jac=true)
sol = solve(prob,Tsit5())

sys2 = louiville_transform(sys)
@variables trJ

u0 = [D(x) => 2.0,
      x => 1.0,
      y => 0.0,
      z => 0.0,
      trJ => 0.0]

prob = ODEProblem(sys2,u0,tspan,p,jac=true)
sol = solve(prob,Tsit5())
@test sol[end,end] ≈ 366.66666666
