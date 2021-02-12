using ModelingToolkit, OrdinaryDiffEq, Test, LinearAlgebra

@parameters t σ ρ β
@variables x(t) y(t) z(t) k(t)
D = Differential(t)

eqs = [D(D(x)) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

sys′ = ODESystem(eqs)
sys = ode_order_lowering(sys′)

eqs2 = [0 ~ x*y - k,
        D(D(x)) ~ σ*(y-x),
        D(y) ~ x*(ρ-z)-y,
        D(z) ~ x*y - β*z]
sys2 = ODESystem(eqs2, t, [x, y, z, k], parameters(sys′))
sys2 = ode_order_lowering(sys2)
# test equation/varible ordering
ModelingToolkit.calculate_massmatrix(sys2) == Diagonal([1, 1, 1, 1, 0])

u0 = [D(x) => 2.0,
      x => 1.0,
      y => 0.0,
      z => 0.0]

p  = [σ => 28.0,
      ρ => 10.0,
      β => 8/3]

tspan = (0.0,100.0)
prob = ODEProblem(sys,u0,tspan,p,jac=true)
probexpr = ODEProblemExpr(sys,u0,tspan,p,jac=true)
sol = solve(prob,Tsit5())
solexpr = solve(eval(prob),Tsit5())
@test all(x->x==0,Array(sol - solexpr))
#using Plots; plot(sol,vars=(:x,:y))

@parameters t σ ρ β
@variables x(t) y(t) z(t)
D = Differential(t)

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

lorenz1 = ODESystem(eqs,name=:lorenz1)
lorenz2 = ODESystem(eqs,name=:lorenz2)

@variables α
@parameters γ
connections = [0 ~ lorenz1.x + lorenz2.y + α*γ]
connected = ODESystem(connections,t,[α],[γ],systems=[lorenz1,lorenz2])

u0 = [lorenz1.x => 1.0,
      lorenz1.y => 0.0,
      lorenz1.z => 0.0,
      lorenz2.x => 0.0,
      lorenz2.y => 1.0,
      lorenz2.z => 0.0,
      α => 2.0]

p  = [lorenz1.σ => 10.0,
      lorenz1.ρ => 28.0,
      lorenz1.β => 8/3,
      lorenz2.σ => 10.0,
      lorenz2.ρ => 28.0,
      lorenz2.β => 8/3,
      γ => 2.0]

tspan = (0.0,100.0)
prob = ODEProblem(connected,u0,tspan,p)
sol = solve(prob,Rodas5())
@test maximum(sol[2,:] + sol[6,:] + 2sol[1,:]) < 1e-12
#using Plots; plot(sol,vars=(:α,Symbol(lorenz1.x),Symbol(lorenz2.y)))
