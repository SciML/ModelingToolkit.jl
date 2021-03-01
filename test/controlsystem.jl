using ModelingToolkit, GalacticOptim, Optim

@variables t x(t) v(t) u(t)
@parameters p[1:2]
D = Differential(t)

loss = (4-x)^2 + 2v^2 + u^2
eqs = [
    D(x) ~ v - p[2]*x
    D(v) ~ p[1]*u^3 + v
]

sys = ControlSystem(loss,eqs,t,[x,v],[u],p)
dt = 0.1
tspan = (0.0,1.0)
sys = runge_kutta_discretize(sys,dt,tspan)

u0 = rand(length(states(sys))) # guess for the state values
prob = OptimizationProblem(sys,u0,[0.1,0.1],grad=true)
sol = solve(prob,BFGS())
