using ModelingToolkit

@variables t x(t) v(t) u(t)
@derivatives D'~t

loss = (4-x)^2 + 2v^2 + u^2
eqs = [
    D(x) ~ v
    D(v) ~ u^3
]

sys = ControlSystem(loss,eqs,t,[x,v],[u],[])
dt = 0.1
tspan = (0.0,1.0)
runge_kutta_discretize(sys,dt,tspan)
