# Solving Dynamic Optimization Problems

Systems in ModelingToolkit.jl can be directly converted to dynamic optimization or optimal control problems. In such systems, one has one or more input variables that are externally controlled to control the dynamics of the system. A dynamic optimization solves for the optimal time trajectory of the input variables in order to maximize or minimize a desired objective function. For example, a car driver might like to know how to step on the accelerator if the goal is to finish a race while using the least gas.

To begin, let us take a rocket launch example. The input variable here is the thrust exerted by the engine. The rocket state is described by its current height, mass, and velocity. The mass decreases as the rocket loses fuel while thrusting.

```@example dynamic_opt
using ModelingToolkit
t = ModelingToolkit.t_nounits
D = ModelingToolkit.D_nounits

@parameters h_c m₀ h₀ g₀ D_c c Tₘ m_c
@variables begin
    h(..)
    v(..)
    m(..), [bounds = (m_c, 1)]
    T(..), [input = true, bounds = (0, Tₘ)]
end

drag(h, v) = D_c * v^2 * exp(-h_c * (h - h₀) / h₀)
gravity(h) = g₀ * (h₀ / h)

eqs = [D(h(t)) ~ v(t),
    D(v(t)) ~ (T(t) - drag(h(t), v(t))) / m(t) - gravity(h(t)),
    D(m(t)) ~ -T(t) / c]

(ts, te) = (0.0, 0.2)
costs = [-h(te)]
cons = [T(te) ~ 0, m(te) ~ m_c]

@named rocket = System(eqs, t; costs, constraints = cons)
rocket = mtkcompile(rocket, inputs = [T(t)])

u0map = [h(t) => h₀, m(t) => m₀, v(t) => 0]
pmap = [
    g₀ => 1, m₀ => 1.0, h_c => 500, c => 0.5 * √(g₀ * h₀), D_c => 0.5 * 620 * m₀ / g₀,
    Tₘ => 3.5 * g₀ * m₀, T(t) => 0.0, h₀ => 1, m_c => 0.6]
```

What we would like to optimize here is the final height of the rocket. We do this by providing a vector of expressions corresponding to the costs. By default, the sense of the optimization is to minimize the provided cost. So to maximize the rocket height at the final time, we write `-h(te)` as the cost.

Now we can construct a problem and solve it. Let us use JuMP as our backend here. Note that the package trigger is actually [InfiniteOpt](https://infiniteopt.github.io/InfiniteOpt.jl/stable/), and not JuMP - this package includes JuMP but is designed for optimization on function spaces. Additionally we need to load the solver package - we will use [Ipopt](https://github.com/jump-dev/Ipopt.jl) here (a good choice in general).

Here we have also loaded DiffEqDevTools because we will need to construct the ODE tableau. This is only needed if one desires a custom ODE tableau for the collocation - by default the solver will use RadauIIA5.

```@example dynamic_opt
using InfiniteOpt, Ipopt, DiffEqDevTools
jprob = JuMPDynamicOptProblem(rocket, [u0map; pmap], (ts, te); dt = 0.001)
jsol = solve(jprob, JuMPCollocation(Ipopt.Optimizer, constructRadauIIA5()));
```

The solution has three fields: `jsol.sol` is the ODE solution for the states, `jsol.input_sol` is the ODE solution for the inputs, and `jsol.model` is the wrapped model that we can use to query things like objective and constraint residuals.

Let's plot the final solution and the controller here:

```@example dynamic_opt
using CairoMakie
fig = Figure(resolution = (800, 400))
ax1 = Axis(fig[1, 1], title = "Rocket trajectory", xlabel = "Time")
ax2 = Axis(fig[1, 2], title = "Control trajectory", xlabel = "Time")

for u in unknowns(rocket)
    lines!(ax1, jsol.sol.t, jsol.sol[u], label = string(u))
end
lines!(ax2, jsol.input_sol, label = "Thrust")
axislegend(ax1)
axislegend(ax2)
fig
```

### Free final time problems

There are additionally a class of dynamic optimization problems where we would like to know how to control our system to achieve something in the least time. Such problems are called free final time problems, since the final time is unknown. To model these problems in ModelingToolkit, we declare the final time as a parameter.

Below we have a model system called the double integrator. We control the acceleration of a block in order to reach a desired destination in the least time.

```@example dynamic_opt
@variables begin
    x(..)
    v(..)
    u(..), [bounds = (-1.0, 1.0), input = true]
end

@parameters tf

constr = [v(tf) ~ 0, x(tf) ~ 0]
cost = [tf] # Minimize time

@named block = System(
    [D(x(t)) ~ v(t), D(v(t)) ~ u(t)], t; costs = cost, constraints = constr)

block = mtkcompile(block; inputs = [u(t)])

u0map = [x(t) => 1.0, v(t) => 0.0]
tspan = (0.0, tf)
parammap = [u(t) => 0.0, tf => 1.0]
```

The `tf` mapping in the parameter map is treated as an initial guess.

Please note that, at the moment, free final time problems cannot support constraints defined at definite time values, like `x(3) ~ 2`.

!!! warning
    
    The Pyomo collocation methods (LagrangeRadau, LagrangeLegendre) currently are bugged for free final time problems. Strongly suggest using BackwardEuler() for such problems when using Pyomo as the backend.

When declaring the problem in this case we need to provide the number of steps, since dt can't be known in advanced. Let's solve plot our final solution and the controller for the block, using InfiniteOpt as the backend:

```@example dynamic_opt
iprob = InfiniteOptDynamicOptProblem(block, [u0map; parammap], tspan; steps = 100)
isol = solve(iprob, InfiniteOptCollocation(Ipopt.Optimizer));
```

Let's plot the final solution and the controller here:

```@example dynamic_opt
fig = Figure(resolution = (800, 400))
ax1 = Axis(fig[1, 1], title = "Block trajectory", xlabel = "Time")
ax2 = Axis(fig[1, 2], title = "Control trajectory", xlabel = "Time")

for u in unknowns(block)
    lines!(ax1, isol.sol.t, isol.sol[u], label = string(u))
end
lines!(ax2, isol.input_sol, label = "Acceleration")
axislegend(ax1)
axislegend(ax2)
fig
```
