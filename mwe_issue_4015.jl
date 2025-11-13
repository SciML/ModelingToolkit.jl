using OrdinaryDiffEq

# Simple first-order linear ODE: dx/dt = (1-x)/τ
function fol!(du, u, p, t)
    τ = p[1]
    x = u[1]
    du[1] = (1 - x) / τ
end

# Parameters and initial conditions
τ = 3.0
u0 = [0.0]
p = [τ]
tspan = (0.0, 10.0)

println("\n=== Without Event ===")
prob_no_event = ODEProblem(fol!, u0, tspan, p)
@time sol_no_event = solve(prob_no_event, Tsit5(); dense=true, save_everystep=false)
println("Solution at t=10: x = ", sol_no_event.u[end][1])

println("\n=== With Continuous Event (terminate at x = 0.6) ===")

# Condition function: trigger when x crosses 0.6
function condition(u, t, integrator)
    u[1] - 0.6
end

# Affect function: terminate integration when event triggers
function affect!(integrator)
    terminate!(integrator)
end

# Create continuous callback
cb = ContinuousCallback(condition, affect!)

prob_with_event = ODEProblem(fol!, u0, tspan, p)
@time sol_with_event = solve(prob_with_event, Tsit5(); callback=cb, dense=true, save_everystep=false)
println("Solution at termination: x = ", sol_with_event.u[end][1], " at t = ", sol_with_event.t[end])

println("\n=== Timing Comparison (second run for fair comparison) ===")
@time solve(prob_no_event, Tsit5(); dense=true, save_everystep=false)
@time solve(prob_with_event, Tsit5(); callback=cb, dense=true, save_everystep=false)
