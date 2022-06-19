using ModelingToolkit, OrdinaryDiffEq, Test
using Random

struct MyController
    rng
    h
    MyController() = new(MersenneTwister(1234), Float32[])
end

function newdirection!(c::MyController, pars, t)
    # random update
    pars.a, = Random.rand(c.rng, 1)

    # use history for update
    pars[:b] = t + sum(c.h)
end

function updatehist!(c::MyController, sts)
    push!(c.h, sts.x)
end

function affect!(u, p, t, ctrl::MyController)
    updatehist!(ctrl, u)
    newdirection!(ctrl, p, t)
end

period = 1.0

@variables t, x(t)
@parameters a, b
D = Differential(t)

eqs = [ D(x) ~ a * x + b]

controller = MyController()

@named model = ODESystem(eqs, t, [x], [a, b], periodic_events=(1.0, affect!, controller))

sys = structural_simplify(model)

prob = ODAEProblem(sys, [x => 1.0], (0, 2.0), [a => 1.0, b => 3.0])
sol = solve(prob, Tsit5())

@test length(controller.h) == 1
@test controller.h[1] ≈ 7.873127
