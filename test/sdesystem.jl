using ModelingToolkit, StaticArrays, LinearAlgebra
using StochasticDiffEq, OrdinaryDiffEq, SparseArrays
using Random, Test
using Statistics

# Define some variables
@parameters t σ ρ β
@variables x(t) y(t) z(t)
D = Differential(t)

eqs = [D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

noiseeqs = [0.1 * x,
    0.1 * y,
    0.1 * z]

# ODESystem -> SDESystem shorthand constructor
@named sys = ODESystem(eqs, t, [x, y, z], [σ, ρ, β])
@test SDESystem(sys, noiseeqs, name = :foo) isa SDESystem

@named de = SDESystem(eqs, noiseeqs, t, [x, y, z], [σ, ρ, β])
f = eval(generate_diffusion_function(de)[1])
@test f(ones(3), rand(3), nothing) == 0.1ones(3)

f = SDEFunction(de)
prob = SDEProblem(SDEFunction(de), f.g, [1.0, 0.0, 0.0], (0.0, 100.0), (10.0, 26.0, 2.33))
sol = solve(prob, SRIW1(), seed = 1)

probexpr = SDEProblem(SDEFunction(de), f.g, [1.0, 0.0, 0.0], (0.0, 100.0),
                      (10.0, 26.0, 2.33))
solexpr = solve(eval(probexpr), SRIW1(), seed = 1)

@test all(x -> x == 0, Array(sol - solexpr))

# Test no error
@test_nowarn SDEProblem(de, nothing, (0, 10.0))

noiseeqs_nd = [0.01*x 0.01*x*y 0.02*x*z
               σ 0.01*y 0.02*x*z
               ρ β 0.01*z]
@named de = SDESystem(eqs, noiseeqs_nd, t, [x, y, z], [σ, ρ, β])
f = eval(generate_diffusion_function(de)[1])
@test f([1, 2, 3.0], [0.1, 0.2, 0.3], nothing) == [0.01*1 0.01*1*2 0.02*1*3
       0.1 0.01*2 0.02*1*3
       0.2 0.3 0.01*3]

f = eval(generate_diffusion_function(de)[2])
du = ones(3, 3)
f(du, [1, 2, 3.0], [0.1, 0.2, 0.3], nothing)
@test du == [0.01*1 0.01*1*2 0.02*1*3
             0.1 0.01*2 0.02*1*3
             0.2 0.3 0.01*3]

f = SDEFunction(de)
prob = SDEProblem(SDEFunction(de), f.g, [1.0, 0.0, 0.0], (0.0, 100.0), (10.0, 26.0, 2.33),
                  noise_rate_prototype = zeros(3, 3))
sol = solve(prob, EM(), dt = 0.001)

u0map = [
    x => 1.0,
    y => 0.0,
    z => 0.0,
]

parammap = [
    σ => 10.0,
    β => 26.0,
    ρ => 2.33,
]

prob = SDEProblem(de, u0map, (0.0, 100.0), parammap)
@test size(prob.noise_rate_prototype) == (3, 3)
@test prob.noise_rate_prototype isa Matrix
sol = solve(prob, EM(), dt = 0.001)

prob = SDEProblem(de, u0map, (0.0, 100.0), parammap, sparsenoise = true)
@test size(prob.noise_rate_prototype) == (3, 3)
@test prob.noise_rate_prototype isa SparseMatrixCSC
sol = solve(prob, EM(), dt = 0.001)

# Test eval_expression=false
function test_SDEFunction_no_eval()
    # Need to test within a function scope to trigger world age issues
    f = SDEFunction(de, eval_expression = false)
    @test f([1.0, 0.0, 0.0], (10.0, 26.0, 2.33), (0.0, 100.0)) ≈ [-10.0, 26.0, 0.0]
end
test_SDEFunction_no_eval()

# modelingtoolkitize and Ito <-> Stratonovich sense
seed = 10
Random.seed!(seed)

# simple 2D diagonal noise
u0 = rand(2)
t = randn()
trange = (0.0, 100.0)
p = [1.01, 0.87]
f1!(du, u, p, t) = (du .= p[1] * u)
σ1!(du, u, p, t) = (du .= p[2] * u)

prob = SDEProblem(f1!, σ1!, u0, trange, p)
# no correction
sys = modelingtoolkitize(prob)
fdrift = eval(generate_function(sys)[1])
fdif = eval(generate_diffusion_function(sys)[1])
@test fdrift(u0, p, t) == p[1] * u0
@test fdif(u0, p, t) == p[2] * u0
fdrift! = eval(generate_function(sys)[2])
fdif! = eval(generate_diffusion_function(sys)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du == p[1] * u0
fdif!(du, u0, p, t)
@test du == p[2] * u0

# Ito -> Strat
sys2 = stochastic_integral_transform(sys, -1 // 2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0, p, t) == p[1] * u0 - 1 // 2 * p[2]^2 * u0
@test fdif(u0, p, t) == p[2] * u0
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du == p[1] * u0 - 1 // 2 * p[2]^2 * u0
fdif!(du, u0, p, t)
@test du == p[2] * u0

# Strat -> Ito
sys2 = stochastic_integral_transform(sys, 1 // 2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0, p, t) == p[1] * u0 + 1 // 2 * p[2]^2 * u0
@test fdif(u0, p, t) == p[2] * u0
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du == p[1] * u0 + 1 // 2 * p[2]^2 * u0
fdif!(du, u0, p, t)
@test du == p[2] * u0

# somewhat complicated 1D without explicit parameters but with explicit time-dependence
f2!(du, u, p, t) = (du[1] = sin(t) + cos(u[1]))
σ2!(du, u, p, t) = (du[1] = pi + atan(u[1]))

u0 = rand(1)
prob = SDEProblem(f2!, σ2!, u0, trange)
# no correction
sys = modelingtoolkitize(prob)
fdrift = eval(generate_function(sys)[1])
fdif = eval(generate_diffusion_function(sys)[1])
@test fdrift(u0, p, t) == @. sin(t) + cos(u0)
@test fdif(u0, p, t) == pi .+ atan.(u0)
fdrift! = eval(generate_function(sys)[2])
fdif! = eval(generate_diffusion_function(sys)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du == @. sin(t) + cos(u0)
fdif!(du, u0, p, t)
@test du == pi .+ atan.(u0)

# Ito -> Strat
sys2 = stochastic_integral_transform(sys, -1 // 2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0, p, t) == @. sin(t) + cos(u0) - 1 // 2 * 1 / (1 + u0^2) * (pi + atan(u0))
@test fdif(u0, p, t) == pi .+ atan.(u0)
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du == @. sin(t) + cos(u0) - 1 // 2 * 1 / (1 + u0^2) * (pi + atan(u0))
fdif!(du, u0, p, t)
@test du == pi .+ atan.(u0)

# Strat -> Ito
sys2 = stochastic_integral_transform(sys, 1 // 2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0, p, t) ≈ @. sin(t) + cos(u0) + 1 // 2 * 1 / (1 + u0^2) * (pi + atan(u0))
@test fdif(u0, p, t) == pi .+ atan.(u0)
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du ≈ @. sin(t) + cos(u0) + 1 // 2 * 1 / (1 + u0^2) * (pi + atan(u0))
fdif!(du, u0, p, t)
@test du == pi .+ atan.(u0)

# 2D diagonal noise with mixing terms (no parameters, no time-dependence)
u0 = rand(2)
t = randn()
function f3!(du, u, p, t)
    du[1] = u[1] / 2
    du[2] = u[2] / 2
    return nothing
end
function σ3!(du, u, p, t)
    du[1] = u[2]
    du[2] = u[1]
    return nothing
end

prob = SDEProblem(f3!, σ3!, u0, trange, p)
# no correction
sys = modelingtoolkitize(prob)
fdrift = eval(generate_function(sys)[1])
fdif = eval(generate_diffusion_function(sys)[1])
@test fdrift(u0, p, t) == u0 / 2
@test fdif(u0, p, t) == reverse(u0)
fdrift! = eval(generate_function(sys)[2])
fdif! = eval(generate_diffusion_function(sys)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du == u0 / 2
fdif!(du, u0, p, t)
@test du == reverse(u0)

# Ito -> Strat
sys2 = stochastic_integral_transform(sys, -1 // 2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0, p, t) == u0 * 0
@test fdif(u0, p, t) == reverse(u0)
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du == u0 * 0
fdif!(du, u0, p, t)
@test du == reverse(u0)

# Strat -> Ito
sys2 = stochastic_integral_transform(sys, 1 // 2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0, p, t) == u0
@test fdif(u0, p, t) == reverse(u0)
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du == u0
fdif!(du, u0, p, t)
@test du == reverse(u0)

# simple 2D diagonal noise oop
u0 = rand(2)
t = randn()
p = [1.01, 0.87]
f1(u, p, t) = p[1] * u
σ1(u, p, t) = p[2] * u

prob = SDEProblem(f1, σ1, u0, trange, p)
# no correction
sys = modelingtoolkitize(prob)
fdrift = eval(generate_function(sys)[1])
fdif = eval(generate_diffusion_function(sys)[1])
@test fdrift(u0, p, t) == p[1] * u0
@test fdif(u0, p, t) == p[2] * u0
fdrift! = eval(generate_function(sys)[2])
fdif! = eval(generate_diffusion_function(sys)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du == p[1] * u0
fdif!(du, u0, p, t)
@test du == p[2] * u0

# Ito -> Strat
sys2 = stochastic_integral_transform(sys, -1 // 2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0, p, t) == p[1] * u0 - 1 // 2 * p[2]^2 * u0
@test fdif(u0, p, t) == p[2] * u0
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du == p[1] * u0 - 1 // 2 * p[2]^2 * u0
fdif!(du, u0, p, t)
@test du == p[2] * u0

# Strat -> Ito
sys2 = stochastic_integral_transform(sys, 1 // 2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0, p, t) == p[1] * u0 + 1 // 2 * p[2]^2 * u0
@test fdif(u0, p, t) == p[2] * u0
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du == p[1] * u0 + 1 // 2 * p[2]^2 * u0
fdif!(du, u0, p, t)
@test du == p[2] * u0

# non-diagonal noise
u0 = rand(2)
t = randn()
p = [1.01, 0.3, 0.6, 1.2, 0.2]
f4!(du, u, p, t) = du .= p[1] * u
function g4!(du, u, p, t)
    du[1, 1] = p[2] * u[1]
    du[1, 2] = p[3] * u[1]
    du[2, 1] = p[4] * u[1]
    du[2, 2] = p[5] * u[2]
    return nothing
end

prob = SDEProblem(f4!, g4!, u0, trange, noise_rate_prototype = zeros(2, 2), p)
# no correction
sys = modelingtoolkitize(prob)
fdrift = eval(generate_function(sys)[1])
fdif = eval(generate_diffusion_function(sys)[1])
@test fdrift(u0, p, t) == p[1] * u0
@test fdif(u0, p, t) == [p[2]*u0[1] p[3]*u0[1]
                         p[4]*u0[1] p[5]*u0[2]]
fdrift! = eval(generate_function(sys)[2])
fdif! = eval(generate_diffusion_function(sys)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du == p[1] * u0
du = similar(u0, size(prob.noise_rate_prototype))
fdif!(du, u0, p, t)
@test du == [p[2]*u0[1] p[3]*u0[1]
             p[4]*u0[1] p[5]*u0[2]]

# Ito -> Strat
sys2 = stochastic_integral_transform(sys, -1 // 2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0, p, t) ≈ [
    p[1] * u0[1] - 1 // 2 * (p[2]^2 * u0[1] + p[3]^2 * u0[1]),
    p[1] * u0[2] - 1 // 2 * (p[2] * p[4] * u0[1] + p[5]^2 * u0[2]),
]
@test fdif(u0, p, t) == [p[2]*u0[1] p[3]*u0[1]
                         p[4]*u0[1] p[5]*u0[2]]
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du ≈ [
    p[1] * u0[1] - 1 // 2 * (p[2]^2 * u0[1] + p[3]^2 * u0[1]),
    p[1] * u0[2] - 1 // 2 * (p[2] * p[4] * u0[1] + p[5]^2 * u0[2]),
]
du = similar(u0, size(prob.noise_rate_prototype))
fdif!(du, u0, p, t)
@test du == [p[2]*u0[1] p[3]*u0[1]
             p[4]*u0[1] p[5]*u0[2]]

# Strat -> Ito
sys2 = stochastic_integral_transform(sys, 1 // 2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0, p, t) ≈ [
    p[1] * u0[1] + 1 // 2 * (p[2]^2 * u0[1] + p[3]^2 * u0[1]),
    p[1] * u0[2] + 1 // 2 * (p[2] * p[4] * u0[1] + p[5]^2 * u0[2]),
]
@test fdif(u0, p, t) == [p[2]*u0[1] p[3]*u0[1]
                         p[4]*u0[1] p[5]*u0[2]]
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du ≈ [
    p[1] * u0[1] + 1 // 2 * (p[2]^2 * u0[1] + p[3]^2 * u0[1]),
    p[1] * u0[2] + 1 // 2 * (p[2] * p[4] * u0[1] + p[5]^2 * u0[2]),
]
du = similar(u0, size(prob.noise_rate_prototype))
fdif!(du, u0, p, t)
@test du == [p[2]*u0[1] p[3]*u0[1]
             p[4]*u0[1] p[5]*u0[2]]

# non-diagonal noise: Torus -- Strat and Ito are identical
u0 = rand(2)
t = randn()
p = rand(1)
f5!(du, u, p, t) = du .= false
function g5!(du, u, p, t)
    du[1, 1] = cos(p[1]) * sin(u[1])
    du[1, 2] = cos(p[1]) * cos(u[1])
    du[1, 3] = -sin(p[1]) * sin(u[2])
    du[1, 4] = -sin(p[1]) * cos(u[2])
    du[2, 1] = sin(p[1]) * sin(u[1])
    du[2, 2] = sin(p[1]) * cos(u[1])
    du[2, 3] = cos(p[1]) * sin(u[2])
    du[2, 4] = cos(p[1]) * cos(u[2])
    return nothing
end

prob = SDEProblem(f5!, g5!, u0, trange, noise_rate_prototype = zeros(2, 4), p)
# no correction
sys = modelingtoolkitize(prob)
fdrift = eval(generate_function(sys)[1])
fdif = eval(generate_diffusion_function(sys)[1])
@test fdrift(u0, p, t) == 0 * u0
@test fdif(u0, p, t) ==
      [cos(p[1])*sin(u0[1]) cos(p[1])*cos(u0[1]) -sin(p[1])*sin(u0[2]) -sin(p[1])*cos(u0[2])
       sin(p[1])*sin(u0[1]) sin(p[1])*cos(u0[1]) cos(p[1])*sin(u0[2]) cos(p[1])*cos(u0[2])]
fdrift! = eval(generate_function(sys)[2])
fdif! = eval(generate_diffusion_function(sys)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du == 0 * u0
du = similar(u0, size(prob.noise_rate_prototype))
fdif!(du, u0, p, t)
@test du ==
      [cos(p[1])*sin(u0[1]) cos(p[1])*cos(u0[1]) -sin(p[1])*sin(u0[2]) -sin(p[1])*cos(u0[2])
       sin(p[1])*sin(u0[1]) sin(p[1])*cos(u0[1]) cos(p[1])*sin(u0[2]) cos(p[1])*cos(u0[2])]

# Ito -> Strat
sys2 = stochastic_integral_transform(sys, -1 // 2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0, p, t) == 0 * u0
@test fdif(u0, p, t) ==
      [cos(p[1])*sin(u0[1]) cos(p[1])*cos(u0[1]) -sin(p[1])*sin(u0[2]) -sin(p[1])*cos(u0[2])
       sin(p[1])*sin(u0[1]) sin(p[1])*cos(u0[1]) cos(p[1])*sin(u0[2]) cos(p[1])*cos(u0[2])]
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du == 0 * u0
du = similar(u0, size(prob.noise_rate_prototype))
fdif!(du, u0, p, t)
@test du ==
      [cos(p[1])*sin(u0[1]) cos(p[1])*cos(u0[1]) -sin(p[1])*sin(u0[2]) -sin(p[1])*cos(u0[2])
       sin(p[1])*sin(u0[1]) sin(p[1])*cos(u0[1]) cos(p[1])*sin(u0[2]) cos(p[1])*cos(u0[2])]

# Strat -> Ito
sys2 = stochastic_integral_transform(sys, 1 // 2)
fdrift = eval(generate_function(sys2)[1])
fdif = eval(generate_diffusion_function(sys2)[1])
@test fdrift(u0, p, t) == 0 * u0
@test fdif(u0, p, t) ==
      [cos(p[1])*sin(u0[1]) cos(p[1])*cos(u0[1]) -sin(p[1])*sin(u0[2]) -sin(p[1])*cos(u0[2])
       sin(p[1])*sin(u0[1]) sin(p[1])*cos(u0[1]) cos(p[1])*sin(u0[2]) cos(p[1])*cos(u0[2])]
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du == 0 * u0
du = similar(u0, size(prob.noise_rate_prototype))
fdif!(du, u0, p, t)
@test du ==
      [cos(p[1])*sin(u0[1]) cos(p[1])*cos(u0[1]) -sin(p[1])*sin(u0[2]) -sin(p[1])*cos(u0[2])
       sin(p[1])*sin(u0[1]) sin(p[1])*cos(u0[1]) cos(p[1])*sin(u0[2]) cos(p[1])*cos(u0[2])]

# issue #819
@testset "Combined system name collisions" begin
    @variables t
    eqs_short = [D(x) ~ σ * (y - x),
        D(y) ~ x * (ρ - z) - y,
    ]
    sys1 = SDESystem(eqs_short, noiseeqs, t, [x, y, z], [σ, ρ, β], name = :sys1)
    sys2 = SDESystem(eqs_short, noiseeqs, t, [x, y, z], [σ, ρ, β], name = :sys1)
    @test_throws ArgumentError SDESystem([sys2.y ~ sys1.z], t, [], [], [],
                                         systems = [sys1, sys2], name = :foo)
end

# observed variable handling
@variables t x(t) RHS(t)
@parameters τ
D = Differential(t)
@named fol = SDESystem([D(x) ~ (1 - x) / τ], [x], t, [x], [τ];
                       observed = [RHS ~ (1 - x) / τ])
@test isequal(RHS, @nonamespace fol.RHS)
RHS2 = RHS
@unpack RHS = fol
@test isequal(RHS, RHS2)

# issue #1644
using ModelingToolkit: rename
@variables t
eqs = [D(x) ~ x]
noiseeqs = [0.1 * x]
@named de = SDESystem(eqs, noiseeqs, t, [x], [])
@test nameof(rename(de, :newname)) == :newname

@testset "observed functionality" begin
    @parameters α β
    @variables t x(t) y(t) z(t)
    @variables weight(t)
    D = Differential(t)

    eqs = [D(x) ~ α * x]
    noiseeqs = [β * x]
    dt = 1 // 2^(7)
    x0 = 0.1

    u0map = [
        x => x0,
    ]

    parammap = [
        α => 1.5,
        β => 1.0,
    ]

    @named de = SDESystem(eqs, noiseeqs, t, [x], [α, β], observed = [weight ~ x * 10])

    prob = SDEProblem(de, u0map, (0.0, 1.0), parammap)
    sol = solve(prob, EM(), dt = dt)
    @test observed(de) == [weight ~ x * 10]
    @test sol[weight] == 10 * sol[x]

    @named ode = ODESystem(eqs, t, [x], [α, β], observed = [weight ~ x * 10])
    odeprob = ODEProblem(ode, u0map, (0.0, 1.0), parammap)
    solode = solve(odeprob, Tsit5())
    @test observed(ode) == [weight ~ x * 10]
    @test solode[weight] == 10 * solode[x]
end

@testset "Measure Transformation for variance reduction" begin
    @parameters α β
    @variables t x(t) y(t) z(t)
    D = Differential(t)

    # Evaluate Exp [(X_T)^2]
    # SDE: X_t = x + \int_0^t α X_z dz + \int_0^t b X_z dW_z
    eqs = [D(x) ~ α * x]
    noiseeqs = [β * x]

    @named de = SDESystem(eqs, noiseeqs, t, [x], [α, β])

    g(x) = x[1]^2
    dt = 1 // 2^(7)
    x0 = 0.1

    ## Standard approach
    # EM with 1`000 trajectories for stepsize 2^-7
    u0map = [
        x => x0,
    ]

    parammap = [
        α => 1.5,
        β => 1.0,
    ]

    prob = SDEProblem(de, u0map, (0.0, 1.0), parammap)

    function prob_func(prob, i, repeat)
        remake(prob, seed = seeds[i])
    end
    numtraj = Int(1e3)
    seed = 100
    Random.seed!(seed)
    seeds = rand(UInt, numtraj)

    ensemble_prob = EnsembleProblem(prob;
                                    output_func = (sol, i) -> (g(sol[end]), false),
                                    prob_func = prob_func)

    sim = solve(ensemble_prob, EM(), dt = dt, trajectories = numtraj)
    μ = mean(sim)
    σ = std(sim) / sqrt(numtraj)

    ## Variance reduction method
    u = x
    demod = ModelingToolkit.Girsanov_transform(de, u; θ0 = 0.1)

    probmod = SDEProblem(demod, u0map, (0.0, 1.0), parammap)

    ensemble_probmod = EnsembleProblem(probmod;
                                       output_func = (sol, i) -> (g(sol[x, end]) *
                                                                  sol[demod.weight, end],
                                                                  false),
                                       prob_func = prob_func)

    simmod = solve(ensemble_probmod, EM(), dt = dt, trajectories = numtraj)
    μmod = mean(simmod)
    σmod = std(simmod) / sqrt(numtraj)

    display("μ = $(round(μ, digits=2)) ± $(round(σ, digits=2))")
    display("μmod = $(round(μmod, digits=2)) ± $(round(σmod, digits=2))")

    @test μ≈μmod atol=2σ
    @test σ > σmod
end
