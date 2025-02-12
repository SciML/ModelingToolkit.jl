using ModelingToolkit, StaticArrays, LinearAlgebra
using StochasticDiffEq, OrdinaryDiffEq, SparseArrays
using Random, Test
using Statistics
# imported as tt because `t` is used extensively below
using ModelingToolkit: t_nounits as tt, D_nounits as D, MTKParameters

# Define some variables
@parameters σ ρ β
@variables x(tt) y(tt) z(tt)

eqs = [D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

noiseeqs = [0.1 * x,
    0.1 * y,
    0.1 * z]

# ODESystem -> SDESystem shorthand constructor
@named sys = ODESystem(eqs, tt, [x, y, z], [σ, ρ, β])
@test SDESystem(sys, noiseeqs, name = :foo) isa SDESystem

@named de = SDESystem(eqs, noiseeqs, tt, [x, y, z], [σ, ρ, β], tspan = (0.0, 10.0))
de = complete(de)
f = eval(generate_diffusion_function(de)[1])
@test f(ones(3), (rand(3),), nothing) == 0.1ones(3)

f = SDEFunction(de)
prob = SDEProblem(de, [1.0, 0.0, 0.0], (0.0, 100.0), [10.0, 26.0, 2.33])
sol = solve(prob, SRIW1(), seed = 1)

probexpr = SDEProblem(de, [1.0, 0.0, 0.0], (0.0, 100.0),
    [10.0, 26.0, 2.33])
solexpr = solve(eval(probexpr), SRIW1(), seed = 1)

@test all(x -> x == 0, Array(sol - solexpr))

noiseeqs_nd = [0.01*x 0.01*x*y 0.02*x*z
               σ 0.01*y 0.02*x*z
               ρ β 0.01*z]
@named de = SDESystem(eqs, noiseeqs_nd, tt, [x, y, z], [σ, ρ, β])
de = complete(de)
f = eval(generate_diffusion_function(de)[1])
p = MTKParameters(de, [σ => 0.1, ρ => 0.2, β => 0.3])
@test f([1, 2, 3.0], p, nothing) == [0.01*1 0.01*1*2 0.02*1*3
       0.1 0.01*2 0.02*1*3
       0.2 0.3 0.01*3]

f = eval(generate_diffusion_function(de)[2])
du = ones(3, 3)
f(du, [1, 2, 3.0], p, nothing)
@test du == [0.01*1 0.01*1*2 0.02*1*3
             0.1 0.01*2 0.02*1*3
             0.2 0.3 0.01*3]

prob = SDEProblem(de, [1.0, 0.0, 0.0], (0.0, 100.0), (σ => 10.0, ρ => 26.0, β => 2.33),
    noise_rate_prototype = zeros(3, 3))
sol = solve(prob, EM(), dt = 0.001)

u0map = [
    x => 1.0,
    y => 0.0,
    z => 0.0
]

parammap = [
    σ => 10.0,
    β => 26.0,
    ρ => 2.33
]

prob = SDEProblem(de, u0map, (0.0, 100.0), parammap)
@test prob.f.sys === de
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
    p = MTKParameters(de, [σ => 10.0, ρ => 26.0, β => 2.33])
    @test f([1.0, 0.0, 0.0], p, (0.0, 100.0)) ≈ [-10.0, 26.0, 0.0]
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
    p[1] * u0[2] - 1 // 2 * (p[2] * p[4] * u0[1] + p[5]^2 * u0[2])
]
@test fdif(u0, p, t) == [p[2]*u0[1] p[3]*u0[1]
                         p[4]*u0[1] p[5]*u0[2]]
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du ≈ [
    p[1] * u0[1] - 1 // 2 * (p[2]^2 * u0[1] + p[3]^2 * u0[1]),
    p[1] * u0[2] - 1 // 2 * (p[2] * p[4] * u0[1] + p[5]^2 * u0[2])
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
    p[1] * u0[2] + 1 // 2 * (p[2] * p[4] * u0[1] + p[5]^2 * u0[2])
]
@test fdif(u0, p, t) == [p[2]*u0[1] p[3]*u0[1]
                         p[4]*u0[1] p[5]*u0[2]]
fdrift! = eval(generate_function(sys2)[2])
fdif! = eval(generate_diffusion_function(sys2)[2])
du = similar(u0)
fdrift!(du, u0, p, t)
@test du ≈ [
    p[1] * u0[1] + 1 // 2 * (p[2]^2 * u0[1] + p[3]^2 * u0[1]),
    p[1] * u0[2] + 1 // 2 * (p[2] * p[4] * u0[1] + p[5]^2 * u0[2])
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
    @independent_variables t
    D = Differential(t)
    eqs_short = [D(x) ~ σ * (y - x),
        D(y) ~ x * (ρ - z) - y
    ]
    noise_eqs = [y - x
                 x - y]
    sys1 = SDESystem(eqs_short, noise_eqs, t, [x, y, z], [σ, ρ, β], name = :sys1)
    sys2 = SDESystem(eqs_short, noise_eqs, t, [x, y, z], [σ, ρ, β], name = :sys1)
    @test_throws ArgumentError SDESystem([sys2.y ~ sys1.z], [sys2.y], t, [], [],
        systems = [sys1, sys2], name = :foo)
end

# observed variable handling
@variables x(tt) RHS(tt)
@parameters τ
@named fol = SDESystem([D(x) ~ (1 - x) / τ], [x], tt, [x], [τ];
    observed = [RHS ~ (1 - x) / τ])
@test isequal(RHS, @nonamespace fol.RHS)
RHS2 = RHS
@unpack RHS = fol
@test isequal(RHS, RHS2)

# issue #1644
using ModelingToolkit: rename
eqs = [D(x) ~ x]
noiseeqs = [0.1 * x]
@named de = SDESystem(eqs, noiseeqs, tt, [x], [])
@test nameof(rename(de, :newname)) == :newname

@testset "observed functionality" begin
    @parameters α β
    @variables x(tt) y(tt) z(tt)
    @variables weight(tt)

    eqs = [D(x) ~ α * x]
    noiseeqs = [β * x]
    dt = 1 // 2^(7)
    x0 = 0.1

    u0map = [
        x => x0
    ]

    parammap = [
        α => 1.5,
        β => 1.0
    ]

    @named de = SDESystem(eqs, noiseeqs, tt, [x], [α, β], observed = [weight ~ x * 10])
    de = complete(de)
    prob = SDEProblem(de, u0map, (0.0, 1.0), parammap)
    sol = solve(prob, EM(), dt = dt)
    @test observed(de) == [weight ~ x * 10]
    @test sol[weight] == 10 * sol[x]

    @named ode = ODESystem(eqs, tt, [x], [α, β], observed = [weight ~ x * 10])
    ode = complete(ode)
    odeprob = ODEProblem(ode, u0map, (0.0, 1.0), parammap)
    solode = solve(odeprob, Tsit5())
    @test observed(ode) == [weight ~ x * 10]
    @test solode[weight] == 10 * solode[x]
end

@testset "Measure Transformation for variance reduction" begin
    @parameters α β
    @variables x(tt) y(tt) z(tt)

    # Evaluate Exp [(X_T)^2]
    # SDE: X_t = x + \int_0^t α X_z dz + \int_0^t b X_z dW_z
    eqs = [D(x) ~ α * x]
    noiseeqs = [β * x]

    @named de = SDESystem(eqs, noiseeqs, tt, [x], [α, β])
    de = complete(de)
    g(x) = x[1]^2
    dt = 1 // 2^(7)
    x0 = 0.1

    ## Standard approach
    # EM with 1`000 trajectories for stepsize 2^-7
    u0map = [
        x => x0
    ]

    parammap = [
        α => 1.5,
        β => 1.0
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
        output_func = (sol, i) -> (g(sol.u[end]), false),
        prob_func = prob_func)

    sim = solve(ensemble_prob, EM(), dt = dt, trajectories = numtraj)
    μ = mean(sim)
    σ = std(sim) / sqrt(numtraj)

    ## Variance reduction method
    u = x
    demod = complete(ModelingToolkit.Girsanov_transform(de, u; θ0 = 0.1))

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

sts = @variables x(tt) y(tt) z(tt)
ps = @parameters σ ρ
@brownian β η
s = 0.001
β *= s
η *= s

eqs = [D(x) ~ σ * (y - x) + x * β,
    D(y) ~ x * (ρ - z) - y + y * β + x * η,
    D(z) ~ x * y - β * z + (x * z) * β]
@named sys1 = System(eqs, tt)
sys1 = structural_simplify(sys1)

drift_eqs = [D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y]

diffusion_eqs = [s*x 0
                 s*y s*x
                 (s * x * z)-s * z 0]

sys2 = SDESystem(drift_eqs, diffusion_eqs, tt, sts, ps, name = :sys1)
sys2 = complete(sys2)
@test sys1 == sys2

prob = SDEProblem(sys1, sts .=> [1.0, 0.0, 0.0],
    (0.0, 100.0), ps .=> (10.0, 26.0))
solve(prob, LambaEulerHeun(), seed = 1)

# Test ill-formed due to more equations than states in noise equations

@independent_variables t
@parameters p d
@variables X(t)
eqs = [D(X) ~ p - d * X]
noise_eqs = [sqrt(p), -sqrt(d * X)]
@test_throws ArgumentError SDESystem(eqs, noise_eqs, t, [X], [p, d]; name = :ssys)

noise_eqs = reshape([sqrt(p), -sqrt(d * X)], 1, 2)
ssys = SDESystem(eqs, noise_eqs, t, [X], [p, d]; name = :ssys)

# SDEProblem construction with StaticArrays
# Issue#2814
@parameters p d
@variables x(tt)
@brownian a
eqs = [D(x) ~ p - d * x + a * sqrt(p)]
@mtkbuild sys = System(eqs, tt)
u0 = @SVector[x => 10.0]
tspan = (0.0, 10.0)
ps = @SVector[p => 5.0, d => 0.5]
sprob = SDEProblem(sys, u0, tspan, ps)
@test !isinplace(sprob)
@test !isinplace(sprob.f)
@test_nowarn solve(sprob, ImplicitEM())

# Ensure diagonal noise generates vector noise function
@variables y(tt)
@brownian b
eqs = [D(x) ~ p - d * x + a * sqrt(p)
       D(y) ~ p - d * y + b * sqrt(d)]
@mtkbuild sys = System(eqs, tt)
u0 = @SVector[x => 10.0, y => 20.0]
tspan = (0.0, 10.0)
ps = @SVector[p => 5.0, d => 0.5]
sprob = SDEProblem(sys, u0, tspan, ps)
@test sprob.f.g(sprob.u0, sprob.p, sprob.tspan[1]) isa SVector{2, Float64}
@test_nowarn solve(sprob, ImplicitEM())

let
    @parameters σ ρ β
    @variables x(t) y(t) z(t)
    @brownian a
    eqs = [D(x) ~ σ * (y - x) + 0.1a * x,
        D(y) ~ x * (ρ - z) - y + 0.1a * y,
        D(z) ~ x * y - β * z + 0.1a * z]

    @mtkbuild de = System(eqs, t)

    u0map = [
        x => 1.0,
        y => 0.0,
        z => 0.0
    ]

    parammap = [
        σ => 10.0,
        β => 26.0,
        ρ => 2.33
    ]
    prob = SDEProblem(de, u0map, (0.0, 100.0), parammap)
    # TODO: re-enable this when we support scalar noise
    @test solve(prob, SOSRI()).retcode == ReturnCode.Success
end

let # test to make sure that scalar noise always receive the same kicks
    @variables x(t) y(t)
    @brownian a
    eqs = [D(x) ~ a,
        D(y) ~ a]

    @mtkbuild de = System(eqs, t)
    prob = SDEProblem(de, [x => 0, y => 0], (0.0, 10.0), [])
    sol = solve(prob, SOSRI())
    @test sol.u[end][1] == sol.u[end][2]
end

let # test that diagonal noise is correctly handled
    @parameters σ ρ β
    @variables x(t) y(t) z(t)
    @brownian a b c
    eqs = [D(x) ~ σ * (y - x) + 0.1a * x,
        D(y) ~ x * (ρ - z) - y + 0.1b * y,
        D(z) ~ x * y - β * z + 0.1c * z]

    @mtkbuild de = System(eqs, t)

    u0map = [
        x => 1.0,
        y => 0.0,
        z => 0.0
    ]

    parammap = [
        σ => 10.0,
        β => 26.0,
        ρ => 2.33
    ]

    prob = SDEProblem(de, u0map, (0.0, 100.0), parammap)
    # SOSRI only works for diagonal and scalar noise
    @test solve(prob, SOSRI()).retcode == ReturnCode.Success
end

@testset "Non-diagonal noise check" begin
    @parameters σ ρ β
    @variables x(tt) y(tt) z(tt)
    @brownian a b c d e f
    eqs = [D(x) ~ σ * (y - x) + 0.1a * x + d,
        D(y) ~ x * (ρ - z) - y + 0.1b * y + e,
        D(z) ~ x * y - β * z + 0.1c * z + f]
    @mtkbuild de = System(eqs, tt)

    u0map = [
        x => 1.0,
        y => 0.0,
        z => 0.0
    ]

    parammap = [
        σ => 10.0,
        β => 26.0,
        ρ => 2.33
    ]

    prob = SDEProblem(de, u0map, (0.0, 100.0), parammap)
    # SOSRI only works for diagonal and scalar noise
    @test_throws ErrorException solve(prob, SOSRI()).retcode==ReturnCode.Success
    # ImplicitEM does work for non-diagonal noise
    @test solve(prob, ImplicitEM()).retcode == ReturnCode.Success
    @test size(ModelingToolkit.get_noiseeqs(de)) == (3, 6)
end

@testset "Diagonal noise, less brownians than equations" begin
    @parameters σ ρ β
    @variables x(tt) y(tt) z(tt)
    @brownian a b
    eqs = [D(x) ~ σ * (y - x) + 0.1a * x,     # One brownian
        D(y) ~ x * (ρ - z) - y + 0.1b * y, # Another brownian
        D(z) ~ x * y - β * z]              # no brownians -- still diagonal
    @mtkbuild de = System(eqs, tt)

    u0map = [
        x => 1.0,
        y => 0.0,
        z => 0.0
    ]

    parammap = [
        σ => 10.0,
        β => 26.0,
        ρ => 2.33
    ]

    prob = SDEProblem(de, u0map, (0.0, 100.0), parammap)
    @test solve(prob, SOSRI()).retcode == ReturnCode.Success
end

@testset "Passing `nothing` to `u0`" begin
    @variables x(t) = 1
    @brownian b
    @mtkbuild sys = System([D(x) ~ x + b], t)
    prob = @test_nowarn SDEProblem(sys, nothing, (0.0, 1.0))
    @test_nowarn solve(prob, ImplicitEM())
end

@testset "Issue#3212: Noise dependent on observed" begin
    sts = @variables begin
        x(t) = 1.0
        input(t)
        [input = true]
    end
    ps = @parameters a = 2
    @brownian η

    eqs = [D(x) ~ -a * x + (input + 1) * η
           input ~ 0.0]

    sys = System(eqs, t, sts, ps; name = :name)
    sys = structural_simplify(sys)
    @test ModelingToolkit.get_noiseeqs(sys) ≈ [1.0]
    prob = SDEProblem(sys, [], (0.0, 1.0), [])
    @test_nowarn solve(prob, RKMil())
end

@testset "Observed variables retained after `structural_simplify`" begin
    @variables x(t) y(t) z(t)
    @brownian a
    @mtkbuild sys = System([D(x) ~ x + a, D(y) ~ y + a, z ~ x + y], t)
    @test sys isa SDESystem
    @test length(observed(sys)) == 1
    prob = SDEProblem(sys, [x => 1.0, y => 1.0], (0.0, 1.0))
    @test prob[z] ≈ 2.0
end

@testset "SDESystem to ODESystem" begin
    @variables x(t) y(t) z(t)
    @testset "Scalar noise" begin
        @named sys = SDESystem([D(x) ~ x, D(y) ~ y, z ~ x + y], [x, y, 3],
            t, [x, y, z], [], is_scalar_noise = true)
        odesys = ODESystem(sys)
        @test odesys isa ODESystem
        vs = ModelingToolkit.vars(equations(odesys))
        nbrownian = count(
            v -> ModelingToolkit.getvariabletype(v) == ModelingToolkit.BROWNIAN, vs)
        @test nbrownian == 3
        for eq in equations(odesys)
            ModelingToolkit.isdiffeq(eq) || continue
            @test length(arguments(eq.rhs)) == 4
        end
    end

    @testset "Non-scalar vector noise" begin
        @named sys = SDESystem([D(x) ~ x, D(y) ~ y, z ~ x + y], [x, y, 0],
            t, [x, y, z], [], is_scalar_noise = false)
        odesys = ODESystem(sys)
        @test odesys isa ODESystem
        vs = ModelingToolkit.vars(equations(odesys))
        nbrownian = count(
            v -> ModelingToolkit.getvariabletype(v) == ModelingToolkit.BROWNIAN, vs)
        @test nbrownian == 1
        for eq in equations(odesys)
            ModelingToolkit.isdiffeq(eq) || continue
            @test length(arguments(eq.rhs)) == 2
        end
    end

    @testset "Matrix noise" begin
        noiseeqs = [x+y y+z z+x
                    2y 2z 2x
                    z+1 x+1 y+1]
        @named sys = SDESystem([D(x) ~ x, D(y) ~ y, D(z) ~ z], noiseeqs, t, [x, y, z], [])
        odesys = ODESystem(sys)
        @test odesys isa ODESystem
        vs = ModelingToolkit.vars(equations(odesys))
        nbrownian = count(
            v -> ModelingToolkit.getvariabletype(v) == ModelingToolkit.BROWNIAN, vs)
        @test nbrownian == 3
        for eq in equations(odesys)
            @test length(arguments(eq.rhs)) == 4
        end
    end
end

@testset "`structural_simplify(::SDESystem)`" begin
    @variables x(t) y(t)
    @mtkbuild sys = SDESystem(
        [D(x) ~ x, y ~ 2x], [x, 0], t, [x, y], []; is_scalar_noise = true)
    @test sys isa SDESystem
    @test length(equations(sys)) == 1
    @test length(ModelingToolkit.get_noiseeqs(sys)) == 1
    @test length(observed(sys)) == 1
end

# Test validating types of states
@testset "Validate input types" begin
    @parameters p d
    @variables X(t)::Int64
    @brownian z
    eq2 = D(X) ~ p - d * X + z
    @test_throws ArgumentError @mtkbuild ssys = System([eq2], t)
    noiseeq = [1]
    @test_throws ArgumentError @named ssys = SDESystem([eq2], [noiseeq], t)
end

@testset "SDEFunctionExpr" begin
    @parameters σ ρ β
    @variables x(tt) y(tt) z(tt)

    eqs = [D(x) ~ σ * (y - x),
        D(y) ~ x * (ρ - z) - y,
        D(z) ~ x * y - β * z]

    noiseeqs = [0.1 * x,
        0.1 * y,
        0.1 * z]

    @named sys = ODESystem(eqs, tt, [x, y, z], [σ, ρ, β])

    @named de = SDESystem(eqs, noiseeqs, tt, [x, y, z], [σ, ρ, β], tspan = (0.0, 10.0))
    de = complete(de)

    f = SDEFunctionExpr(de)
    @test f isa Expr

    @testset "Configuration Tests" begin
        # Test with `tgrad`
        f_tgrad = SDEFunctionExpr(de; tgrad = true)
        @test f_tgrad isa Expr

        # Test with `jac`
        f_jac = SDEFunctionExpr(de; jac = true)
        @test f_jac isa Expr

        # Test with sparse Jacobian
        f_sparse = SDEFunctionExpr(de; sparse = true)
        @test f_sparse isa Expr
    end

    @testset "Ordering Tests" begin
        dvs = [z, y, x]
        ps = [β, ρ, σ]
        f_order = SDEFunctionExpr(de, dvs, ps)
        @test f_order isa Expr
    end
end

@testset "SDESystem Equality with events" begin
    @variables X(t)
    @parameters p d
    @brownian a
    seq = D(X) ~ p - d * X + a
    @mtkbuild ssys1 = System([seq], t; name = :ssys)
    @mtkbuild ssys2 = System([seq], t; name = :ssys)
    @test ssys1 == ssys2 # true

    continuous_events = [[X ~ 1.0] => [X ~ X + 5.0]]
    discrete_events = [5.0 => [d ~ d / 2.0]]

    @mtkbuild ssys1 = System([seq], t; name = :ssys, continuous_events)
    @mtkbuild ssys2 = System([seq], t; name = :ssys)
    @test ssys1 !== ssys2

    @mtkbuild ssys1 = System([seq], t; name = :ssys, discrete_events)
    @mtkbuild ssys2 = System([seq], t; name = :ssys)
    @test ssys1 !== ssys2

    @mtkbuild ssys1 = System([seq], t; name = :ssys, continuous_events)
    @mtkbuild ssys2 = System([seq], t; name = :ssys, discrete_events)
    @test ssys1 !== ssys2
end

@testset "Error when constructing SDESystem without `structural_simplify`" begin
    @parameters σ ρ β
    @variables x(tt) y(tt) z(tt)
    @brownian a
    eqs = [D(x) ~ σ * (y - x) + 0.1a * x,
        D(y) ~ x * (ρ - z) - y + 0.1a * y,
        D(z) ~ x * y - β * z + 0.1a * z]

    @named de = System(eqs, t)
    de = complete(de)

    u0map = [x => 1.0, y => 0.0, z => 0.0]
    parammap = [σ => 10.0, β => 26.0, ρ => 2.33]

    @test_throws ErrorException("SDESystem constructed by defining Brownian variables with @brownian must be simplified by calling `structural_simplify` before a SDEProblem can be constructed.") SDEProblem(
        de, u0map, (0.0, 100.0), parammap)
    de = structural_simplify(de)
    @test SDEProblem(de, u0map, (0.0, 100.0), parammap) isa SDEProblem
end
