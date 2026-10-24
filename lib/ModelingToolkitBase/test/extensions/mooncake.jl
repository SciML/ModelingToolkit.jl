using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using OrdinaryDiffEqTsit5
using SciMLSensitivity
using SymbolicIndexingInterface: setp_oop
using ForwardDiff
using Mooncake

function mooncake_gradient(f, x)
    cache = Mooncake.prepare_gradient_cache(f, x)
    _, (_, dx) = Mooncake.value_and_gradient!!(cache, f, x)
    return dx
end

@parameters α q[1:2]
@variables x(t) y(t)
@mtkcompile sys = System([D(x) ~ -α * x + q[1] * y, D(y) ~ -q[2] * y], t)
prob = ODEProblem(sys, [x => 1.0, y => 0.5, α => 1.0, q => [0.3, 1.0]], (0.0, 1.0))
loss(prob) = sum(abs2, Array(solve(prob, Tsit5(); saveat = 0.1)))

@testset "Parameters set with `setp_oop`" begin # https://github.com/SciML/SciMLSensitivity.jl/issues/1712
    set_αq1 = setp_oop(prob, [α, q[1]])
    loss_αq1(p) = loss(remake(prob; p = set_αq1(prob, p)))
    p0 = [1.0, 0.3]
    @test mooncake_gradient(loss_αq1, p0) ≈ ForwardDiff.gradient(loss_αq1, p0) rtol = 1.0e-6

    # An array parameter set as a whole.
    set_q = setp_oop(prob, q)
    loss_q(p) = loss(remake(prob; p = set_q(prob, p)))
    p0 = [0.3, 1.0]
    @test mooncake_gradient(loss_q, p0) ≈ ForwardDiff.gradient(loss_q, p0) rtol = 1.0e-6
end

@testset "`u0` set with `remake`" begin # https://github.com/SciML/ModelingToolkit.jl/issues/4180
    loss_u0(u0) = loss(remake(prob; u0))
    u0 = [1.0, 0.5]
    @test mooncake_gradient(loss_u0, u0) ≈ ForwardDiff.gradient(loss_u0, u0) rtol = 1.0e-6
end

# The value of a callable parameter sits in the non-numeric portion; a Float64 field in it,
# as in a ModelingToolkitNeuralNets network, gives it a Mooncake tangent.
struct Gain
    k::Float64
end
(gain::Gain)(x) = gain.k * x

@testset "Callable parameter" begin
    @parameters β (gain::Gain)(..) = Gain(2.0)
    @variables z(t)
    @mtkcompile gain_sys = System([D(z) ~ -β * gain(z)], t)
    gain_prob = ODEProblem(gain_sys, [z => 1.0, β => 1.0], (0.0, 1.0))
    set_β = setp_oop(gain_prob, [β])
    loss_β(p) = loss(remake(gain_prob; p = set_β(gain_prob, p)))
    @test mooncake_gradient(loss_β, [1.0]) ≈ ForwardDiff.gradient(loss_β, [1.0]) rtol = 1.0e-6
    loss_z(u0) = loss(remake(gain_prob; u0))
    @test mooncake_gradient(loss_z, [1.0]) ≈ ForwardDiff.gradient(loss_z, [1.0]) rtol = 1.0e-6
end
