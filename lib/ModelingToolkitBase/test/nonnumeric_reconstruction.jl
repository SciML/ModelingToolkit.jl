using ModelingToolkitBase, SciMLBase, SymbolicIndexingInterface, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D, NonnumericBuffersGetter

struct Medium
    rho::Float64
end
struct Gauge
    k::Float64
end
density(m::Medium) = m.rho
scale(g::Gauge) = g.k
@register_symbolic density(m::Medium)
@register_symbolic scale(g::Gauge)

# Two struct-valued parameters (two nonnumeric buffers of different types) and a parameter `p0`
# that initialization solves for, so the problem's parameters are rebuilt for the initialization
# problem and back.
@parameters med::Medium = Medium(2.0) g::Gauge = Gauge(3.0) p0 = missing [guess = 1.0]
@variables x(t) = 1.0
@named model = System([D(x) ~ -density(med) * scale(g) * p0 * x], t; initialization_eqs = [p0 ~ 2x])
sys = mtkcompile(model)
prob = ODEProblem(sys, [], (0.0, 1.0))
data = prob.f.initialization_data
iprob = data.initializeprob
# the initialization problem's parameters from the problem's, and the problem's from the
# initialization problem's (as from an initialization solution)
to_init = data.metadata.oop_reconstruct_u0_p.pgetter
to_prob = data.initializeprobpmap.getter

@testset "nonnumeric buffers are rebuilt with static types" begin
    @test length(prob.p.nonnumeric) == 2
    for (recon, src, dst) in ((to_init, prob, iprob), (to_prob, iprob, prob))
        @test recon.nonnumerics_fn isa NonnumericBuffersGetter
        nonnumerics = @inferred recon.nonnumerics_fn(src)
        @test typeof(nonnumerics) == typeof(parameter_values(dst).nonnumeric)
        # (the buffers of the two systems need not be in the same order)
        @test Set(nonnumerics) == Set(parameter_values(src).nonnumeric)
        @test recon(src, dst).nonnumeric == nonnumerics
    end
    @test Set(only.(to_init(prob, iprob).nonnumeric)) == Set([Medium(2.0), Gauge(3.0)])
end

@testset "nonnumeric values follow the problem" begin
    prob2 = remake(prob; p = [med => Medium(4.0)])
    @test Set(only.(to_init(prob2, iprob).nonnumeric)) == Set([Medium(4.0), Gauge(3.0)])
end
