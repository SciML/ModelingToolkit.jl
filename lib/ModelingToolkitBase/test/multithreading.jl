using ModelingToolkitBase, OrdinaryDiffEqTsit5, SciMLBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D, SymbolicT
using SciMLBase: EnsembleProblem, EnsembleThreads

const TRAJECTORIES = 100

@variables x(t) y(t) z(t)
@parameters sigma rho beta

eqs = [
    D(x) ~ sigma * (y - x),
    D(y) ~ x * (rho - z) - y,
    D(z) ~ x * y - beta * z,
]

@named raw_lorenz = System(eqs, t)
lorenz = mtkcompile(raw_lorenz)
u0_p0 = Dict{SymbolicT, Float64}(
    x => 1.0,
    y => 0.0,
    z => 0.0,
    sigma => 10.0,
    rho => 28.0,
    beta => 8 / 3,
)

base_prob = ODEProblem(lorenz, u0_p0, (0.0, 100.0))

function prob_func(prob, ctx)
    # Intentional symbolic `remake`
    rho_i = 20.0 + 0.01 * ctx.sim_id
    return remake(
        prob; p = Dict{SymbolicT, Float64}(
            sigma => 10.0,
            rho => rho_i,
            beta => 8 / 3,
        )
    )
end

ensemble = EnsembleProblem(
    base_prob;
    prob_func = prob_func,
    safetycopy = false,
)

@test_nowarn solve(
    ensemble, Tsit5(), EnsembleThreads();
    trajectories = TRAJECTORIES,
)

@testset "Observed functions are cached safely from several threads" begin
    @variables w(t) = 1.0
    @mtkcompile obs_sys = System([D(w) ~ -w], t)
    # A fresh problem per round, so that every round grows a new cache from several threads
    for _ in 1:20
        prob = ODEProblem(obs_sys, [], (0.0, 1.0))
        vals = zeros(200)
        Threads.@threads for i in 1:200
            vals[i] = prob[w^2 + i]
        end
        @test vals == 1.0 .+ (1:200)
    end
    prob = ODEProblem(obs_sys, [], (0.0, 1.0))
    @test ModelingToolkitBase.EnzymeCore.EnzymeRules.inactive_type(typeof(prob.f.observed))
end
