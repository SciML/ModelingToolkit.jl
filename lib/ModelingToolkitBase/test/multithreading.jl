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
