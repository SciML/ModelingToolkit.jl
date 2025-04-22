using ModelingToolkit, DelayDiffEq, StaticArrays, Test
using SymbolicIndexingInterface: is_markovian
using ModelingToolkit: t_nounits as t, D_nounits as D

p0 = 0.2;
q0 = 0.3;
v0 = 1;
d0 = 5;
p1 = 0.2;
q1 = 0.3;
v1 = 1;
d1 = 1;
d2 = 1;
beta0 = 1;
beta1 = 1;
tau = 1;
function bc_model(du, u, h, p, t)
    du[1] = (v0 / (1 + beta0 * (h(p, t - tau)[3]^2))) * (p0 - q0) * u[1] - d0 * u[1]
    du[2] = (v0 / (1 + beta0 * (h(p, t - tau)[3]^2))) * (1 - p0 + q0) * u[1] +
            (v1 / (1 + beta1 * (h(p, t - tau)[3]^2))) * (p1 - q1) * u[2] - d1 * u[2]
    du[3] = (v1 / (1 + beta1 * (h(p, t - tau)[3]^2))) * (1 - p1 + q1) * u[2] - d2 * u[3]
end
lags = [tau]
h(p, t) = ones(3)
h2(p, t) = ones(3) .- t * q1 * 10
tspan = (0.0, 10.0)
u0 = [1.0, 1.0, 1.0]
prob = DDEProblem(bc_model, u0, h, tspan, constant_lags = lags)
alg = MethodOfSteps(Vern9())
sol = solve(prob, alg, reltol = 1e-7, abstol = 1e-10)
prob2 = DDEProblem(bc_model, u0, h2, tspan, constant_lags = lags)
sol2 = solve(prob2, alg, reltol = 1e-7, abstol = 1e-10)

@parameters p0=0.2 p1=0.2 q0=0.3 q1=0.3 v0=1 v1=1 d0=5 d1=1 d2=1 beta0=1 beta1=1
@variables x₀(t) x₁(t) x₂(..)
tau = 1
eqs = [D(x₀) ~ (v0 / (1 + beta0 * (x₂(t - tau)^2))) * (p0 - q0) * x₀ - d0 * x₀
       D(x₁) ~ (v0 / (1 + beta0 * (x₂(t - tau)^2))) * (1 - p0 + q0) * x₀ +
               (v1 / (1 + beta1 * (x₂(t - tau)^2))) * (p1 - q1) * x₁ - d1 * x₁
       D(x₂(t)) ~ (v1 / (1 + beta1 * (x₂(t - tau)^2))) * (1 - p1 + q1) * x₁ - d2 * x₂(t)]
@mtkbuild sys = System(eqs, t)
@test ModelingToolkit.is_dde(sys)
@test !is_markovian(sys)
prob = DDEProblem(sys,
    [x₀ => 1.0, x₁ => 1.0, x₂(t) => 1.0],
    tspan,
    constant_lags = [tau])
sol_mtk = solve(prob, alg, reltol = 1e-7, abstol = 1e-10)
@test sol_mtk.u[end] ≈ sol.u[end]
prob2 = DDEProblem(sys,
    [x₀ => 1.0 - t * q1 * 10, x₁ => 1.0 - t * q1 * 10, x₂(t) => 1.0 - t * q1 * 10],
    tspan,
    constant_lags = [tau])
sol2_mtk = solve(prob2, alg, reltol = 1e-7, abstol = 1e-10)
@test sol2_mtk.u[end] ≈ sol2.u[end]
@test_nowarn sol2_mtk[[x₀, x₁, x₂(t)]]
@test_nowarn sol2_mtk[[x₀, x₁, x₂(t - 0.1)]]

using StochasticDelayDiffEq
function hayes_modelf(du, u, h, p, t)
    τ, a, b, c, α, β, γ = p
    du .= a .* u .+ b .* h(p, t - τ) .+ c
end
function hayes_modelg(du, u, h, p, t)
    τ, a, b, c, α, β, γ = p
    du .= α .* u .+ γ
end
h(p, t) = (ones(1) .+ t);
tspan = (0.0, 10.0)

pmul = [1.0,
    -4.0, -2.0, 10.0,
    -1.3, -1.2, 1.1]

prob = SDDEProblem(hayes_modelf, hayes_modelg, [1.0], h, tspan, pmul;
    constant_lags = (pmul[1],));
sol = solve(prob, RKMil(), seed = 100)

@variables x(..) delx(t)
@parameters a=-4.0 b=-2.0 c=10.0 α=-1.3 β=-1.2 γ=1.1
@brownian η
τ = 1.0
eqs = [D(x(t)) ~ a * x(t) + b * x(t - τ) + c + (α * x(t) + γ) * η, delx ~ x(t - τ)]
@mtkbuild sys = System(eqs, t)
@test ModelingToolkit.has_observed_with_lhs(sys, delx)
@test ModelingToolkit.is_dde(sys)
@test !is_markovian(sys)
@test equations(sys) == [D(x(t)) ~ a * x(t) + b * x(t - τ) + c]
@test isequal(ModelingToolkit.get_noiseeqs(sys), [α * x(t) + γ])
prob_mtk = SDDEProblem(sys, [x(t) => 1.0 + t], tspan; constant_lags = (τ,));
@test_nowarn sol_mtk = solve(prob_mtk, RKMil(), seed = 100)

prob_sa = SDDEProblem(
    sys, [x(t) => 1.0 + t], tspan; constant_lags = (τ,), u0_constructor = SVector{1})
@test prob_sa.u0 isa SVector{1, Float64}

@parameters x(..) a

function oscillator(; name, k = 1.0, τ = 0.01)
    @parameters k=k τ=τ
    @variables x(..)=0.1 y(t)=0.1 jcn(t)=0.0 delx(t)
    eqs = [D(x(t)) ~ y,
        D(y) ~ -k * x(t - τ) + jcn,
        delx ~ x(t - τ)]
    return System(eqs, t; name = name)
end

systems = @named begin
    osc1 = oscillator(k = 1.0, τ = 0.01)
    osc2 = oscillator(k = 2.0, τ = 0.04)
end
eqs = [osc1.jcn ~ osc2.delx,
    osc2.jcn ~ osc1.delx]
@named coupledOsc = System(eqs, t)
@named coupledOsc = compose(coupledOsc, systems)
@test ModelingToolkit.is_dde(coupledOsc)
@test !is_markovian(coupledOsc)
@named coupledOsc2 = System(eqs, t; systems)
@test ModelingToolkit.is_dde(coupledOsc2)
@test !is_markovian(coupledOsc2)
for coupledOsc in [coupledOsc, coupledOsc2]
    local sys = structural_simplify(coupledOsc)
    @test length(equations(sys)) == 4
    @test length(unknowns(sys)) == 4
end
sys = structural_simplify(coupledOsc)
prob = DDEProblem(sys, [], (0.0, 10.0); constant_lags = [sys.osc1.τ, sys.osc2.τ])
sol = solve(prob, MethodOfSteps(Tsit5()))
obsfn = ModelingToolkit.build_explicit_observed_function(
    sys, [sys.osc1.delx, sys.osc2.delx])
@test_nowarn sol[[sys.osc1.delx, sys.osc2.delx]]
@test sol[sys.osc1.delx] ≈ sol(sol.t .- 0.01; idxs = sys.osc1.x).u

prob_sa = DDEProblem(sys, [], (0.0, 10.0); constant_lags = [sys.osc1.τ, sys.osc2.τ],
    u0_constructor = SVector{4})
@test prob_sa.u0 isa SVector{4, Float64}

@testset "DDE observed with array variables" begin
    @component function valve(; name)
        @parameters begin
            open(t)::Bool = false
            Kp = 2
            Ksnap = 1.1
            τ = 0.1
        end
        @variables begin
            opening(..)
            lag_opening(t)
            snap_opening(t)
        end
        eqs = [D(opening(t)) ~ Kp * (open - opening(t))
               lag_opening ~ opening(t - τ)
               snap_opening ~ clamp(Ksnap * lag_opening - 1 / Ksnap, 0, 1)]
        return System(eqs, t; name = name)
    end

    @component function veccy(; name)
        @parameters dx[1:3] = ones(3)
        @variables begin
            x(t)[1:3] = zeros(3)
        end
        return System([D(x) ~ dx], t; name = name)
    end

    @mtkbuild ssys = System(
        Equation[], t; systems = [valve(name = :valve), veccy(name = :vvecs)])
    prob = DDEProblem(ssys, [ssys.valve.opening => 1.0], (0.0, 1.0))
    sol = solve(prob, MethodOfSteps(Tsit5()))
    obsval = @test_nowarn sol[ssys.valve.lag_opening + sum(ssys.vvecs.x)]
    @test obsval ≈
          sol(sol.t .- prob.ps[ssys.valve.τ]; idxs = ssys.valve.opening).u .+
          sum.(sol[ssys.vvecs.x])
end

@testset "Issue#3165 DDEs with non-tunables" begin
    @variables x(..) = 1.0
    @parameters w=1.0 [tunable = false] τ=0.5
    eqs = [D(x(t)) ~ -w * x(t - τ)]

    @named sys = System(eqs, t)
    sys = structural_simplify(sys)

    prob = DDEProblem(sys,
        [],
        (0.0, 10.0),
        constant_lags = [τ])

    alg = MethodOfSteps(Vern7())
    @test_nowarn solve(prob, alg)

    @brownian r
    eqs = [D(x(t)) ~ -w * x(t - τ) + r]
    @named sys = System(eqs, t)
    sys = structural_simplify(sys)
    prob = SDDEProblem(sys,
        [],
        (0.0, 10.0),
        constant_lags = [τ])

    @test_nowarn solve(prob, RKMil())
end
