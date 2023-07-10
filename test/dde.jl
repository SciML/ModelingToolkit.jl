using ModelingToolkit, DelayDiffEq, Test
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
tspan = (0.0, 10.0)
u0 = [1.0, 1.0, 1.0]
prob = DDEProblem(bc_model, u0, h, tspan, constant_lags = lags)
alg = MethodOfSteps(Tsit5())
sol = solve(prob, alg)

@parameters p0=0.2 p1=0.2 q0=0.3 q1=0.3 v0=1 v1=1 d0=5 d1=1 d2=1 beta0=1 beta1=1 tau=1
@variables t x₀(t) x₁(t) x₂(..)
D = Differential(t)
eqs = [D(x₀) ~ (v0 / (1 + beta0 * (x₂(t - tau)^2))) * (p0 - q0) * x₀ - d0 * x₀
    D(x₁) ~ (v0 / (1 + beta0 * (x₂(t - tau)^2))) * (1 - p0 + q0) * x₀ +
            (v1 / (1 + beta1 * (x₂(t - tau)^2))) * (p1 - q1) * x₁ - d1 * x₁
    D(x₂(t)) ~ (v1 / (1 + beta1 * (x₂(t - tau)^2))) * (1 - p1 + q1) * x₁ - d2 * x₂(t)]
@named sys = System(eqs)
prob = DDEProblem(sys,
    [x₀ => 1.0, x₁ => 1.0, x₂(t) => 1.0],
    h,
    tspan,
    constant_lags = [tau])
sol_mtk = solve(prob, alg)
@test Array(sol_mtk) ≈ Array(sol)
