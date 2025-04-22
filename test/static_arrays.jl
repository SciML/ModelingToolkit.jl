using ModelingToolkit, SciMLBase, StaticArrays, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters σ ρ β
@variables x(t) y(t) z(t)

eqs = [D(D(x)) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

@named sys = ODESystem(eqs, t)
sys = structural_simplify(sys)

u0 = @SVector [D(x) => 2.0,
    x => 1.0,
    y => 0.0,
    z => 0.0]

p = @SVector [σ => 28.0,
    ρ => 10.0,
    β => 8 / 3]

tspan = (0.0, 100.0)
prob_mtk = ODEProblem(sys, u0, tspan, p)

@test !SciMLBase.isinplace(prob_mtk)
@test prob_mtk.u0 isa SArray
