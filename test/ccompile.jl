using ModelingToolkit, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters a
@variables x y
eqs = [D(x) ~ a * x - x * y,
    D(y) ~ -3y + x * y]
f = build_function([x.rhs for x in eqs], [x, y], [a], t, expression = Val{false},
    target = ModelingToolkit.CTarget())
f2 = eval(build_function([x.rhs for x in eqs], [x, y], [a], t)[2])
du = rand(2);
du2 = rand(2);
u = rand(2)
p = rand(1)
_t = rand()
f(du, u, p, _t)
f2(du2, u, p, _t)
@test du ≈ du2
