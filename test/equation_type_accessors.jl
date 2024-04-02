# Fetch packages.
using ModelingToolkit
import ModelingToolkit: get_systems, namespace_equations
import ModelingToolkit: is_alg_equation, is_diff_equation
import ModelingToolkit: t_nounits as t, D_nounits as D, wrap, get_eqs

# Creates equations.
@variables X(t) Y(t) Z(t)
@parameters a b c d
eq1 = X^Z - Z^(X + 1) ~ log(X - a + b) * Y
eq2 = X ~ Y^(X + 1)
eq3 = a + b + c + d ~ X * (Y + d * (Y + Z))
eq4 = X ~ sqrt(a + Z) + t
eq5 = D(D(X)) ~ a^(2Y) + 3Z * t - 6
eq6 = X * (Z - Z * (b + X)) ~ c^(X + D(Y))
eq7 = sqrt(X + c) ~ 2 * (Y + log(a + D(Z)))
eq8 = -0.1 ~ D(Z) + X

@test is_alg_equation(eq1)
@test is_alg_equation(eq2)
@test is_alg_equation(eq3)
@test is_alg_equation(eq4)
@test !is_alg_equation(eq5)
@test !is_alg_equation(eq6)
@test !is_alg_equation(eq7)
@test !is_alg_equation(eq8)

@test !is_diff_equation(eq1)
@test !is_diff_equation(eq2)
@test !is_diff_equation(eq3)
@test !is_diff_equation(eq4)
@test is_diff_equation(eq5)
@test is_diff_equation(eq6)
@test is_diff_equation(eq7)
@test is_diff_equation(eq8)

# Creates systems.
eqs1 = [X * Y + a ~ Z^3 - X * log(b + Y)
        X ~ Z * Y * X + a + b
        c * sin(X) + sin(Y) ~ d * (a + X * (b + Y * (c + Z)))]
eqs2 = [X + Y + c ~ b * X^(X + Z + a)
        D(X) ~ a * Y + b * X + c * Z
        D(Z) + Z * Y ~ X - log(Z)]
eqs3 = [D(X) ~ sqrt(X + b) + sqrt(Z + c)
        2Z * (Z + Y) ~ D(Y) * log(a)
        D(Z) + c * X ~ b / (X + Y^d) + D(Z)]
@named osys1 = ODESystem(eqs1, t)
@named osys2 = ODESystem(eqs2, t)
@named osys3 = ODESystem(eqs3, t)

# Test `has...` for non-composed systems. 
@test has_alg_equations(osys1)
@test has_alg_equations(osys2)
@test !has_alg_equations(osys3)
@test has_alg_eqs(osys1)
@test has_alg_eqs(osys2)
@test !has_alg_eqs(osys3)
@test !has_diff_equations(osys1)
@test has_diff_equations(osys2)
@test has_diff_equations(osys3)
@test !has_diff_eqs(osys1)
@test has_diff_eqs(osys2)
@test has_diff_eqs(osys3)

# Test getters for non-composed systems.
isequal(alg_equations(osys1), eqs1)
isequal(alg_equations(osys2), eqs2[1:1])
isequal(alg_equations(osys3), [])
isequal(get_alg_eqs(osys1), eqs1)
isequal(get_alg_eqs(osys2), eqs2[1:1])
isequal(get_alg_eqs(osys3), [])
isequal(diff_equations(osys1), [])
isequal(diff_equations(osys2), eqs2[2:3])
isequal(diff_equations(osys3), eqs3)
isequal(get_diff_eqs(osys1), [])
isequal(get_diff_eqs(osys2), eqs2[2:3])
isequal(get_diff_eqs(osys3), eqs3)

# Creates composed systems.
osys1_1 = compose(osys1, [osys1])
osys1_12 = compose(osys1, [osys1, osys2])
osys1_12_1 = compose(osys1, [osys1, compose(osys2, [osys1])])
osys3_2 = compose(osys3, [osys2])
osys3_33 = compose(osys3, [osys3, osys3])

# Test `has...` for composed systems.
@test has_alg_equations(osys1_1)
@test !has_diff_equations(osys1_1)
@test has_alg_eqs(osys1_1)
@test !has_diff_eqs(osys1_1)
@test has_alg_equations(get_systems(osys1_1)[1])
@test !has_diff_equations(get_systems(osys1_1)[1])
@test has_alg_eqs(get_systems(osys1_1)[1])
@test !has_diff_eqs(get_systems(osys1_1)[1])

@test has_alg_equations(osys1_12)
@test has_diff_equations(osys1_12)
@test has_alg_eqs(osys1_12)
@test !has_diff_eqs(osys1_12)
@test has_alg_equations(get_systems(osys1_12)[1])
@test !has_diff_equations(get_systems(osys1_12)[1])
@test has_alg_eqs(get_systems(osys1_12)[1])
@test !has_diff_eqs(get_systems(osys1_12)[1])
@test has_alg_equations(get_systems(osys1_12)[2])
@test has_diff_equations(get_systems(osys1_12)[2])
@test has_alg_eqs(get_systems(osys1_12)[2])
@test has_diff_eqs(get_systems(osys1_12)[2])

@test has_alg_equations(osys1_12_1)
@test has_diff_equations(osys1_12_1)
@test has_alg_eqs(osys1_12_1)
@test !has_diff_eqs(osys1_12_1)
@test has_alg_equations(get_systems(osys1_12_1)[1])
@test !has_diff_equations(get_systems(osys1_12_1)[1])
@test has_alg_eqs(get_systems(osys1_12_1)[1])
@test !has_diff_eqs(get_systems(osys1_12_1)[1])
@test has_alg_equations(get_systems(osys1_12_1)[2])
@test has_diff_equations(get_systems(osys1_12_1)[2])
@test has_alg_eqs(get_systems(osys1_12_1)[2])
@test has_diff_eqs(get_systems(osys1_12_1)[2])
@test has_alg_equations(get_systems(get_systems(osys1_12_1)[2])[1])
@test !has_diff_equations(get_systems(get_systems(osys1_12_1)[2])[1])
@test has_alg_eqs(get_systems(get_systems(osys1_12_1)[2])[1])
@test !has_diff_eqs(get_systems(get_systems(osys1_12_1)[2])[1])

@test has_alg_equations(osys3_2)
@test has_diff_equations(osys3_2)
@test !has_alg_eqs(osys3_2)
@test has_diff_eqs(osys3_2)
@test has_alg_equations(get_systems(osys3_2)[1])
@test has_diff_equations(get_systems(osys3_2)[1])
@test has_alg_eqs(get_systems(osys3_2)[1])
@test has_diff_eqs(get_systems(osys3_2)[1])

@test !has_alg_equations(osys3_33)
@test has_diff_equations(osys3_33)
@test !has_alg_eqs(osys3_33)
@test has_diff_eqs(osys3_33)
@test !has_alg_equations(get_systems(osys3_33)[1])
@test has_diff_equations(get_systems(osys3_33)[1])
@test !has_alg_eqs(get_systems(osys3_33)[1])
@test has_diff_eqs(get_systems(osys3_33)[1])
@test !has_alg_equations(get_systems(osys3_33)[2])
@test has_diff_equations(get_systems(osys3_33)[2])
@test !has_alg_eqs(get_systems(osys3_33)[2])
@test has_diff_eqs(get_systems(osys3_33)[2])

# Test getters for composed systems.
ns_eqs1 = namespace_equations(osys1)
ns_eqs2 = namespace_equations(osys2)
ns_eqs3 = namespace_equations(osys3)

isequal(alg_equations(osys1_1), vcat(eqs1, ns_eqs1))
isequal(diff_equations(osys1_1), [])
isequal(get_alg_eqs(osys1_1), eqs1)
isequal(get_diff_eqs(osys1_1), [])
isequal(alg_equations(get_systems(osys1_1)[1]), eqs1)
isequal(diff_equations(get_systems(osys1_1)[1]), [])
isequal(get_alg_eqs(get_systems(osys1_1)[1]), eqs1)
isequal(get_diff_eqs(get_systems(osys1_1)[1]), [])

isequal(alg_equations(osys1_12), vcat(eqs1, ns_eqs1, filter(is_alg_equation, ns_eqs2)))
isequal(diff_equations(osys1_12), filter(is_diff_equation, ns_eqs2))
isequal(get_alg_eqs(osys1_12), eqs1)
isequal(get_diff_eqs(osys1_12), [])
isequal(alg_equations(get_systems(osys1_12)[1]), eqs1)
isequal(diff_equations(get_systems(osys1_12)[1]), [])
isequal(get_alg_eqs(get_systems(osys1_12)[1]), eqs1)
isequal(get_diff_eqs(get_systems(osys1_12)[1]), [])
isequal(alg_equations(get_systems(osys1_12)[2]), eqs2[1:1])
isequal(diff_equations(get_systems(osys1_12)[2]), eqs2[2:3])
isequal(get_alg_eqs(get_systems(osys1_12)[2]), eqs2[1:1])
isequal(get_diff_eqs(get_systems(osys1_12)[2]), eqs2[2:3])

isequal(alg_equations(osys3_2), vcat(filter(is_alg_equation, ns_eqs2)))
isequal(diff_equations(osys3_2), vcat(eqs3, filter(is_diff_equation, ns_eqs2)))
isequal(get_alg_eqs(osys3_2), [])
isequal(get_diff_eqs(osys3_2), eqs3)
isequal(alg_equations(get_systems(osys3_2)[1]), eqs2[1:1])
isequal(diff_equations(get_systems(osys3_2)[1]), eqs2[2:3])
isequal(get_alg_eqs(get_systems(osys3_2)[1]), eqs2[1:1])
isequal(get_diff_eqs(get_systems(osys3_2)[1]), eqs2[2:3])
