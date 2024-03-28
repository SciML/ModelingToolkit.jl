using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D, get_u0
using SymbolicIndexingInterface: getu

@variables x(t)[1:3]=[1.0, 2.0, 3.0] y(t) z(t)[1:2]

@mtkbuild sys=ODESystem([D(x) ~ t * x], t) simplify=false
@test get_u0(sys, [])[1] == [1.0, 2.0, 3.0]
@test get_u0(sys, [x => [2.0, 3.0, 4.0]])[1] == [2.0, 3.0, 4.0]
@test get_u0(sys, [x[1] => 2.0, x[2] => 3.0, x[3] => 4.0])[1] == [2.0, 3.0, 4.0]
@test get_u0(sys, [2.0, 3.0, 4.0])[1] == [2.0, 3.0, 4.0]

@mtkbuild sys=ODESystem([
        D(x) ~ 3x,
        D(y) ~ t,
        D(z[1]) ~ z[2] + t,
        D(z[2]) ~ y + z[1]
    ], t) simplify=false

@test_throws ModelingToolkit.MissingVariablesError get_u0(sys, [])
getter = getu(sys, [x..., y, z...])
@test getter(get_u0(sys, [y => 4.0, z => [5.0, 6.0]])[1]) == collect(1.0:6.0)
@test getter(get_u0(sys, [y => 4.0, z => [3y, 4y]])[1]) == [1.0, 2.0, 3.0, 4.0, 12.0, 16.0]
@test getter(get_u0(sys, [y => 3.0, z[1] => 3y, z[2] => 2x[1]])[1]) ==
      [1.0, 2.0, 3.0, 3.0, 9.0, 2.0]

@variables w(t)
@parameters p1 p2

@test getter(get_u0(sys, [y => 2p1, z => [3y, 2p2]], [p1 => 5.0, p2 => 6.0])[1]) ==
      [1.0, 2.0, 3.0, 10.0, 30.0, 12.0]
@test_throws Any getter(get_u0(sys, [y => 2w, w => 3.0, z[1] => 2p1, z[2] => 3p2]))
@test getter(get_u0(
    sys, [y => 2w, w => 3.0, z[1] => 2p1, z[2] => 3p2], [p1 => 3.0, p2 => 4.0])[1]) ==
      [1.0, 2.0, 3.0, 6.0, 6.0, 12.0]

# Issue#2566
@variables X(t)
@parameters p1 p2 p3

p_vals = [p1 => 1.0, p2 => 2.0]
u_vals = [X => 3.0]

var_vals = [p1 => 1.0, p2 => 2.0, X => 3.0]
desired_values = [p1, p2, p3]
defaults = Dict([p3 => X])
vals = ModelingToolkit.varmap_to_vars(var_vals, desired_values; defaults = defaults)
@test vals == [1.0, 2.0, 3.0]

# Issue#2565
# Create ODESystem.
@variables X1(t) X2(t)
@parameters k1 k2 Γ[1:1]=X1 + X2
eq = D(X1) ~ -k1 * X1 + k2 * (-X1 + Γ[1])
obs = X2 ~ Γ[1] - X1
@mtkbuild osys_m = ODESystem([eq], t, [X1], [k1, k2, Γ[1]]; observed = [X2 ~ Γ[1] - X1])

# Creates ODEProblem.
u0 = [X1 => 1.0, X2 => 2.0]
tspan = (0.0, 1.0)
ps = [k1 => 1.0, k2 => 5.0]
@test_nowarn oprob = ODEProblem(osys_m, u0, tspan, ps)

# Make sure it doesn't error on array variables with unspecified size
@parameters p::Vector{Real} q[1:3]
varmap = Dict(p => ones(3), q => 2ones(3))
cvarmap = ModelingToolkit.canonicalize_varmap(varmap)
target_varmap = Dict(p => ones(3), q => 2ones(3), q[1] => 2.0, q[2] => 2.0, q[3] => 2.0)
@test cvarmap == target_varmap
