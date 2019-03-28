using ModelingToolkit
using Test

# Define some variables
@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

# Define a differential equation
eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]
de = DiffEqSystem(eqs,t,[x,y,z],[σ,ρ,β])
generate_function(de)
generate_function(de;version=ModelingToolkit.SArrayFunction)
jac_expr = generate_jacobian(de)
jac = calculate_jacobian(de)
f = ODEFunction(de)
ModelingToolkit.generate_ode_iW(de)

# Differential equation with automatic extraction of variables
de2 = DiffEqSystem(eqs, t)

function test_vars_extraction(de, de2)
    @test isequal(de.iv, de2.iv)
    for el in (:dvs, :ps)
        names2 = sort(collect(var.name for var in getfield(de2,el)))
        names = sort(collect(var.name for var in getfield(de,el)))
        @test names2 == names
    end
end
test_vars_extraction(de, de2)

@testset "time-varying parameters" begin
    @parameters σ′(t)
    eqs = [D(x) ~ σ′*(y-x),
           D(y) ~ x*(ρ-z)-y,
           D(z) ~ x*y - β*z]
    de = DiffEqSystem(eqs,t,[x,y,z],[σ′,ρ,β])
    @test begin
        f = eval(generate_function(de))
        du = [0.0,0.0,0.0]
        f(du, [1.0,2.0,3.0], [x->x+7,2,3], 5.0)
        du ≈ [12, -3, -7]
    end
end

# Conversion to first-order ODEs #17
@derivatives D3'''~t
@derivatives D2''~t
@variables u(t) u_tt(t) u_t(t) x_t(t)
eqs = [D3(u) ~ 2(D2(u)) + D(u) + D(x) + 1
       D2(x) ~ D(x) + 2]
de = DiffEqSystem(eqs, t)
de1 = ode_order_lowering(de)
lowered_eqs = [D(u_tt) ~ 2u_tt + u_t + x_t + 1
               D(x_t)  ~ x_t + 2
               D(u_t)  ~ u_tt
               D(u)    ~ u_t
               D(x)    ~ x_t]
@test de1.eqs == convert.(ModelingToolkit.DiffEq, lowered_eqs)
@test isequal(de1.dvs, [u_tt, x_t, u_t, u, x])
@test isequal(de1.iv, t)
du = zeros(5)
ODEFunction(de1)(du, ones(5), nothing, 0.1)
@test du == [5.0, 3.0, 1.0, 1.0, 1.0]

# Internal calculations
a = y - x
eqs = [D(x) ~ σ*a,
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]
de = DiffEqSystem(eqs,t,[x,y,z],[σ,ρ,β])
generate_function(de)
jac = calculate_jacobian(de)
f = ODEFunction(de)

# Define a nonlinear system
eqs = [0 ~ σ*(y-x),
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs,[x,y,z],[t,σ,ρ,β])
ns2 = NonlinearSystem(eqs)
for el in (:vs, :ps)
    names2 = sort(collect(var.name for var in getfield(ns2,el)))
    names = sort(collect(var.name for var in getfield(ns,el)))
    @test names2 == names
end

generate_function(ns)

@derivatives D'~t
@parameters A B C
_x = y / C
eqs = [D(x) ~ -A*x,
       D(y) ~ A*x - B*_x]
de = DiffEqSystem(eqs,t,[x,y],[A,B,C])
test_vars_extraction(de, DiffEqSystem(eqs,t))
test_vars_extraction(de, DiffEqSystem(eqs))
@test begin
    f = eval(generate_function(de))
    du = [0.0,0.0]
    f(du, [1.0,2.0], [1,2,3], 0.0)
    du ≈ [-1, -1/3]
end

# Now nonlinear system with only variables
@variables x y z
@parameters σ ρ β

# Define a nonlinear system
eqs = [0 ~ σ*(y-x),
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs, [x,y,z], [σ,ρ,β])
jac = calculate_jacobian(ns)
@testset "nlsys jacobian" begin
    @test isequal(jac[1,1], σ * -1)
    @test isequal(jac[1,2], σ)
    @test isequal(jac[1,3], 0)
    @test isequal(jac[2,1], ρ - z)
    @test isequal(jac[2,2], -1)
    @test isequal(jac[2,3], x * -1)
    @test isequal(jac[3,1], y)
    @test isequal(jac[3,2], x)
    @test isequal(jac[3,3], -1 * β)
end
nlsys_func = generate_function(ns)
jac_func = generate_jacobian(ns)
f = @eval eval(nlsys_func)

# Intermediate calculations
# Define a nonlinear system
eqs = [0 ~ σ*a,
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs,[x,y,z],[σ,ρ,β])
nlsys_func = generate_function(ns)
jac = calculate_jacobian(ns)
jac = generate_jacobian(ns)
