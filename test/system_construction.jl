using ModelingToolkit
using Test

# Define some variables
@Param t σ ρ β
@Unknown x(t) y(t) z(t)
@Deriv D'~t # Default of first derivative, Derivative(t,1)
@Const c=0

# Define a differential equation
eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]
de = DiffEqSystem(eqs,[t],[x,y,z],[σ,ρ,β])
ModelingToolkit.generate_ode_function(de)
ModelingToolkit.generate_ode_function(de;version=ModelingToolkit.SArrayFunction)
jac_expr = ModelingToolkit.generate_ode_jacobian(de)
jac = ModelingToolkit.calculate_jacobian(de)
f = ODEFunction(de)
ModelingToolkit.generate_ode_iW(de)

# Differential equation with automatic extraction of variables
de2 = DiffEqSystem(eqs, [t])

function test_vars_extraction(de, de2)
    for el in (:ivs, :dvs, :ps)
        names2 = sort(collect(var.name for var in getfield(de2,el)))
        names = sort(collect(var.name for var in getfield(de,el)))
        @test names2 == names
    end
end
test_vars_extraction(de, de2)

# Conversion to first-order ODEs #17
@Deriv D3'''~t
@Deriv D2''~t
@Unknown u(t) u_tt(t) u_t(t) x_t(t)
eqs = [D3(u) ~ 2(D2(u)) + D(u) + D(x) + 1
       D2(x) ~ D(x) + 2]
de = DiffEqSystem(eqs, [t])
de1 = ode_order_lowering(de)
lowered_eqs = [D(u_tt) ~ 2u_tt + u_t + x_t + 1
               D(x_t)  ~ x_t + 2
               D(u_t)  ~ u_tt
               D(u)    ~ u_t
               D(x)    ~ x_t]
function test_eqs(eqs1, eqs2)
    length(eqs1) == length(eqs2) || return false
    eq = true
    for (eq1, eq2) in zip(eqs1, eqs2)
        lhs1, lhs2 = eq1.lhs, eq2.lhs
        typeof(lhs1) === typeof(lhs2) || return false
        for f in fieldnames(typeof(lhs1))
            eq = eq & isequal(getfield(lhs1, f), getfield(lhs2, f))
        end
        eq = eq & isequal(eq1.rhs, eq2.rhs)
    end
    eq
end
@test test_eqs(de1.eqs, lowered_eqs)

# Internal calculations
a = y - x
eqs = [D(x) ~ σ*a,
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]
de = DiffEqSystem(eqs,[t],[x,y,z],[σ,ρ,β])
ModelingToolkit.generate_ode_function(de)
jac = ModelingToolkit.calculate_jacobian(de)
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

ModelingToolkit.generate_nlsys_function(ns)

@Deriv D'~t
@Param A B C
_x = y / C
eqs = [D(x) ~ -A*x,
       D(y) ~ A*x - B*_x]
de = DiffEqSystem(eqs,[t],[x,y],[A,B,C])
test_vars_extraction(de, DiffEqSystem(eqs,[t]))
test_vars_extraction(de, DiffEqSystem(eqs))
@test eval(ModelingToolkit.generate_ode_function(de))([0.0,0.0],[1.0,2.0],[1,2,3],0.0) ≈ -1/3

# Now nonlinear system with only variables
@Unknown x y z
@Param σ ρ β

# Define a nonlinear system
eqs = [0 ~ σ*(y-x),
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs)
nlsys_func = ModelingToolkit.generate_nlsys_function(ns)
jac = ModelingToolkit.generate_nlsys_jacobian(ns)
f = @eval eval(nlsys_func)

# Intermediate calculations
# Define a nonlinear system
eqs = [a ~ y-x,
       0 ~ σ*a,
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs,[x,y,z],[σ,ρ,β])
nlsys_func = ModelingToolkit.generate_nlsys_function(ns)
jac = ModelingToolkit.calculate_jacobian(ns)
jac = ModelingToolkit.generate_nlsys_jacobian(ns)
