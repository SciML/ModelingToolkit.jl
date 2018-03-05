using SciCompDSL
using Base.Test

# Define some variables
@DVar x y z
@IVar t
@Deriv D'~t # Default of first derivative, Derivative(t,1)
@Param σ ρ β
@Const c=0

# Define a differential equation
eqs = [D*x == σ*(y-x),
       D*y == x*(ρ-z)-y,
       D*z == x*y - β*z]
de = DiffEqSystem(eqs,[t],[x,y,z],Variable[],[σ,ρ,β])
SciCompDSL.generate_ode_function(de)
jac = SciCompDSL.generate_ode_jacobian(de,false)
jac = SciCompDSL.generate_ode_jacobian(de)
f = DiffEqFunction(de)

# Differential equation with automatic extraction of variables on rhs
de2 = DiffEqSystem(eqs, [t])
for el in (:ivs, :dvs, :vs, :ps)
    names2 = sort(collect(var.name for var in getfield(de2,el)))
    names = sort(collect(var.name for var in getfield(de,el)))
    @test names2 == names
end

# Define a nonlinear system
eqs = [0 == σ*(y-x),
       0 == x*(ρ-z)-y,
       0 == x*y - β*z]
ns = NonlinearSystem(eqs,[x,y,z],[σ,ρ,β])
ns2 = NonlinearSystem(eqs)
for el in (:vs, :ps)
    names2 = sort(collect(var.name for var in getfield(ns2,el)))
    names = sort(collect(var.name for var in getfield(ns,el)))
    @test names2 == names
end

SciCompDSL.generate_nlsys_function(ns)

# Now nonlinear system with only variables
@Var x y z
@Param σ ρ β

# Define a nonlinear system
eqs = [0 == σ*(y-x),
       0 == x*(ρ-z)-y,
       0 == x*y - β*z]
ns = NonlinearSystem(eqs)
nlsys_func = SciCompDSL.generate_nlsys_function(ns)
jac = SciCompDSL.generate_nlsys_jacobian(ns,false)
jac = SciCompDSL.generate_nlsys_jacobian(ns)
f = @eval eval(nlsys_func)
