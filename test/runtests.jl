using SciCompDSL
using Base.Test

# Define some variables
@dvars x y z
@idvars t
D = Differential(t) # Default of first derivative, Derivative(t,1)
@paras σ ρ β
c = Constant(0)
σ*(y-x)
D*x
D*x == -σ*(y-x)
D*y == x*(ρ-z)-sin(y)
s = JumpVariable(:s,3)
n = NoiseVariable(:n)

@test isequal(D*t,Constant(1))
null_op = 0*t
@test isequal(simplify_constants(null_op),Constant(0))

# Define a differential equation
eqs = [D*x == σ*(y-x),
       D*y == x*(ρ-z)-y,
       D*z == x*y - β*z]
de = DiffEqSystem(eqs,[t],[x,y,z],Variable[],[σ,ρ,β])
map(eq->:($(eq.args[1].name) = $(eq.args[2])),de.eqs)
SciCompDSL.generate_ode_function(de)
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


# Default values
p = Parameter(:p, 1)
u = DependentVariable(:u, [1])

# Derivatives
dsin = D*sin(t)
isequal(expand_derivatives(dsin),cos(t))
