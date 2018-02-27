using SciCompDSL
using Base.Test

# Define some variables
x = DependentVariable(:x)
y = DependentVariable(:y)
z = DependentVariable(:z)
t = IndependentVariable(:t)
D = Differential(t) # Default of first derivative, Derivative(t,1)
σ = Parameter(:σ)
ρ = Parameter(:ρ)
β = Parameter(:β)

σ*(y-x)
D*x
D*x == σ*(y-x)
D*y == x*(ρ-z)-y

# Define a differential equation
eqs = [D*x == σ*(y-x),
       D*y == x*(ρ-z)-y,
       D*z == x*y - β*z]
de = DiffEqSystem(eqs,IndependentVariable[],DependentVariable[],Variable[],Parameter[])

# Define a nonlinear system
eqs = [Constant(0) == σ*(y-x),
       Constant(0) == x*(ρ-z)-y,
       Constant(0) == x*y - β*z]
ns = NonlinearSystem(eqs,Variable[],Parameter[])
