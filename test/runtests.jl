using SciCompDSL
using Base.Test

@testset "Parsing Test" begin
    @DVar x y=sin(1)+exp(1) z
    x1 = DependentVariable(:x)
    y1 = DependentVariable(:y, sin(1) + exp(1))
    z1 = DependentVariable(:z)
    @test isequal(x1, x)
    @test isequal(y1, y)
    @test isequal(z1, z)
    @IVar begin
        t
        s = cos(2.5)
    end
    t1 = IndependentVariable(:t)
    s1 = IndependentVariable(:s, cos(2.5))
    @test isequal(t1, t)
    @test isequal(s1, s)
    @Deriv D''~t
    D1 = Differential(t, 2)
    @test isequal(D1, D)
    @Const c=0 v=2
    c1 = Constant(0)
    v1 = Constant(2)
    @test isequal(c1, c)
    @test isequal(v1, v)
end

# Define some variables
@DVar x y z
@IVar t
@Deriv D'~t # Default of first derivative, Derivative(t,1)
@Param σ ρ β
@Const c=0
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
@testset "Derivatives Test" begin
    dsin = D*sin(t)
    @test isequal(expand_derivatives(dsin),cos(t))
    dcsch = D*csch(t)
    @test isequal(expand_derivatives(dcsch),Operation(-coth(t)*csch(t)))
    # Chain rule
    dsinsin = D*sin(sin(t))
    @test isequal(expand_derivatives(dsinsin),cos(sin(t))*cos(t))
end
