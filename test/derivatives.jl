using ModelingToolkit
using Test

# Derivatives
@parameters t σ ρ β
@variables x y z
@derivatives D'~t D2''~t Dx'~x

@test @macroexpand(@derivatives D'~t D2''~t) == @macroexpand(@derivatives (D'~t), (D2''~t))

@test isequal(expand_derivatives(D(t)), 1)
@test isequal(expand_derivatives(D(D(t))), 0)

dsin = D(sin(t))
@test isequal(expand_derivatives(dsin), cos(t))

dcsch = D(csch(t))
@test isequal(expand_derivatives(dcsch), simplify_constants(coth(t) * csch(t) * -1))

@test isequal(expand_derivatives(D(-7)), 0)
@test isequal(expand_derivatives(D(sin(2t))), simplify_constants(cos(2t) * 2))
@test isequal(expand_derivatives(D2(sin(t))), simplify_constants(-sin(t)))
@test isequal(expand_derivatives(D2(sin(2t))), simplify_constants(sin(2t) * -4))
@test isequal(expand_derivatives(D2(t)), 0)
@test isequal(expand_derivatives(D2(5)), 0)

# Chain rule
dsinsin = D(sin(sin(t)))
@test isequal(expand_derivatives(dsinsin), cos(sin(t))*cos(t))

d1 = D(sin(t)*t)
d2 = D(sin(t)*cos(t))
@test isequal(expand_derivatives(d1), t*cos(t)+sin(t))
@test isequal(expand_derivatives(d2), simplify_constants(cos(t)*cos(t)+(sin(t)*-1)*sin(t)))

eqs = [0 ~ σ*(y-x),
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
sys = NonlinearSystem(eqs, [x,y,z], [σ,ρ,β])
jac = calculate_jacobian(sys)
@test isequal(jac[1,1], σ*-1)
@test isequal(jac[1,2], σ)
@test isequal(jac[1,3], 0)
@test isequal(jac[2,1], ρ-z)
@test isequal(jac[2,2], -1)
@test isequal(jac[2,3], x*-1)
@test isequal(jac[3,1], y)
@test isequal(jac[3,2], x)
@test isequal(jac[3,3], -1*β)

# Variable dependence checking in differentiation
@variables a(t) b(a)
@test !isequal(D(b), 0)
@test isequal(expand_derivatives(D(t)), 1)
@test isequal(expand_derivatives(Dx(x)), 1)

@variables x(t) y(t) z(t)

@test isequal(expand_derivatives(D(x * y)), simplify_constants(y*D(x) + x*D(y)))
@test_broken isequal(expand_derivatives(D(x * y)), simplify_constants(D(x)*y + x*D(y)))

@test isequal(expand_derivatives(D(2t)), 2)
@test isequal(expand_derivatives(D(2x)), 2D(x))
@test isequal(expand_derivatives(D(x^2)), simplify_constants(2 * x * D(x)))

# n-ary * and +
isequal(ModelingToolkit.derivative(Operation(*, [x, y, z*ρ]), 1), y*(z*ρ))
isequal(ModelingToolkit.derivative(Operation(+, [x*y, y, z]), 1), 1)

@test iszero(ModelingToolkit.derivative(ModelingToolkit.Constant(42), x))
@test all(iszero, ModelingToolkit.gradient(ModelingToolkit.Constant(42), [t, x, y, z]))
@test all(iszero, ModelingToolkit.hessian(ModelingToolkit.Constant(42), [t, x, y, z]))
@test isequal(ModelingToolkit.jacobian([t, x, ModelingToolkit.Constant(42)], [t, x]),
              Expression[ModelingToolkit.Constant(1)  ModelingToolkit.Constant(0)
                         Differential(t)(x)           ModelingToolkit.Constant(1)
                         ModelingToolkit.Constant(0)  ModelingToolkit.Constant(0)])

# issue 252
@variables beta, alpha, delta
@variables x1, x2, x3

# expression
tmp = beta * (alpha * exp(x1) * x2 ^ (alpha - 1) + 1 - delta) / x3
# derivative w.r.t. x1 and x2
t1 = ModelingToolkit.gradient(tmp, [x1, x2])
@test_nowarn ModelingToolkit.gradient(t1[1], [beta])
