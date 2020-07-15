using ModelingToolkit
using Test

@parameters t
@variables x(t) y(t) z(t)

null_op = 0*t
@test isequal(simplify(null_op), 0)

one_op = 1*t
@test isequal(simplify(one_op), t)

identity_op = Num(Term(identity,[x.val]))
@test isequal(simplify(identity_op), x)

minus_op = -x
@test isequal(simplify(minus_op), -1x)
simplify(minus_op)

@variables x

@test toexpr(expand_derivatives(Differential(x)((x-2)^2))) == :(2 * (-2 + x))
@test toexpr(expand_derivatives(Differential(x)((x-2)^3))) == :(3 * (-2 + x)^2)
@test toexpr(simplify(x+2+3)) == :(5 + x)

d1 = Differential(x)((-2 + x)^2)
d2 = Differential(x)(d1)
d3 = Differential(x)(d2)

@test toexpr(expand_derivatives(d3)) == :(0)
@test toexpr(simplify(x^0)) == :(1)

@test ModelingToolkit.substitute(2x + y == 1, Dict(x => 0.0, y => 0.0)) === false
@test ModelingToolkit.substitute(2x + y == 1, Dict(x => 0.0, y => 1.0)) === true
