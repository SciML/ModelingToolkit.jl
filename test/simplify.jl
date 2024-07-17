using ModelingToolkit
using ModelingToolkit: value
using Test

@independent_variables t
@variables x(t) y(t) z(t)

null_op = 0 * t
@test isequal(simplify(null_op), 0)

one_op = 1 * t
@test isequal(simplify(one_op), t)

identity_op = Num(Term(identity, [value(x)]))
@test isequal(simplify(identity_op), x)

minus_op = -x
@test isequal(simplify(minus_op), -1x)
simplify(minus_op)

@variables x

@test toexpr(expand_derivatives(Differential(x)((x - 2)^2))) == :($(*)(2, $(+)(-2, x)))
@test toexpr(expand_derivatives(Differential(x)((x - 2)^3))) ==
      :($(*)(3, $(^)($(+)(-2, x), 2)))
@test toexpr(simplify(x + 2 + 3)) == :($(+)(5, x))

d1 = Differential(x)((-2 + x)^2)
d2 = Differential(x)(d1)
d3 = Differential(x)(d2)

@test toexpr(expand_derivatives(d3)) == :(0)
@test toexpr(simplify(x^0)) == :(1)

@test ModelingToolkit.substitute(value(2x + y == 1), Dict(x => 0.0, y => 0.0)) === false
@test ModelingToolkit.substitute(value(2x + y == 1), Dict(x => 0.0, y => 1.0)) === true

# 699
using SymbolicUtils: substitute
@independent_variables t
@parameters a(t) b(t)

# back and forth substitution does not work for parameters with dependencies
term = value(a)
term2 = substitute(term, a => b)
@test ModelingToolkit.isparameter(term2)
@test isequal(term2, b)
term3 = substitute(term2, b => a)
@test ModelingToolkit.isparameter(term3)
@test isequal(term3, a)
