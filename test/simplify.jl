using ModelingToolkit
using ModelingToolkit: value, get_defaults, get_default_u0, get_default_p
using Test
using DifferentialEquations

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

@test toexpr(expand_derivatives(Differential(x)((x-2)^2))) == :($(*)(2, $(+)(-2, x)))
@test toexpr(expand_derivatives(Differential(x)((x-2)^3))) == :($(*)(3, $(^)($(+)(-2, x), 2)))
@test toexpr(simplify(x+2+3)) == :($(+)(5, x))

d1 = Differential(x)((-2 + x)^2)
d2 = Differential(x)(d1)
d3 = Differential(x)(d2)

@test toexpr(expand_derivatives(d3)) == :(0)
@test toexpr(simplify(x^0)) == :(1)

@test ModelingToolkit.substitute(value(2x + y == 1), Dict(x => 0.0, y => 0.0)) === false
@test ModelingToolkit.substitute(value(2x + y == 1), Dict(x => 0.0, y => 1.0)) === true

# 699
using SymbolicUtils: substitute
@parameters t a(t) b(t)

# back and forth substitution does not work for parameters with dependencies
term = value(a)
term2 = substitute(term, a=>b)
@test term2 isa Term{ModelingToolkit.Parameter{Real}}
@test isequal(term2, b)
term3 = substitute(term2, b=>a)
@test term3 isa Term{ModelingToolkit.Parameter{Real}}
@test isequal(term3, a)




@parameters t
D = Differential(t)
@variables x(t),y(t)
y0 = x/2
@named sub = ODESystem([D(y) ~ 2x],t,[x,y],[],defaults=Dict(y=>y0))

@variables a(t)
@named model = ODESystem([D(a) ~ 1, sub.x ~ a],t,[a],[],defaults=Dict(a=>0.1),systems=[sub])

sys = structural_simplify(model)
defs = get_defaults(sys)

_u0 = filter(e->!ModelingToolkit.isparameter(e.first),defs)
_p = filter(e->ModelingToolkit.isparameter(e.first),defs)
@test isequal(_u0, Dict([a=>0.1, sub.y=>0.5sub.x]))
prob = ODEProblem(sys,collect(_u0),(0.0,160.0),collect(_p),jac=true)

# SHOULD NOT HAPPEN; TODO: remove this test
@test_throws StackOverflowError solve(prob)

u0 = get_default_u0(sys)
p = get_default_p(sys)
@test isequal(u0, Dict([a=>0.1, sub.y=>0.05]))
@test isequal(p, Dict())
prob = ODEProblem(sys,collect(u0),(0.0,160.0),collect(p),jac=true)
sol = solve(prob)
@test isequal(sol[a], sol[sub.x])

