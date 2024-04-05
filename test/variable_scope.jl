using ModelingToolkit
using ModelingToolkit: SymScope, t_nounits as t, D_nounits as D
using Symbolics: arguments, value
using Test

@variables a b(t) c d e(t)

b = ParentScope(b)
c = ParentScope(ParentScope(c))
d = GlobalScope(d)
@test all(x -> x isa Num, [b, c, d])

# ensure it works on Term too
LocalScope(e.val)
ParentScope(e.val)
GlobalScope(e.val)

ie = ParentScope(1 / e)
@test getmetadata(arguments(value(ie))[2], SymScope) === ParentScope(LocalScope())

eqs = [0 ~ a
       0 ~ b
       0 ~ c
       0 ~ d]
@named sub4 = NonlinearSystem(eqs, [a, b, c, d], [])
@named sub3 = NonlinearSystem(eqs, [a, b, c, d], [])
@named sub2 = NonlinearSystem([], [], [], systems = [sub3, sub4])
@named sub1 = NonlinearSystem([], [], [], systems = [sub2])
@named sys = NonlinearSystem([], [], [], systems = [sub1])

names = ModelingToolkit.getname.(unknowns(sys))
@test :d in names
@test Symbol("sub1₊c") in names
@test Symbol("sub1₊sub2₊b") in names
@test Symbol("sub1₊sub2₊sub3₊a") in names
@test Symbol("sub1₊sub2₊sub4₊a") in names

@named foo = NonlinearSystem(eqs, [a, b, c, d], [])
@named bar = NonlinearSystem(eqs, [a, b, c, d], [])
@test ModelingToolkit.getname(ModelingToolkit.namespace_expr(
    ModelingToolkit.namespace_expr(b,
        foo),
    bar)) == Symbol("bar₊b")

function renamed(nss, sym)
    ModelingToolkit.getname(foldr(ModelingToolkit.renamespace, nss, init = sym))
end

@test renamed([:foo :bar :baz], a) == Symbol("foo₊bar₊baz₊a")
@test renamed([:foo :bar :baz], b) == Symbol("foo₊bar₊b")
@test renamed([:foo :bar :baz], c) == Symbol("foo₊c")
@test renamed([:foo :bar :baz], d) == :d

@parameters a b c d e f
p = [a
     ParentScope(b)
     ParentScope(ParentScope(c))
     DelayParentScope(d)
     DelayParentScope(e, 2)
     GlobalScope(f)]

level0 = ODESystem(Equation[], t, [], p; name = :level0)
level1 = ODESystem(Equation[], t, [], []; name = :level1) ∘ level0
level2 = ODESystem(Equation[], t, [], []; name = :level2) ∘ level1
level3 = ODESystem(Equation[], t, [], []; name = :level3) ∘ level2

ps = ModelingToolkit.getname.(parameters(level3))

@test isequal(ps[1], :level2₊level1₊level0₊a)
@test isequal(ps[2], :level2₊level1₊b)
@test isequal(ps[3], :level2₊c)
@test isequal(ps[4], :level2₊level0₊d)
@test isequal(ps[5], :level1₊level0₊e)
@test isequal(ps[6], :f)

@variables x(t) y(t)[1:2]
@parameters p q[1:2]

@test_throws ["Symbol", "x(t)", "does not occur"] ODESystem(
    [D(x) ~ p], t, [], [p]; name = :foo)
@test_nowarn ODESystem([D(x) ~ p], t, [x], [p]; name = :foo)
@test_throws ["Symbol", "y(t)", "does not occur"] ODESystem(
    D(y) ~ q, t, [], [q]; name = :foo)
@test_nowarn ODESystem(D(y) ~ q, t, [y], [q]; name = :foo)
@test_throws ["Symbol", "y(t)", "[1]", "does not occur"] ODESystem(
    D(y[1]) ~ x, t, [x], []; name = :foo)
@test_nowarn ODESystem(D(y[1]) ~ x, t, [x, y], []; name = :foo)
@test_throws ["Symbol", "p", "does not occur"] ODESystem(D(x) ~ p, t, [x], []; name = :foo)
@test_nowarn ODESystem(D(x) ~ p, t, [x], [p]; name = :foo)
@test_throws ["Symbol", "q", "does not occur"] ODESystem(D(y) ~ q, t, [y], []; name = :foo)
@test_nowarn ODESystem(D(y) ~ q, t, [y], [q]; name = :foo)
@test_throws ["Symbol", "q", "[1]", "does not occur"] ODESystem(
    D(y[1]) ~ q[1], t, [y], []; name = :foo)
@test_nowarn ODESystem(D(y[1]) ~ q[1], t, [y], [q]; name = :foo)
@test_throws ["Symbol", "x(t)", "does not occur"] ODESystem(
    Equation[], t, [], [p]; name = :foo, continuous_events = [[x ~ 0.0] => [p ~ 1.0]])
@test_nowarn ODESystem(
    Equation[], t, [x], [p]; name = :foo, continuous_events = [[x ~ 0.0] => [p ~ 1.0]])

@named sys1 = ODESystem(Equation[], t, [x, y], [p, q])
@test_throws ["Unexpected", "sys1₊x(t)", "subsystem with name sys1"] ODESystem(
    [D(x) ~ sys1.x], t; name = :sys2)
@test_nowarn ODESystem([D(x) ~ sys1.x], t; name = :sys2, systems = [sys1])
@test_throws ["Unexpected", "sys1₊y(t)", "subsystem with name sys1"] ODESystem(
    [D(x) ~ sum(sys1.y)], t; name = :sys2)
@test_nowarn ODESystem([D(x) ~ sum(sys1.y)], t; name = :sys2, systems = [sys1])
@test_throws ["Unexpected", "sys1₊y(t)", "[1]", "subsystem with name sys1"] ODESystem(
    D(x) ~ sys1.y[1], t; name = :sys2)
@test_nowarn ODESystem(D(x) ~ sys1.y[1], t; name = :sys2, systems = [sys1])
@test_throws ["Unexpected", "sys1₊p", "subsystem with name sys1"] ODESystem(
    D(x) ~ sys1.p, t; name = :sys2)
@test_nowarn ODESystem(D(x) ~ sys1.p, t; name = :sys2, systems = [sys1])
@test_throws ["Unexpected", "sys1₊q", "subsystem with name sys1"] ODESystem(
    D(y) ~ sys1.q, t; name = :sys2)
@test_nowarn ODESystem(D(y) ~ sys1.q, t; name = :sys2, systems = [sys1])
@test_throws ["Unexpected", "sys1₊q", "[1]", "subsystem with name sys1"] ODESystem(
    D(x) ~ sys1.q[1], t; name = :sys2)
@test_nowarn ODESystem(D(x) ~ sys1.q[1], t; name = :sys2, systems = [sys1])
@test_throws ["Unexpected", "sys1₊x(t)", "subsystem with name sys1"] ODESystem(
    Equation[], t, [], [p]; name = :sys2, continuous_events = [[sys1.x ~ 0] => [p ~ 1.0]])
@test_nowarn ODESystem(Equation[], t, [], [p]; name = :sys2,
    continuous_events = [[sys1.x ~ 0] => [p ~ 1.0]], systems = [sys1])

# Ensure SDESystem checks noise eqs as well
@test_throws ["Symbol", "x(t)", "does not occur"] SDESystem(
    Equation[], [0.1x], t, [], []; name = :foo)
@test_nowarn SDESystem(Equation[], [0.1x], t, [x], []; name = :foo)
@named sys1 = SDESystem(Equation[], [], t, [x], [])
@test_throws ["Unexpected", "sys1₊x(t)", "subsystem with name sys1"] SDESystem(
    Equation[], [0.1sys1.x], t, [], []; name = :foo)
@test_nowarn SDESystem(Equation[], [0.1sys1.x], t, [], []; name = :foo, systems = [sys1])

# Ensure DiscreteSystem checks work
k = ShiftIndex(t)
@test_throws ["Symbol", "x(t)", "does not occur"] DiscreteSystem(
    [x ~ x(k - 1) + x(k - 2)], t, [], []; name = :foo)
@test_nowarn DiscreteSystem([x ~ x(k - 1) + x(k - 2)], t; name = :foo)
@named sys1 = DiscreteSystem(Equation[], t, [x], [])
@test_throws ["Unexpected", "sys1₊x(t)", "subsystem with name sys1"] DiscreteSystem(
    [x ~ x(k - 1) + sys1.x(k - 2)], t, [x], []; name = :sys2)
@test_nowarn DiscreteSystem(
    [x ~ x(k - 1) + sys1.x(k - 2)], t, [x], []; name = :sys2, systems = [sys1])

# Ensure NonlinearSystem checks work
@variables x
@test_throws ["Symbol", "x", "does not occur"] NonlinearSystem(
    [0 ~ 2x + 3], [], []; name = :foo)
@test_nowarn NonlinearSystem([0 ~ 2x + 3], [x], []; name = :foo)
@named sys1 = NonlinearSystem(Equation[], [x], [])
@test_throws ["Unexpected", "sys1₊x", "subsystem with name sys1"] NonlinearSystem(
    [0 ~ sys1.x + 3], [], []; name = :foo)
@test_nowarn NonlinearSystem([0 ~ sys1.x + 3], [], []; name = :foo, systems = [sys1])

# Ensure ConstraintsSystem checks work
@test_throws ["Symbol", "x", "does not occur"] ConstraintsSystem(
    [0 ~ x^2 - 3], [], []; name = :foo)
@test_nowarn ConstraintsSystem([0 ~ x^2 - 3], [x], []; name = :foo)
@test_throws ["Symbol", "x", "does not occur"] ConstraintsSystem(
    [Inequality(x^2, 3, <)], [], []; name = :foo)
@test_nowarn ConstraintsSystem([Inequality(x^2, 3, <)], [x], []; name = :foo)
@named sys1 = ConstraintsSystem(Equation[], [x], [])
@test_throws ["Unexpected", "sys1₊x", "subsystem with name sys1"] ConstraintsSystem(
    [0 ~ sys1.x^2 - 2], [], []; name = :sys2)
@test_nowarn ConstraintsSystem([0 ~ sys1.x^2 - 2], [], []; name = :sys2, systems = [sys1])

# Ensure OptimizationSystem checks work
@test_throws ["Symbol", "x", "does not occur"] OptimizationSystem(
    y[1], [y[1]], []; constraints = [x ~ 3], name = :foo)
@test_nowarn OptimizationSystem(y[1], [y[1], x], []; constraints = [x ~ 3], name = :foo)
@named sys1 = OptimizationSystem(x, [x], [])
@test_throws ["Unexpected", "sys1₊x", "subsystem with name sys1"] OptimizationSystem(
    sys1.x^2 - 2, [], []; name = :sys2)
@test_nowarn OptimizationSystem(sys1.x^2 - 2, [], []; name = :sys2, systems = [sys1])
