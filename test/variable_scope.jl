using ModelingToolkit
using ModelingToolkit: SymScope
using Symbolics: arguments, value
using Test

@independent_variables t
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

@independent_variables t
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

# Issue@2252
# Tests from PR#2354
@parameters xx[1:2]
arr_p = [ParentScope(xx[1]), xx[2]]
arr0 = ODESystem(Equation[], t, [], arr_p; name = :arr0)
arr1 = ODESystem(Equation[], t, [], []; name = :arr1) ∘ arr0
arr_ps = ModelingToolkit.getname.(parameters(arr1))
@test isequal(arr_ps[1], Symbol("xx"))
@test isequal(arr_ps[2], Symbol("arr0₊xx"))
