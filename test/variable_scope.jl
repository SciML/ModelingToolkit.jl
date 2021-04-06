using ModelingToolkit
using Test

@parameters t
@variables a b(t) c d

b = ParentScope(b)
c = ParentScope(ParentScope(c))
d = GlobalScope(d)

renamed(nss, sym) = ModelingToolkit.getname(foldr(ModelingToolkit.renamespace, nss, init=sym))

@test renamed([:foo :bar :baz], a) == :foo₊bar₊baz₊a
@test renamed([:foo :bar :baz], b) == :foo₊bar₊b
@test renamed([:foo :bar :baz], c) == :foo₊c
@test renamed([:foo :bar :baz], d) == :d

eqs = [
    0 ~ a
    0 ~ b
    0 ~ c
    0 ~ d
]
@named sub4 = NonlinearSystem(eqs, [a, b, c, d], [])
@named sub3 = NonlinearSystem(eqs, [a, b, c, d], [])
@named sub2 = NonlinearSystem([], [], [], systems=[sub3, sub4])
@named sub1 = NonlinearSystem([], [], [], systems=[sub2])
@named sys = NonlinearSystem([], [], [], systems=[sub1])

names = ModelingToolkit.getname.(states(sys))
@test :d in names
@test :sub1₊c in names
@test :sub1₊sub2₊b in names
@test :sub1₊sub2₊sub3₊a in names
@test :sub1₊sub2₊sub4₊a in names