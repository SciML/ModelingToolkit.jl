using ModelingToolkit
using Test

@variables a b c

b = ParentScope(ParentScope(b))
c = GlobalScope(c)

renamed(nss, sym) = nameof(foldr(ModelingToolkit.renamespace, nss, init=sym))

@test renamed([:foo :bar :baz], a) == :foo₊bar₊baz₊a
@test renamed([:foo :bar :baz], b) == :foo₊b
@test renamed([:foo :bar :baz], c) == :c