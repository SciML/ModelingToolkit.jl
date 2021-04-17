using Test, ModelingToolkit

@parameters t

@connector function Foo(;name)
    @variables x(t)
    ODESystem(Equation[], t, [x], [], defaults=Dict(x=>1.0))
end

@connector function Goo(;name)
    @variables x(t)
    @parameters p
    ODESystem(Equation[], t, [x], [p], defaults=Dict(x=>1.0, p=>1.0))
end

ModelingToolkit.connect(::Type{<:Foo}, sys1, sys2) = [sys1.x ~ sys2.x]
@named f1 = Foo()
@named f2 = Foo()
@named g = Goo()

@test isequal(connect(f1, f2), [f1.x ~ f2.x])
@test_throws ArgumentError connect(f1, g)

# Note that since there're overloadings, these tests are not re-runable.
ModelingToolkit.promote_connect_rule(::Type{<:Foo}, ::Type{<:Goo}) = Foo
@test isequal(connect(f1, g), [f1.x ~ g.x])
ModelingToolkit.promote_connect_rule(::Type{<:Goo}, ::Type{<:Foo}) = Foo
@test isequal(connect(f1, g), [f1.x ~ g.x])
# test conflict
ModelingToolkit.promote_connect_rule(::Type{<:Goo}, ::Type{<:Foo}) = Goo
@test_throws ArgumentError connect(f1, g)
