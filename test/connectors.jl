using Test, ModelingToolkit

@parameters t

@connector function Foo(;name)
    @variables x(t)
    ODESystem(Equation[], t, [x], [], defaults=Dict(x=>1.0), name=name)
end

@connector function Goo(;name)
    @variables x(t)
    @parameters p
    ODESystem(Equation[], t, [x], [p], defaults=Dict(x=>1.0, p=>1.0), name=name)
end

function ModelingToolkit.connect(::Type{<:Foo}, ss...)
    n = length(ss)-1
    eqs = Vector{Equation}(undef, n)
    for i in 1:n
        eqs[i] = ss[i].x ~ ss[i+1].x
    end
    eqs
end

@named f1 = Foo()
@named f2 = Foo()
@named f3 = Foo()
@named f4 = Foo()
@named g = Goo()

@test isequal(connect(f1, f2), [f1.x ~ f2.x])
@test_throws ArgumentError connect(f1, g)

# Note that since there're overloadings, these tests are not re-runable.
ModelingToolkit.promote_connect_rule(::Type{<:Foo}, ::Type{<:Goo}) = Foo
@test isequal(connect(f1, g), [f1.x ~ g.x])
@test isequal(connect(f1, f2, g), [f1.x ~ f2.x; f2.x ~ g.x])
@test isequal(connect(f1, f2, g, f3), [f1.x ~ f2.x; f2.x ~ g.x; g.x ~ f3.x])
@test isequal(connect(f1, f2, g, f3, f4), [f1.x ~ f2.x; f2.x ~ g.x; g.x ~ f3.x; f3.x ~ f4.x])
ModelingToolkit.promote_connect_rule(::Type{<:Goo}, ::Type{<:Foo}) = Foo
@test isequal(connect(f1, g), [f1.x ~ g.x])
# test conflict
ModelingToolkit.promote_connect_rule(::Type{<:Goo}, ::Type{<:Foo}) = Goo
@test_throws ArgumentError connect(f1, g)
