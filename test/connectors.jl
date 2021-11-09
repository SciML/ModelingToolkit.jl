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

function ModelingToolkit.connect(::Type{<:Foo}, c::Connection)
    @show c.inners
    ss = c.inners
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

function connection_eqs(eqs, subsys)
    @named sys = ODESystem(eqs, t)
    @named newsys = compose(sys, subsys)
    equations(expand_connections(newsys))
end

connection_eqs(subsys) = Base.Fix2(connection_eqs, subsys)
ceqs = connection_eqs([f1, f2, f3, f4, g])

@test isequal(ceqs(connect(f1, f2)), [f1.x ~ f2.x])
@test_throws ArgumentError ceqs(connect(f1, g))

# Note that since there're overloadings, these tests are not re-runable.
ModelingToolkit.promote_connect_rule(::Type{<:Foo}, ::Type{<:Goo}) = Foo
@test isequal(ceqs(connect(f1, g)), [f1.x ~ g.x])
@test isequal(ceqs(connect(f1, f2, g)), [f1.x ~ f2.x; f2.x ~ g.x])
@test isequal(ceqs(connect(f1, f2, g, f3)), [f1.x ~ f2.x; f2.x ~ g.x; g.x ~ f3.x])
@test isequal(ceqs(connect(f1, f2, g, f3, f4)), [f1.x ~ f2.x; f2.x ~ g.x; g.x ~ f3.x; f3.x ~ f4.x])
ModelingToolkit.promote_connect_rule(::Type{<:Goo}, ::Type{<:Foo}) = Foo
@test isequal(ceqs(connect(f1, g)), [f1.x ~ g.x])
# test conflict
ModelingToolkit.promote_connect_rule(::Type{<:Goo}, ::Type{<:Foo}) = Goo
@test_throws ArgumentError ceqs(connect(f1, g))

@connector Hoo(;name) = ODESystem(Equation[], t, [], [], name=name)
function ModelingToolkit.connect(::Type{<:Hoo}, c::Connection)
    ss = c.inners
    nameof.(ss) ~ 0
end
@named hs[1:8] = Hoo()
ceqs = connection_eqs(hs)

@test ceqs([connect(hs[1], hs[2]),
            connect(hs[1], hs[3])]) == [[:hs_1, :hs_2, :hs_3] ~ 0]

@test ceqs([connect(hs[1], hs[2]),
            connect(hs[2], hs[3])]) == [[:hs_1, :hs_2, :hs_3] ~ 0]

@test ceqs([connect(hs[1], hs[2]),
            connect(hs[4], hs[3])]) == [[:hs_1, :hs_2] ~ 0, [:hs_4, :hs_3] ~ 0]
@test_throws Any ceqs([connect(hs[1], hs[2]),
                       connect(hs[1], hs[2])])
@test ceqs([connect(hs[1], hs[2]),
            connect(hs[3], hs[2]),
            connect(hs[1], hs[4]),
            connect(hs[8], hs[4]),
            connect(hs[7], hs[5]),]) == [[:hs_1, :hs_2, :hs_3, :hs_4, :hs_8] ~ 0, [:hs_7, :hs_5] ~ 0]
@test_throws Any ceqs([connect(hs[1], hs[2]),
                       connect(hs[3], hs[2]),
                       connect(hs[1], hs[4]),
                       connect(hs[8], hs[4]),
                       connect(hs[2], hs[8])])

# Outer/inner connectors

