using ModelingToolkitBase
using Test

using ModelingToolkitBase: value, Flow
using Symbolics: SSym
using SymbolicUtils: FnType, ShapeVecT

@independent_variables t
@variables x(t) y(t) # test multi-arg
@variables z(t) # test single-arg

x1 = Num(SSym(:x; type = FnType{Tuple, Real, Nothing}, shape = ShapeVecT())(value(t)))
y1 = Num(SSym(:y; type = FnType{Tuple, Real, Nothing}, shape = ShapeVecT())(value(t)))
z1 = Num(SSym(:z; type = FnType{Tuple, Real, Nothing}, shape = ShapeVecT())(value(t)))

@test isequal(x1, x)
@test isequal(y1, y)
@test isequal(z1, z)

@parameters begin
    t
    s
end
@parameters σ(..)

t1 = Num(SSym(:t; type = Real, shape = ShapeVecT()))
s1 = Num(SSym(:s; type = Real, shape = ShapeVecT()))
σ1 = SSym(:σ; type = FnType{Tuple, Real, Nothing}, shape = ShapeVecT())
@test isequal(t1, t)
@test isequal(s1, s)
@test isequal(σ1(t), σ(t))

@test ModelingToolkitBase.isparameter(t)
@test ModelingToolkitBase.isparameter(s)
@test ModelingToolkitBase.isparameter(σ)

@test @macroexpand(@parameters x, y, z(t)) == @macroexpand(@parameters x y z(t))
@test @macroexpand(@variables x, y, z(t)) == @macroexpand(@variables x y z(t))

# Test array expressions
@parameters begin
    t[1:2]
    s[1:4, 1:2]
end
@parameters σ(..)[1:2]

@test all(ModelingToolkitBase.isparameter, collect(t))
@test all(ModelingToolkitBase.isparameter, collect(s))
@test all(ModelingToolkitBase.isparameter, Any[σ(t)[1], σ(t)[2]])

@variables a[1:11, 1:2]
@variables a()

using Symbolics: value, VariableDefaultValue, getdefaultval
using ModelingToolkitBase: VariableConnectType, VariableUnit, rename
using DynamicQuantities

vals = [1, 2, 3, 4]
@variables x = 1 xs[1:4] = vals ys[1:5] = ones(5)

@test getmetadata(x, VariableDefaultValue) === 1
@test getdefaultval(xs) == vals
@test getdefaultval(ys) == ones(Int, 5)

u = u"m^3/s"
@variables begin
    x = [1, 2], [connect = Flow, unit = u]
    y = 2
end

@test getmetadata(x, VariableDefaultValue) == [1, 2]
@test getmetadata(x, VariableConnectType) == Flow
@test getmetadata(x, VariableUnit) == u
@test getmetadata(y, VariableDefaultValue) === 2

@variables x = [1, 2] [connect = Flow, unit = u] y = 2

@test getmetadata(x, VariableDefaultValue) == [1, 2]
@test getmetadata(x, VariableConnectType) == Flow
@test getmetadata(x, VariableUnit) == u
@test getmetadata(y, VariableDefaultValue) === 2

@variables begin
    x, [connect = Flow, unit = u]
    y = 2, [connect = Flow]
end

@test_throws ErrorException ModelingToolkitBase.getdefault(x)
@test !hasmetadata(x, VariableDefaultValue)
@test getmetadata(x, VariableConnectType) == Flow
@test getmetadata(x, VariableUnit) == u
@test getmetadata(y, VariableDefaultValue) === 2
@test getmetadata(y, VariableConnectType) == Flow

a = rename(value(x), :a)
@test_throws ErrorException ModelingToolkitBase.getdefault(a)
@test !hasmetadata(a, VariableDefaultValue)
@test getmetadata(a, VariableConnectType) == Flow
@test getmetadata(a, VariableUnit) == u

@independent_variables t
@variables x(t) = 1 [connect = Flow, unit = u]

@test getmetadata(x, VariableDefaultValue) == 1
@test getmetadata(x, VariableConnectType) == Flow
@test getmetadata(x, VariableUnit) == u

a = rename(value(x), :a)
@test getmetadata(a, VariableDefaultValue) == 1
@test getmetadata(a, VariableConnectType) == Flow
@test getmetadata(a, VariableUnit) == u

@parameters p = 2 [unit = u"m"]
@test getmetadata(p, VariableDefaultValue) == 2
@test !hasmetadata(p, VariableConnectType)
@test getmetadata(p, VariableUnit) == u"m"
@test ModelingToolkitBase.isparameter(p)

@test_throws Any (@macroexpand @parameters p = 2 [unit = u"m", abc = 2])

@testset "Parameters cannot be dependent" begin
    @test_throws ["cannot create time-dependent"] @parameters p(t)
    @test_throws ["cannot create time-independent"] @discretes p
end
