using ModelingToolkit
using Test

using ModelingToolkit: value, Flow
using SymbolicUtils: FnType

@independent_variables t
@variables x(t) y(t) # test multi-arg
@variables z(t) # test single-arg

x1 = Num(Sym{FnType{Tuple{Any}, Real}}(:x)(value(t)))
y1 = Num(Sym{FnType{Tuple{Any}, Real}}(:y)(value(t)))
z1 = Num(Sym{FnType{Tuple{Any}, Real}}(:z)(value(t)))

@test isequal(x1, x)
@test isequal(y1, y)
@test isequal(z1, z)

@parameters begin
    t
    s
end
@parameters σ(..)

t1 = Num(Sym{Real}(:t))
s1 = Num(Sym{Real}(:s))
σ1 = Num(Sym{FnType{Tuple, Real}}(:σ))
@test isequal(t1, t)
@test isequal(s1, s)
@test isequal(σ1(t), σ(t))

@test ModelingToolkit.isparameter(t)
@test ModelingToolkit.isparameter(s)
@test ModelingToolkit.isparameter(σ)

@test @macroexpand(@parameters x, y, z(t)) == @macroexpand(@parameters x y z(t))
@test @macroexpand(@variables x, y, z(t)) == @macroexpand(@variables x y z(t))

# Test array expressions
@parameters begin
    t[1:2]
    s[1:4, 1:2]
end
@parameters σ(..)[1:2]

@test all(ModelingToolkit.isparameter, collect(t))
@test all(ModelingToolkit.isparameter, collect(s))
@test all(ModelingToolkit.isparameter, Any[σ(t)[1], σ(t)[2]])

# fntype(n, T) = FnType{NTuple{n, Any}, T}
# t1 = Num[Variable{Real}(:t, 1), Variable{Real}(:t, 2)]
# s1 = Num[Variable{Real}(:s, 1, 1) Variable{Real}(:s, 1, 2);
#         Variable{Real}(:s, 3, 1) Variable{Real}(:s, 3, 2)]
# σ1 = [Num(Variable{fntype(1, Real)}(:σ, 1)), Num(Variable{fntype(1, Real)}(:σ, 2))]
# @test isequal(t1, collect(t))
# @test isequal(s1, collect(s))
# @test isequal(σ1, σ)

#@independent_variables t
#@variables x[1:2](t)
#x1 = Num[Variable{FnType{Tuple{Any}, Real}}(:x, 1)(t.val),
#      Variable{FnType{Tuple{Any}, Real}}(:x, 2)(t.val)]
#
#@test isequal(x1, x)

@variables a[1:11, 1:2]
@variables a()

using Symbolics: value, VariableDefaultValue
using ModelingToolkit: VariableConnectType, VariableUnit, rename
using Unitful

vals = [1, 2, 3, 4]
@variables x=1 xs[1:4]=vals ys[1:5]=1

@test getmetadata(x, VariableDefaultValue) === 1
@test getmetadata.(collect(xs), (VariableDefaultValue,)) == vals
@test getmetadata.(collect(ys), (VariableDefaultValue,)) == ones(Int, 5)

u = u"m^3/s"
@variables begin
    x = [1, 2], [connect = Flow, unit = u]
    y = 2
end

@test getmetadata(x, VariableDefaultValue) == [1, 2]
@test getmetadata(x, VariableConnectType) == Flow
@test getmetadata(x, VariableUnit) == u
@test getmetadata(y, VariableDefaultValue) === 2

@variables x=[1, 2] [connect = Flow, unit = u] y=2

@test getmetadata(x, VariableDefaultValue) == [1, 2]
@test getmetadata(x, VariableConnectType) == Flow
@test getmetadata(x, VariableUnit) == u
@test getmetadata(y, VariableDefaultValue) === 2

@variables begin
    x, [connect = Flow, unit = u]
    y = 2, [connect = Flow]
end

@test_throws ErrorException ModelingToolkit.getdefault(x)
@test !hasmetadata(x, VariableDefaultValue)
@test getmetadata(x, VariableConnectType) == Flow
@test getmetadata(x, VariableUnit) == u
@test getmetadata(y, VariableDefaultValue) === 2
@test getmetadata(y, VariableConnectType) == Flow

a = rename(value(x), :a)
@test_throws ErrorException ModelingToolkit.getdefault(a)
@test !hasmetadata(a, VariableDefaultValue)
@test getmetadata(a, VariableConnectType) == Flow
@test getmetadata(a, VariableUnit) == u

@independent_variables t
@variables x(t)=1 [connect = Flow, unit = u]

@test getmetadata(x, VariableDefaultValue) == 1
@test getmetadata(x, VariableConnectType) == Flow
@test getmetadata(x, VariableUnit) == u

a = rename(value(x), :a)
@test getmetadata(a, VariableDefaultValue) == 1
@test getmetadata(a, VariableConnectType) == Flow
@test getmetadata(a, VariableUnit) == u

@parameters p=2 [unit = u"m"]
@test getmetadata(p, VariableDefaultValue) == 2
@test !hasmetadata(p, VariableConnectType)
@test getmetadata(p, VariableUnit) == u"m"
@test ModelingToolkit.isparameter(p)

@test_throws Any (@macroexpand @parameters p=2 [unit = u"m", abc = 2])
