using ModelingToolkit
using Test

using ModelingToolkit: value
using SymbolicUtils: FnType

const PType = ModelingToolkit.Parameter{Real}

@parameters t
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

t1 = Num(Sym{PType}(:t))
s1 = Num(Sym{PType}(:s))
σ1 = Num(Sym{FnType{Tuple, PType}}(:σ))
@test isequal(t1, t)
@test isequal(s1, s)
@test isequal(σ1, σ)

@derivatives D'~t
D1 = Differential(t)
@test D1 == D

@test @macroexpand(@parameters x, y, z(t)) == @macroexpand(@parameters x y z(t))
@test @macroexpand(@variables x, y, z(t)) == @macroexpand(@variables x y z(t))

# Test array expressions
@parameters begin
    t[1:2]
    s[1:2:4,1:2]
end
@parameters σ[1:2](..)

fntype(n, T) = FnType{NTuple{n, Any}, T}
t1 = Num[Variable{PType}(:t, 1), Variable{PType}(:t, 2)]
s1 = Num[Variable{PType}(:s, 1, 1) Variable{PType}(:s, 1, 2);
        Variable{PType}(:s, 3, 1) Variable{PType}(:s, 3, 2)]
σ1 = [Num(Variable{fntype(1, PType)}(:σ, 1)), Num(Variable{fntype(1, PType)}(:σ, 2))]
@test isequal(t1, t)
@test isequal(s1, s)
@test isequal(σ1, σ)

@parameters t
@variables x[1:2](t)
x1 = Num[Variable{FnType{Tuple{Any}, Real}}(:x, 1)(t.val),
      Variable{FnType{Tuple{Any}, Real}}(:x, 2)(t.val)]

@test isequal(x1, x)

@variables a[1:11,1:2]
@variables a()

using Symbolics: value, VariableDefaultValue
using ModelingToolkit: VariableConnectType, VariableUnit, rename
using Unitful

vals = [1,2,3,4]
@variables x=1 xs[1:4]=vals ys[1:5]=1

@test getmetadata(x, VariableDefaultValue) === 1
@test getmetadata.(xs, (VariableDefaultValue,)) == vals
@test getmetadata.(ys, (VariableDefaultValue,)) == ones(Int, 5)

struct Flow end
u = u"m^3/s"
@variables begin
    x = [1, 2], [connect=Flow,unit=u]
    y = 2
end

@test getmetadata(x, VariableDefaultValue) == [1, 2]
@test getmetadata(x, VariableConnectType) == Flow
@test getmetadata(x, VariableUnit) == u
@test getmetadata(y, VariableDefaultValue) === 2

@variables begin
    x, [connect=Flow,unit=u]
    y = 2, [connect=Flow]
end

@test !hasmetadata(x, VariableDefaultValue)
@test getmetadata(x, VariableConnectType) == Flow
@test getmetadata(x, VariableUnit) == u
@test getmetadata(y, VariableDefaultValue) === 2
@test getmetadata(y, VariableConnectType) == Flow

a = rename(value(x), :a)
@test !hasmetadata(x, VariableDefaultValue)
@test getmetadata(x, VariableConnectType) == Flow
@test getmetadata(x, VariableUnit) == u

@variables t x(t)=1 [connect=Flow,unit=u]

@test getmetadata(x, VariableDefaultValue) == 1
@test getmetadata(x, VariableConnectType) == Flow
@test getmetadata(x, VariableUnit) == u

a = rename(value(x), :a)
@test getmetadata(a, VariableDefaultValue) == 1
@test getmetadata(a, VariableConnectType) == Flow
@test getmetadata(a, VariableUnit) == u

@parameters p=2 [unit=u"m",]
@test getmetadata(p, VariableDefaultValue) == 2
@test !hasmetadata(p, VariableConnectType)
@test getmetadata(p, VariableUnit) == u"m"

@test_throws Any (@macroexpand @parameters p=2 [unit=u"m",abc=2])
