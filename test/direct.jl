using ModelingToolkit, StaticArrays, LinearAlgebra, SparseArrays
using DiffEqBase
using Test

using ModelingToolkit: getdefault, getmetadata, SymScope

canonequal(a, b) = isequal(simplify(a), simplify(b))

# Calculus
@parameters t σ ρ β
@variables x y z
@test isequal((Differential(z) * Differential(y) * Differential(x))(t),
    Differential(z)(Differential(y)(Differential(x)(t))))

@test canonequal(ModelingToolkit.derivative(sin(cos(x)), x),
    -sin(x) * cos(cos(x)))

@register_symbolic no_der(x)
@test canonequal(ModelingToolkit.derivative([sin(cos(x)), hypot(x, no_der(x))], x),
    [
        -sin(x) * cos(cos(x)),
        x / hypot(x, no_der(x)) +
        no_der(x) * Differential(x)(no_der(x)) / hypot(x, no_der(x)),
    ])

@register_symbolic intfun(x)::Int
@test ModelingToolkit.symtype(intfun(x)) === Int

eqs = [σ * (y - x),
    x * (ρ - z) - y,
    x * y - β * z]

simpexpr = [:($(*)(σ, $(+)(y, $(*)(-1, x))))
    :($(+)($(*)(x, $(+)(ρ, $(*)(-1, z))), $(*)(-1, y)))
    :($(+)($(*)(x, y), $(*)(-1, z, β)))]

σ, β, ρ = 2 // 3, 3 // 4, 4 // 5
x, y, z = 6 // 7, 7 // 8, 8 // 9
for i in 1:3
    @test eval(ModelingToolkit.toexpr.(eqs)[i]) == eval(simpexpr[i])
    @test eval(ModelingToolkit.toexpr.(eqs)[i]) == eval(simpexpr[i])
end

@parameters t σ ρ β
@variables x y z
∂ = ModelingToolkit.jacobian(eqs, [x, y, z])
for i in 1:3
    ∇ = ModelingToolkit.gradient(eqs[i], [x, y, z])
    @test canonequal(∂[i, :], ∇)
end

@test all(canonequal.(ModelingToolkit.gradient(eqs[1], [x, y, z]), [σ * -1, σ, 0]))
@test all(canonequal.(ModelingToolkit.hessian(eqs[1], [x, y, z]), 0))

du = [x^2, y^3, x^4, sin(y), x + y, x + z^2, z + x, x + y^2 + sin(z)]
reference_jac = sparse(ModelingToolkit.jacobian(du, [x, y, z]))

@test findnz(ModelingToolkit.jacobian_sparsity(du, [x, y, z]))[[1, 2]] ==
      findnz(reference_jac)[[1, 2]]

let
    @variables t x(t) y(t) z(t)
    @test ModelingToolkit.exprs_occur_in([x, y, z], x^2 * y) == [true, true, false]
end

@test isequal(ModelingToolkit.sparsejacobian(du, [x, y, z]), reference_jac)

using ModelingToolkit

rosenbrock(X) =
    sum(1:(length(X) - 1)) do i
        100 * (X[i + 1] - X[i]^2)^2 + (1 - X[i])^2
    end

@variables a, b
X = [a, b]

spoly(x) = simplify(x, expand = true)
rr = rosenbrock(X)

reference_hes = ModelingToolkit.hessian(rr, X)
@test findnz(sparse(reference_hes))[1:2] ==
      findnz(ModelingToolkit.hessian_sparsity(rr, X))[1:2]

sp_hess = ModelingToolkit.sparsehessian(rr, X)
@test findnz(sparse(reference_hes))[1:2] == findnz(sp_hess)[1:2]
@test isequal(map(spoly, findnz(sparse(reference_hes))[3]), map(spoly, findnz(sp_hess)[3]))

Joop, Jiip = eval.(ModelingToolkit.build_function(∂, [x, y, z], [σ, ρ, β], t))
J = Joop([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], 1.0)
@test J isa Matrix
J2 = copy(J)
Jiip(J2, [1.0, 2.0, 3.0], [1.0, 2.0, 3.0], 1.0)
@test J2 == J

Joop, Jiip = eval.(ModelingToolkit.build_function(vcat(∂, ∂), [x, y, z], [σ, ρ, β], t))
J = Joop([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], 1.0)
@test J isa Matrix
J2 = copy(J)
Jiip(J2, [1.0, 2.0, 3.0], [1.0, 2.0, 3.0], 1.0)
@test J2 == J

Joop, Jiip = eval.(ModelingToolkit.build_function(hcat(∂, ∂), [x, y, z], [σ, ρ, β], t))
J = Joop([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], 1.0)
@test J isa Matrix
J2 = copy(J)
Jiip(J2, [1.0, 2.0, 3.0], [1.0, 2.0, 3.0], 1.0)
@test J2 == J

∂3 = cat(∂, ∂, dims = 3)
Joop, Jiip = eval.(ModelingToolkit.build_function(∂3, [x, y, z], [σ, ρ, β], t))
J = Joop([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], 1.0)
@test size(J) == (3, 3, 2)
J2 = copy(J)
Jiip(J2, [1.0, 2.0, 3.0], [1.0, 2.0, 3.0], 1.0)
@test J2 == J

s∂ = sparse(∂)
@test nnz(s∂) == 8
Joop, Jiip = eval.(ModelingToolkit.build_function(s∂, [x, y, z], [σ, ρ, β], t,
    linenumbers = true))
J = Joop([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], 1.0)
@test length(nonzeros(s∂)) == 8
J2 = copy(J)
Jiip(J2, [1.0, 2.0, 3.0], [1.0, 2.0, 3.0], 1.0)
@test J2 == J

# Function building

@parameters σ ρ β
@variables x y z
eqs = [σ * (y - x),
    x * (ρ - z) - y,
    x * y - β * z]
f1, f2 = ModelingToolkit.build_function(eqs, [x, y, z], [σ, ρ, β])
f = eval(f1)
out = [1.0, 2, 3]
o1 = f([1.0, 2, 3], [1.0, 2, 3])
f = eval(f2)
f(out, [1.0, 2, 3], [1.0, 2, 3])
@test all(o1 .== out)

function test_worldage()
    @parameters σ ρ β
    @variables x y z
    eqs = [σ * (y - x),
        x * (ρ - z) - y,
        x * y - β * z]
    f, f_iip = ModelingToolkit.build_function(eqs, [x, y, z], [σ, ρ, β];
        expression = Val{false})
    out = [1.0, 2, 3]
    o1 = f([1.0, 2, 3], [1.0, 2, 3])
    f_iip(out, [1.0, 2, 3], [1.0, 2, 3])
end
test_worldage()

## No parameters
@variables x y z
eqs = [(y - x)^2,
    x * (x - z) - y,
    x * y - y * z]
f1, f2 = ModelingToolkit.build_function(eqs, [x, y, z])
f = eval(f1)
out = zeros(3)
o1 = f([1.0, 2, 3])
f = eval(f2)
f(out, [1.0, 2, 3])
@test all(out .== o1)

# y ^ -1 test
g = let
    f(x, y) = x / y
    @variables x y
    ex = expand_derivatives(Differential(x)(f(x, y)))
    func_ex = build_function(ex, x, y)
    eval(func_ex)
end

@test g(42, 4) == 1 / 4

function test_worldage()
    @variables x y z
    eqs = [(y - x)^2,
        x * (x - z) - y,
        x * y - y * z]
    f, f_iip = ModelingToolkit.build_function(eqs, [x, y, z]; expression = Val{false})
    out = zeros(3)
    o1 = f([1.0, 2, 3])
    f_iip(out, [1.0, 2, 3])
end
test_worldage()

@test_nowarn muladd(x, y, 0)
@test promote(x, 0) == (x, identity(0))
@test_nowarn [x, y, z]'

let
    @register_symbolic foo(x)
    @variables t
    D = Differential(t)

    @test isequal(expand_derivatives(D(foo(t))), D(foo(t)))
    @test isequal(expand_derivatives(D(sin(t) * foo(t))),
        cos(t) * foo(t) + sin(t) * D(foo(t)))
end

foo(; kw...) = kw
foo(args...; kw...) = args, kw
pp = :name => :cool_name

@named cool_name = foo()
@test collect(cool_name) == [pp]

@named cool_name = foo(42)
@test cool_name[1] == (42,)
@test collect(cool_name[2]) == [pp]

@named cool_name = foo(42; a = 2)
@test cool_name[1] == (42,)
@test collect(cool_name[2]) == [pp; :a => 2]

@named cool_name = foo(a = 2)
@test collect(cool_name) == [pp; :a => 2]

@named cool_name = foo(; a = 2)
@test collect(cool_name) == [pp; :a => 2]

@named cool_name = foo(name = 2)
@test collect(cool_name) == [:name => 2]

@named cool_name = foo(42; name = 3)
@test cool_name[1] == (42,)
@test collect(cool_name[2]) == [:name => 3]

kwargs = (; name = 3)
@named cool_name = foo(42; kwargs...)
@test cool_name[1] == (42,)
@test collect(cool_name[2]) == [:name => 3]

if VERSION >= v"1.5"
    name = 3
    @named cool_name = foo(42; name)
    @test cool_name[1] == (42,)
    @test collect(cool_name[2]) == [:name => name]
    @named cool_name = foo(; name)
    @test collect(cool_name) == [:name => name]

    ff = 3
    @named cool_name = foo(42; ff)
    @test cool_name[1] == (42,)
    @test collect(cool_name[2]) == [pp; :ff => ff]

    @named cool_name = foo(; ff)
    @test collect(cool_name) == [pp; :ff => ff]
end

foo(i; name) = (; i, name)
@named goo[1:3] = foo(10)
@test isequal(goo, [(i = 10, name = Symbol(:goo_, i)) for i in 1:3])
@named koo 1:3 i->foo(10i)
@test isequal(koo, [(i = 10i, name = Symbol(:koo_, i)) for i in 1:3])
xys = @named begin
    x = foo(12)
    y[1:3] = foo(13)
end
@test isequal(x, (i = 12, name = :x))
@test isequal(y, [(i = 13, name = Symbol(:y_, i)) for i in 1:3])
@test isequal(xys, [x; y])

@variables x [misc = "wow"]
@test SymbolicUtils.getmetadata(Symbolics.unwrap(x), ModelingToolkit.VariableMisc,
    nothing) == "wow"
@parameters x [misc = "wow"]
@test SymbolicUtils.getmetadata(Symbolics.unwrap(x), ModelingToolkit.VariableMisc,
    nothing) == "wow"

# Scope of defaults in the systems generated by @named
@mtkmodel MoreThanOneArg begin
    @variables begin
        x(t)
        y(t)
        z(t)
    end
end

@parameters begin
    l
    m
    n
end

@named model = MoreThanOneArg(x = l, y = m, z = n)

@test getmetadata(getdefault(model.x), SymScope) == ParentScope(LocalScope())
@test getmetadata(getdefault(model.y), SymScope) == ParentScope(LocalScope())
@test getmetadata(getdefault(model.z), SymScope) == ParentScope(LocalScope())
