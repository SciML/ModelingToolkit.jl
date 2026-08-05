using ModelingToolkitBase, Test
using ModelingToolkitBase: value, parse_variable
using OrderedCollections: OrderedSet
using SymbolicUtils: <ₑ
import SymbolicUtils as SU
import Symbolics

@variables x
@test ModelingToolkitBase.float_type_from_varmap(
    [x => Union{Nothing, BigFloat}[nothing]]
) == BigFloat
@test ModelingToolkitBase.float_type_from_varmap(
    [x => Union{Nothing, Float32}[nothing]]
) == Float32

@parameters α β δ
expr = (((1 / β - 1) + δ) / α)^(1 / (α - 1))
ref = sort([β, δ, α], lt = <ₑ)
sol = sort(Num.(ModelingToolkitBase.get_variables(expr)), lt = <ₑ)
@test all(x -> x isa Num, sol[i] == ref[i] for i in 1:3)
@test all(isequal(sol[i], ref[i]) for i in 1:3)

@testset "collect_expression_vars!" begin
    @independent_variables τ
    @variables X(τ) Y(τ) Z(τ) = α
    @parameters k_obs k_ratio eps_val

    us = OrderedSet{Symbolics.SymbolicT}()
    ps = OrderedSet{Symbolics.SymbolicT}()
    ModelingToolkitBase.collect_expression_vars!(
        us, ps, k_obs + k_ratio * (X + eps_val) + τ, Symbolics.unwrap(τ)
    )
    @test issetequal(us, Symbolics.unwrap.([X]))
    @test issetequal(ps, Symbolics.unwrap.([k_obs, k_ratio, eps_val]))

    empty!(us)
    empty!(ps)
    ModelingToolkitBase.collect_expression_vars!(us, ps, Y + Z, Symbolics.unwrap(τ))
    @test issetequal(us, Symbolics.unwrap.([Y, Z]))
    @test Symbolics.unwrap(α) in ps
end

@parameters γ
s = α => γ
expr = (((1 / β - 1) + δ) / α)^(1 / (α - 1))
sol = ModelingToolkitBase.substitute(expr, s)
new = (((1 / β - 1) + δ) / γ)^(1 / (γ - 1))
@test iszero(sol - new)

# Continuous
using ModelingToolkitBase: isdifferential, collect_differential_variables,
    collect_ivs
@independent_variables t
@variables u(t) y(t)
D = Differential(t)
eq = D(y) ~ u
v = SU.search_variables(eq)
@test v == Set([D(y), u])

ov = collect_differential_variables(eq)
@test ov == Set(Any[y])

aov = ModelingToolkitBase.collect_applied_operators(eq, Differential)
@test aov == Set(Any[D(y)])

ts = collect_ivs([eq])
@test ts == Set([t])

@testset "parse_variable with scalarized arrays" begin
    @variables scalarized_x(t)[1:2]
    @parameters scalarized_p[1:2]
    scalarized_sys = System(
        Equation[], t, collect(scalarized_x), collect(scalarized_p); name = :scalarized
    )
    @test isequal(parse_variable(scalarized_sys, "scalarized_x"), scalarized_x)
    @test isequal(parse_variable(scalarized_sys, "scalarized_x[2]"), scalarized_x[2])
    @test isequal(parse_variable(scalarized_sys, "scalarized_p"), scalarized_p)
    @test isequal(parse_variable(scalarized_sys, "scalarized_p[2]"), scalarized_p[2])
end

@testset "parse_variable with iv: $iv" for iv in [t, only(@independent_variables tt)]
    D = Differential(iv)
    function Lorenz(; name)
        @variables begin
            x(iv)
            y(iv)
            z(iv)
        end
        @parameters begin
            σ
            ρ
            β
        end
        sys = System(
            [
                D(D(x)) ~ σ * (y - x)
                D(y) ~ x * (ρ - z) - y
                D(z) ~ x * y - β * z
            ], iv; name
        )
    end
    function ArrSys(; name)
        @variables begin
            x(iv)[1:2]
        end
        @parameters begin
            p[1:2, 1:2]
        end
        sys = System([D(D(x)) ~ p * x], iv; name)
    end
    function Outer(; name)
        @named 😄 = Lorenz()
        @named arr = ArrSys()
        sys = System(Equation[], iv; name, systems = [😄, arr])
    end

    @mtkcompile sys = Outer()
    for (str, var) in [
            # unicode system, scalar variable
            ("😄.x", sys.😄.x),
            ("😄.x($iv)", sys.😄.x),
            ("😄₊x", sys.😄.x),
            ("😄₊x($iv)", sys.😄.x),
            # derivative
            ("D(😄.x)", D(sys.😄.x)),
            ("D(😄.x($iv))", D(sys.😄.x)),
            ("D(😄₊x)", D(sys.😄.x)),
            ("D(😄₊x($iv))", D(sys.😄.x)),
            ("Differential($iv)(😄.x)", D(sys.😄.x)),
            ("Differential($iv)(😄.x($iv))", D(sys.😄.x)),
            ("Differential($iv)(😄₊x)", D(sys.😄.x)),
            ("Differential($iv)(😄₊x($iv))", D(sys.😄.x)),
            # other derivative
            ("😄.xˍ$iv", D(sys.😄.x)),
            ("😄.x($iv)ˍ$iv", D(sys.😄.x)),
            ("😄₊xˍ$iv", D(sys.😄.x)),
            ("😄₊x($iv)ˍ$iv", D(sys.😄.x)),
            # scalar parameter
            ("😄.σ", sys.😄.σ),
            ("😄₊σ", sys.😄.σ),
            # array variable
            ("arr.x", sys.arr.x),
            ("arr₊x", sys.arr.x),
            ("arr.x($iv)", sys.arr.x),
            ("arr₊x($iv)", sys.arr.x),
            # getindex
            ("arr.x[1]", sys.arr.x[1]),
            ("arr₊x[1]", sys.arr.x[1]),
            ("arr.x($iv)[1]", sys.arr.x[1]),
            ("arr₊x($iv)[1]", sys.arr.x[1]),
            # derivative
            ("D(arr.x($iv))", D(sys.arr.x)),
            ("D(arr₊x($iv))", D(sys.arr.x)),
            ("D(arr.x[1])", D(sys.arr.x[1])),
            ("D(arr₊x[1])", D(sys.arr.x[1])),
            ("D(arr.x($iv)[1])", D(sys.arr.x[1])),
            ("D(arr₊x($iv)[1])", D(sys.arr.x[1])),
            ("Differential($iv)(arr.x($iv))", D(sys.arr.x)),
            ("Differential($iv)(arr₊x($iv))", D(sys.arr.x)),
            ("Differential($iv)(arr.x[1])", D(sys.arr.x[1])),
            ("Differential($iv)(arr₊x[1])", D(sys.arr.x[1])),
            ("Differential($iv)(arr.x($iv)[1])", D(sys.arr.x[1])),
            ("Differential($iv)(arr₊x($iv)[1])", D(sys.arr.x[1])),
            # other derivative
            ("arr.xˍ$iv", D(sys.arr.x)),
            ("arr₊xˍ$iv", D(sys.arr.x)),
            ("arr.xˍ$iv($iv)", D(sys.arr.x)),
            ("arr₊xˍ$iv($iv)", D(sys.arr.x)),
            ("arr.xˍ$iv[1]", D(sys.arr.x[1])),
            ("arr₊xˍ$iv[1]", D(sys.arr.x[1])),
            ("arr.xˍ$iv($iv)[1]", D(sys.arr.x[1])),
            ("arr₊xˍ$iv($iv)[1]", D(sys.arr.x[1])),
            ("arr.x($iv)ˍ$iv", D(sys.arr.x)),
            ("arr₊x($iv)ˍ$iv", D(sys.arr.x)),
            ("arr.x($iv)ˍ$iv[1]", D(sys.arr.x[1])),
            ("arr₊x($iv)ˍ$iv[1]", D(sys.arr.x[1])),
            # array parameter
            ("arr.p", sys.arr.p),
            ("arr₊p", sys.arr.p),
            ("arr.p[1, 2]", sys.arr.p[1, 2]),
            ("arr₊p[1, 2]", sys.arr.p[1, 2]),
        ]
        @test isequal(parse_variable(sys, str), var)
    end
end

@testset "isinitial" begin
    t = ModelingToolkitBase.t_nounits
    @variables x(t) z(t)[1:5]
    @parameters a b c[1:4]
    @test isinitial(Initial(z))
    @test isinitial(Initial(x))
    @test isinitial(Initial(a))
    @test isinitial(Initial(z[1]))
    @test isinitial(Initial(c[4]))
    @test !isinitial(c)
    @test !isinitial(x)
end

@testset "At" begin
    @independent_variables u
    @variables x(t) v(..) w(t)[1:3]
    @parameters y r[1:3]
    @discretes z(u, t)

    @test EvalAt(1)(x) isa Num
    @test isequal(EvalAt(1)(y), y)
    @test_throws ErrorException EvalAt(1)(z)
    @test isequal(EvalAt(1)(v), v(1))
    @test isequal(EvalAt(1)(v(t)), v(1))
    @test isequal(EvalAt(1)(v(2)), v(2))

    arr = EvalAt(1)(w)
    var = EvalAt(1)(w[1])
    @test arr isa Symbolics.Arr
    @test var isa Num

    @test isequal(EvalAt(1)(r), r)
    @test isequal(EvalAt(1)(r[2]), r[2])

    _x = ModelingToolkitBase.unwrap(x)
    @test EvalAt(1)(_x) isa Symbolics.BasicSymbolic
    @test value(only(arguments(EvalAt(1)(_x)))) == 1
    @test EvalAt(1)(D(x)) isa Num
end
