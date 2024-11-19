using ModelingToolkit, Test
using ModelingToolkit: value, vars, parse_variable
using SymbolicUtils: <ₑ

@parameters α β δ
expr = (((1 / β - 1) + δ) / α)^(1 / (α - 1))
ref = sort([β, δ, α], lt = <ₑ)
sol = sort(Num.(ModelingToolkit.get_variables(expr)), lt = <ₑ)
@test all(x -> x isa Num, sol[i] == ref[i] for i in 1:3)
@test all(simplify ∘ value, sol[i] == ref[i] for i in 1:3)

@parameters γ
s = α => γ
expr = (((1 / β - 1) + δ) / α)^(1 / (α - 1))
sol = ModelingToolkit.substitute(expr, s)
new = (((1 / β - 1) + δ) / γ)^(1 / (γ - 1))
@test iszero(sol - new)

# Continuous
using ModelingToolkit: isdifferential, vars, collect_differential_variables,
                       collect_ivs
@independent_variables t
@variables u(t) y(t)
D = Differential(t)
eq = D(y) ~ u
v = vars(eq)
@test v == Set([D(y), u])

ov = collect_differential_variables(eq)
@test ov == Set(Any[y])

aov = ModelingToolkit.collect_applied_operators(eq, Differential)
@test aov == Set(Any[D(y)])

ts = collect_ivs([eq])
@test ts == Set([t])

@testset "vars searching through array of symbolics" begin
    fn(x, y) = sum(x) + y
    @register_symbolic fn(x::AbstractArray, y)
    @variables x y z
    res = vars(fn([x, y], z))
    @test length(res) == 3
end

@testset "parse_variable" begin
    @mtkmodel Lorenz begin
        @variables begin
            x(t)
            y(t)
            z(t)
        end
        @parameters begin
            σ
            ρ
            β
        end
        @equations begin
            D(D(x)) ~ σ * (y - x)
            D(y) ~ x * (ρ - z) - y
            D(z) ~ x * y - β * z
        end
    end
    @mtkmodel ArrSys begin
        @variables begin
            x(t)[1:2]
        end
        @parameters begin
            p[1:2, 1:2]
        end
        @equations begin
            D(D(x)) ~ p * x
        end
    end
    @mtkmodel Outer begin
        @components begin
            😄 = Lorenz()
            arr = ArrSys()
        end
    end

    @mtkbuild sys = Outer()
    for (str, var) in [
        # unicode system, scalar variable
        ("😄.x", sys.😄.x),
        ("😄.x(t)", sys.😄.x),
        ("😄₊x", sys.😄.x),
        ("😄₊x(t)", sys.😄.x),
        # derivative
        ("D(😄.x)", D(sys.😄.x)),
        ("D(😄.x(t))", D(sys.😄.x)),
        ("D(😄₊x)", D(sys.😄.x)),
        ("D(😄₊x(t))", D(sys.😄.x)),
        # other derivative
        ("😄.xˍt", D(sys.😄.x)),
        ("😄.x(t)ˍt", D(sys.😄.x)),
        ("😄₊xˍt", D(sys.😄.x)),
        ("😄₊x(t)ˍt", D(sys.😄.x)),
        # scalar parameter
        ("😄.σ", sys.😄.σ),
        ("😄₊σ", sys.😄.σ),
        # array variable
        ("arr.x", sys.arr.x),
        ("arr₊x", sys.arr.x),
        ("arr.x(t)", sys.arr.x),
        ("arr₊x(t)", sys.arr.x),
        # getindex
        ("arr.x[1]", sys.arr.x[1]),
        ("arr₊x[1]", sys.arr.x[1]),
        ("arr.x(t)[1]", sys.arr.x[1]),
        ("arr₊x(t)[1]", sys.arr.x[1]),
        # derivative
        ("D(arr.x(t))", D(sys.arr.x)),
        ("D(arr₊x(t))", D(sys.arr.x)),
        ("D(arr.x[1])", D(sys.arr.x[1])),
        ("D(arr₊x[1])", D(sys.arr.x[1])),
        ("D(arr.x(t)[1])", D(sys.arr.x[1])),
        ("D(arr₊x(t)[1])", D(sys.arr.x[1])),
        # other derivative
        ("arr.xˍt", D(sys.arr.x)),
        ("arr₊xˍt", D(sys.arr.x)),
        ("arr.xˍt(t)", D(sys.arr.x)),
        ("arr₊xˍt(t)", D(sys.arr.x)),
        ("arr.xˍt[1]", D(sys.arr.x[1])),
        ("arr₊xˍt[1]", D(sys.arr.x[1])),
        ("arr.xˍt(t)[1]", D(sys.arr.x[1])),
        ("arr₊xˍt(t)[1]", D(sys.arr.x[1])),
        ("arr.x(t)ˍt", D(sys.arr.x)),
        ("arr₊x(t)ˍt", D(sys.arr.x)),
        ("arr.x(t)ˍt[1]", D(sys.arr.x[1])),
        ("arr₊x(t)ˍt[1]", D(sys.arr.x[1])),
        # array parameter
        ("arr.p", sys.arr.p),
        ("arr₊p", sys.arr.p),
        ("arr.p[1, 2]", sys.arr.p[1, 2]),
        ("arr₊p[1, 2]", sys.arr.p[1, 2])
    ]
        isequal(parse_variable(sys, str), var)
    end
end
