using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D
import SymbolicUtils as SU
using Symbolics: SConst

@testset "`state_priorities` and `irreducibles` kwargs take priority over metadata" begin
    @variables x(t) [state_priority = 3] y(t) [irreducible = false]
    @named sys = System([D(x) ~ 2x + y], t; state_priorities = [x => -4], irreducibles = [y])
    @test state_priorities(sys)[x] == -4
    @test y in irreducibles(sys)
end

@testset "Issue#4138: expressions of `Initial` allowed as parameter defaults" begin
    @variables X1(t) X2(t)
    @parameters k1 k2 Γ = Initial(X1) + Initial(X2)
    eqs = [
        D(X1) ~ -k1 * X1 + k2 * (-X1 + Γ),
        X2 ~ -X1 + Γ,
    ]
    @mtkcompile sys = System(eqs, t)
    @test isequal(bindings(sys)[Γ], Initial(X1) + Initial(X2))
end

@testset "Variable discovery from `initial_conditions` and `bindings` kwargs (time-dependent)" begin
    @variables x(t) y(t)   # y is NOT in any equation
    @parameters p q         # q is NOT in any equation

    # Only x and p appear in the equations
    eqs = [D(x) ~ p * x]

    # y appears only as the RHS of an initial_conditions pair
    # q appears only as the RHS of a bindings pair
    # p appears only as the LHS of a bindings pair (already in eqs, but verifying no breakage)
    @named sys = System(eqs, t; initial_conditions = [x => y], bindings = [p => q])

    @test y ∈ Set(unknowns(sys))
    @test q ∈ Set(parameters(sys))
end

@testset "Variable discovery from `initial_conditions` and `bindings` kwargs (time-independent)" begin
    @variables x y          # y is NOT in any equation
    @parameters p q         # q is NOT in any equation

    # Only x and p appear in the equations
    eqs = [0 ~ p * x + 1]

    @named sys = System(eqs; initial_conditions = [x => y], bindings = [p => q])

    @test y ∈ Set(unknowns(sys))
    @test q ∈ Set(parameters(sys))
end

@testset "`Adjoint` in equations" begin
    @variables x(t)[1:3, 1:3]
    @named sys = System([D(x) ~ x'], t)
    @test_nowarn mtkcompile(sys)
end

@testset "Necessary initial conditions" begin
    @variables x(t) y(t)
    @mtkcomplete sys = System([D(x) ~ t, x ~ y], t)
    ics = ModelingToolkitBase.get_necessary_initial_conditions(sys)
    ics[y] = "Because I want to"
    sys = ModelingToolkitBase.set_necessary_initial_conditions(sys, ics)
    @test_throws ModelingToolkitBase.MissingNecessaryInitialConditionsError ODEProblem(sys, [x => 1], (0.0, 1.0))
end

struct MyOp <: SU.Operator end
ModelingToolkitBase.validate_operator(::MyOp, args...; kws...) = nothing

@testset "Const argument of operator is not considered a variable" begin
    @variables x(t)
    @named sys = System([D(x) ~ x + SU.term(MyOp(), [0.0]; type = Real, shape = SU.ShapeVecT())], t)
    @test issetequal(unknowns(sys), [x])
end

@testset "Variable discovery `x[i]` where `i` is a variable" begin
    @variables x(t)[1:2] i(t)::Int
    @named sys = System([x[i] ~ 0], t)
    @test issetequal(unknowns(sys), [x, i])
end

@testset "`set_defaults`" begin
    @testset "Literal values become initial conditions, symbolic values become bindings" begin
        @variables x(t) y(t)
        @parameters p q
        @named sys = System([D(x) ~ p * x, D(y) ~ q * y], t)

        sys2 = set_defaults(sys, [:x => 1.0, :p => 2.0, :y => x, :q => p])
        @test isequal(initial_conditions(sys2), Dict(x => SConst(1.0), p => SConst(2.0)))
        @test isequal(bindings(sys2), Dict(y => x, q => p))
        # `sys` is untouched
        @test isempty(initial_conditions(sys))
        @test isempty(bindings(sys))
    end

    @testset "Symbolic variables work directly as keys, not just Symbol/String names" begin
        @variables x(t) y(t)
        @parameters p q
        @named sys = System([D(x) ~ p * x, D(y) ~ q * y], t)
        sys2 = set_defaults(sys, [x => 1.0, p => y])
        @test isequal(initial_conditions(sys2), Dict(x => SConst(1.0)))
        @test isequal(bindings(sys2), Dict(p => y))
    end

    @testset "Overriding pre-existing metadata defaults, including flipping IC <-> binding" begin
        @variables x(t) = 1 y(t) = x
        @parameters p = 1 q = p
        @named sys = System([D(x) ~ p * x, D(y) ~ q * y], t)
        @test isequal(bindings(sys), Dict(y => x, q => p))
        @test isequal(initial_conditions(sys), Dict(x => SConst(1), p => SConst(1)))

        sys2 = set_defaults(sys, Dict(:x => y, :y => 2.0))
        @test isequal(bindings(sys2), Dict(x => y, q => p))
        @test isequal(initial_conditions(sys2), Dict(y => SConst(2.0), p => SConst(1)))
    end

    @testset "`nothing` clears a default, `missing` sets a solvable-parameter binding" begin
        @variables x(t) = 1
        @parameters p = 1
        @named sys = System(D(x) ~ p * x, t)
        sys2 = set_defaults(sys, [:x => nothing, :p => missing])
        @test isempty(initial_conditions(sys2))
        @test isequal(bindings(sys2), Dict(p => ModelingToolkitBase.COMMON_MISSING))
    end

    @testset "`nothing`/`missing` are recognized in their already-symbolic (wrapped) form too" begin
        @variables x(t) = 1 y(t)
        @parameters p = 1
        @named sys = System([D(x) ~ -x, D(y) ~ p], t)
        # Mixing `nothing`/`missing` with a symbolic value in a homogeneous collection can
        # promote them to `Num`-wrapped `COMMON_NOTHING`/`COMMON_MISSING` nodes instead of
        # the raw singletons -- both forms must classify identically.
        v = [:x => nothing, :y => y, :p => missing]
        sys2 = set_defaults(sys, v)
        @test !haskey(initial_conditions(sys2), x)
        @test !haskey(bindings(sys2), x)
        @test bindings(sys2)[p] === ModelingToolkitBase.COMMON_MISSING
    end

    @testset "Works with a keyword-argument pairs iterator" begin
        @variables x(t)
        @parameters p
        @named sys = System(D(x) ~ p * x, t)
        f(sys; kwargs...) = set_defaults(sys, kwargs)
        sys2 = f(sys; x = 1.0, p = 2.0)
        @test isequal(initial_conditions(sys2), Dict(x => SConst(1.0), p => SConst(2.0)))
    end

    @testset "String names, and an error when a name cannot be resolved" begin
        @variables x(t)
        @parameters p
        @named sys = System(D(x) ~ p * x, t)
        sys2 = set_defaults(sys, ["x" => 1.0])
        @test isequal(initial_conditions(sys2), Dict(x => SConst(1.0)))

        @test_throws ArgumentError set_defaults(sys, [:nonexistent => 1.0])
    end

    @component function SubComp(; k, name)
        @parameters k = k
        return System(Equation[], t, [], [k]; name)
    end

    @testset "Namespaced names resolve variables owned by a subsystem" begin
        @named sub = SubComp(; k = 1.0)
        @named sys = System(Equation[], t; systems = [sub])
        sys2 = set_defaults(sys, ["sub₊k" => 2.0])
        nnsys = toggle_namespacing(sys2, false)
        @test isequal(initial_conditions(sys2), Dict(nnsys.sub.k => SConst(2.0)))
    end

    @testset "A symbolic key from a completed, nested system resolves to itself" begin
        @named sub = SubComp(; k = 1.0)
        @named sys = System(Equation[], t; systems = [sub])
        csys = complete(sys)
        sys2 = set_defaults(csys, [csys.sub.k => 2.0])
        @test isequal(initial_conditions(sys2)[csys.sub.k], SConst(2.0))
    end

    @testset "A Vector eltype promoted to `Num` is still classified by its actual value" begin
        @variables x(t) y(t)
        @parameters p q
        @named sys = System([D(x) ~ p * x, D(y) ~ q * y], t)
        # This Vector literal has a common eltype of `Pair{Symbol, Num}`, silently
        # promoting `2.0` into `Num(2.0)` -- it must still be classified as an initial
        # condition, not a binding, despite arriving with static type `Num`.
        v = [:x => y, :y => 2.0]
        @test eltype(v) <: Pair{Symbol, Num}
        sys2 = set_defaults(sys, v)
        @test isequal(bindings(sys2), Dict(x => y))
        @test isequal(initial_conditions(sys2), Dict(y => SConst(2.0)))
    end

    @testset "Array-valued variables" begin
        @variables x(t)[1:2] y(t)[1:2]
        @parameters p[1:2]
        @named sys = System(Equation[D(x) ~ p .* x; D(y) ~ x], t, [x, y], [p])
        sys2 = set_defaults(sys, [:x => [1.0, 2.0], :p => [3.0, 4.0]])
        @test isequal(
            initial_conditions(sys2), Dict(x => SConst([1.0, 2.0]), p => SConst([3.0, 4.0]))
        )
    end

    @testset "Empty pairs is a no-op, returning the same system" begin
        @variables x(t) = 1
        @named sys = System(D(x) ~ -x, t)
        sys2 = set_defaults(sys, Pair{Symbol, Any}[])
        @test sys2 === sys
    end

    @testset "Flattened hierarchical system" begin
        @variables x(t)
        @named inner = System([D(x) ~ x], t)
        @named outer = System([D(x) ~ t], t; systems = [inner])
        outer2 = set_defaults(outer, ["inner₊x" => 1.0])
        @test SU.unwrap_const(initial_conditions(outer2)[inner.x])::Float64 == 1.0
    end
end
