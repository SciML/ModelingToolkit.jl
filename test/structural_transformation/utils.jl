using Test
using ModelingToolkit
using Graphs
using SparseArrays
using UnPack
using ModelingToolkit: t_nounits as t, D_nounits as D, default_toterm
using Symbolics: unwrap
using DataInterpolations
const ST = StructuralTransformations

# Define some variables
@parameters L g
@variables x(t) y(t) w(t) z(t) T(t)

# Simple pendulum in cartesian coordinates
eqs = [D(x) ~ w,
    D(y) ~ z,
    D(w) ~ T * x,
    D(z) ~ T * y - g,
    0 ~ x^2 + y^2 - L^2]
pendulum = ODESystem(eqs, t, [x, y, w, z, T], [L, g], name = :pendulum)
state = TearingState(pendulum)
StructuralTransformations.find_solvables!(state)
sss = state.structure
@unpack graph, solvable_graph, var_to_diff = sss
@test graph.fadjlist == [[1, 7], [2, 8], [3, 5, 9], [4, 6, 9], [5, 6]]
@test length(graph.badjlist) == 9
@test ne(graph) == nnz(incidence_matrix(graph)) == 12
@test nv(solvable_graph) == 9 + 5
let N = nothing
    @test var_to_diff == [N, N, N, N, 1, 2, 3, 4, N]
end

se = collect(StructuralTransformations.edges(graph))
@test se == mapreduce(vcat, enumerate(graph.fadjlist)) do (s, d)
    StructuralTransformations.BipartiteEdge.(s, d)
end

@testset "observed2graph handles unknowns inside callable parameters" begin
    @variables x(t) y(t)
    @parameters p(..)
    g, _ = ModelingToolkit.observed2graph([y ~ p(x), x ~ 0], [y, x])
    @test ModelingToolkit.𝑠neighbors(g, 1) == [2]
    @test ModelingToolkit.𝑑neighbors(g, 2) == [1]
end

@testset "array observed used unscalarized in another observed" begin
    @variables x(t) y(t)[1:2] z(t)[1:2]
    @parameters foo(::AbstractVector)[1:2]
    _tmp_fn(x) = 2x
    @mtkbuild sys = ODESystem(
        [D(x) ~ z[1] + z[2] + foo(z)[1], y[1] ~ 2t, y[2] ~ 3t, z ~ foo(y)], t)
    @test length(equations(sys)) == 1
    @test length(observed(sys)) == 7
    @test any(obs -> isequal(obs, y), observables(sys))
    @test any(obs -> isequal(obs, z), observables(sys))
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0), [foo => _tmp_fn])
    @test_nowarn prob.f(prob.u0, prob.p, 0.0)

    isys = ModelingToolkit.generate_initializesystem(sys)
    @test length(unknowns(isys)) == 5
    @test length(equations(isys)) == 4
    @test !any(equations(isys)) do eq
        iscall(eq.rhs) && operation(eq.rhs) in [StructuralTransformations.getindex_wrapper,
            StructuralTransformations.change_origin]
    end
end

@testset "scalarized array observed calling same function multiple times" begin
    @variables x(t) y(t)[1:2]
    @parameters foo(::Real)[1:2]
    val = Ref(0)
    function _tmp_fn2(x)
        val[] += 1
        return [x, 2x]
    end
    @mtkbuild sys = ODESystem([D(x) ~ y[1] + y[2], y ~ foo(x)], t)
    @test length(equations(sys)) == 1
    @test length(observed(sys)) == 3
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0), [foo => _tmp_fn2])
    val[] = 0
    @test_nowarn prob.f(prob.u0, prob.p, 0.0)
    @test val[] == 1

    isys = ModelingToolkit.generate_initializesystem(sys)
    @test length(unknowns(isys)) == 3
    @test length(equations(isys)) == 2
    @test !any(equations(isys)) do eq
        iscall(eq.rhs) && operation(eq.rhs) in [StructuralTransformations.getindex_wrapper,
            StructuralTransformations.change_origin]
    end

    @testset "CSE hack in equations(sys)" begin
        val[] = 0
        @variables z(t)[1:2]
        @mtkbuild sys = ODESystem(
            [D(y) ~ foo(x), D(x) ~ sum(y), zeros(2) ~ foo(prod(z))], t)
        @test length(equations(sys)) == 5
        @test length(observed(sys)) == 2
        prob = ODEProblem(
            sys, [y => ones(2), z => 2ones(2), x => 3.0], (0.0, 1.0), [foo => _tmp_fn2])
        val[] = 0
        @test_nowarn prob.f(prob.u0, prob.p, 0.0)
        @test val[] == 2

        isys = ModelingToolkit.generate_initializesystem(sys)
        @test length(unknowns(isys)) == 5
        @test length(equations(isys)) == 2
        @test !any(equations(isys)) do eq
            iscall(eq.rhs) &&
                operation(eq.rhs) in [StructuralTransformations.getindex_wrapper,
                    StructuralTransformations.change_origin]
        end
    end
end

@testset "array and cse hacks can be disabled" begin
    @testset "fully_determined = true" begin
        @variables x(t) y(t)[1:2] z(t)[1:2]
        @parameters foo(::AbstractVector)[1:2]
        _tmp_fn(x) = 2x
        @named sys = ODESystem(
            [D(x) ~ z[1] + z[2] + foo(z)[1], y[1] ~ 2t, y[2] ~ 3t, z ~ foo(y)], t)

        sys1 = structural_simplify(sys; cse_hack = false)
        @test length(observed(sys1)) == 6
        @test !any(observed(sys1)) do eq
            iscall(eq.rhs) &&
                operation(eq.rhs) == StructuralTransformations.getindex_wrapper
        end

        sys2 = structural_simplify(sys; array_hack = false)
        @test length(observed(sys2)) == 5
        @test !any(observed(sys2)) do eq
            iscall(eq.rhs) && operation(eq.rhs) == StructuralTransformations.change_origin
        end
    end

    @testset "fully_determined = false" begin
        @variables x(t) y(t)[1:2] z(t)[1:2] w(t)
        @parameters foo(::AbstractVector)[1:2]
        _tmp_fn(x) = 2x
        @named sys = ODESystem(
            [D(x) ~ z[1] + z[2] + foo(z)[1] + w, y[1] ~ 2t, y[2] ~ 3t, z ~ foo(y)], t)

        sys1 = structural_simplify(sys; cse_hack = false, fully_determined = false)
        @test length(observed(sys1)) == 6
        @test !any(observed(sys1)) do eq
            iscall(eq.rhs) &&
                operation(eq.rhs) == StructuralTransformations.getindex_wrapper
        end

        sys2 = structural_simplify(sys; array_hack = false, fully_determined = false)
        @test length(observed(sys2)) == 5
        @test !any(observed(sys2)) do eq
            iscall(eq.rhs) && operation(eq.rhs) == StructuralTransformations.change_origin
        end
    end
end

@testset "additional passes" begin
    @variables x(t) y(t)
    @named sys = ODESystem([D(x) ~ x, y ~ x + t], t)
    value = Ref(0)
    pass(sys; kwargs...) = (value[] += 1; return sys)
    structural_simplify(sys; additional_passes = [pass])
    @test value[] == 1
end

@testset "Distribute shifts" begin
    @variables x(t) y(t) z(t)
    @parameters a b c
    k = ShiftIndex(t)

    # Expand shifts
    @test isequal(
        ST.distribute_shift(Shift(t, -1)(x + y)), Shift(t, -1)(x) + Shift(t, -1)(y))

    expr = a * Shift(t, -2)(x) + Shift(t, 2)(y) + b
    @test isequal(ST.simplify_shifts(ST.distribute_shift(Shift(t, 2)(expr))),
        a * x + Shift(t, 4)(y) + b)
    @test isequal(ST.distribute_shift(Shift(t, 2)(exp(z))), exp(Shift(t, 2)(z)))
    @test isequal(ST.distribute_shift(Shift(t, 2)(exp(a) + b)), exp(a) + b)

    expr = a^x - log(b * y) + z * x
    @test isequal(ST.distribute_shift(Shift(t, -3)(expr)),
        a^(Shift(t, -3)(x)) - log(b * Shift(t, -3)(y)) + Shift(t, -3)(z) * Shift(t, -3)(x))

    expr = x(k + 1) ~ x + x(k - 1)
    @test isequal(ST.distribute_shift(Shift(t, -1)(expr)), x ~ x(k - 1) + x(k - 2))
end

@testset "`map_variables_to_equations`" begin
    @testset "Not supported for systems without `.tearing_state`" begin
        @variables x
        @mtkbuild sys = OptimizationSystem(x^2)
        @test_throws ArgumentError map_variables_to_equations(sys)
    end
    @testset "Requires simplified system" begin
        @variables x(t) y(t)
        @named sys = ODESystem([D(x) ~ x, y ~ 2x], t)
        sys = complete(sys)
        @test_throws ArgumentError map_variables_to_equations(sys)
    end
    @testset "`ODESystem`" begin
        @variables x(t) y(t) z(t)
        @mtkbuild sys = ODESystem([D(x) ~ 2x + y, y ~ x + z, z^3 + x^3 ~ 12], t)
        mapping = map_variables_to_equations(sys)
        @test mapping[x] == (D(x) ~ 2x + y)
        @test mapping[y] == (y ~ x + z)
        @test mapping[z] == (0 ~ 12 - z^3 - x^3)
        @test length(mapping) == 3

        @testset "With dummy derivatives" begin
            @parameters g
            @variables x(t) y(t) [state_priority = 10] λ(t)
            eqs = [D(D(x)) ~ λ * x
                   D(D(y)) ~ λ * y - g
                   x^2 + y^2 ~ 1]
            @mtkbuild sys = ODESystem(eqs, t)
            mapping = map_variables_to_equations(sys)

            yt = default_toterm(unwrap(D(y)))
            xt = default_toterm(unwrap(D(x)))
            xtt = default_toterm(unwrap(D(D(x))))
            @test mapping[x] == (0 ~ 1 - x^2 - y^2)
            @test mapping[y] == (D(y) ~ yt)
            @test mapping[D(y)] == (D(yt) ~ -g + y * λ)
            @test mapping[D(x)] == (0 ~ -2xt * x - 2yt * y)
            @test mapping[D(D(x))] == (xtt ~ x * λ)
            @test length(mapping) == 5

            @testset "`rename_dummy_derivatives = false`" begin
                mapping = map_variables_to_equations(sys; rename_dummy_derivatives = false)

                @test mapping[x] == (0 ~ 1 - x^2 - y^2)
                @test mapping[y] == (D(y) ~ yt)
                @test mapping[yt] == (D(yt) ~ -g + y * λ)
                @test mapping[xt] == (0 ~ -2xt * x - 2yt * y)
                @test mapping[xtt] == (xtt ~ x * λ)
                @test length(mapping) == 5
            end
        end
        @testset "DDEs" begin
            function oscillator(; name, k = 1.0, τ = 0.01)
                @parameters k=k τ=τ
                @variables x(..)=0.1 y(t)=0.1 jcn(t)=0.0 delx(t)
                eqs = [D(x(t)) ~ y,
                    D(y) ~ -k * x(t - τ) + jcn,
                    delx ~ x(t - τ)]
                return System(eqs, t; name = name)
            end

            systems = @named begin
                osc1 = oscillator(k = 1.0, τ = 0.01)
                osc2 = oscillator(k = 2.0, τ = 0.04)
            end
            eqs = [osc1.jcn ~ osc2.delx,
                osc2.jcn ~ osc1.delx]
            @named coupledOsc = System(eqs, t)
            @mtkbuild sys = compose(coupledOsc, systems)
            mapping = map_variables_to_equations(sys)
            x1 = operation(unwrap(osc1.x))
            x2 = operation(unwrap(osc2.x))
            @test mapping[osc1.x] == (D(osc1.x) ~ osc1.y)
            @test mapping[osc1.y] == (D(osc1.y) ~ osc1.jcn - osc1.k * x1(t - osc1.τ))
            @test mapping[osc1.delx] == (osc1.delx ~ x1(t - osc1.τ))
            @test mapping[osc1.jcn] == (osc1.jcn ~ osc2.delx)
            @test mapping[osc2.x] == (D(osc2.x) ~ osc2.y)
            @test mapping[osc2.y] == (D(osc2.y) ~ osc2.jcn - osc2.k * x2(t - osc2.τ))
            @test mapping[osc2.delx] == (osc2.delx ~ x2(t - osc2.τ))
            @test mapping[osc2.jcn] == (osc2.jcn ~ osc1.delx)
            @test length(mapping) == 8
        end
    end
    @testset "`NonlinearSystem`" begin
        @variables x y z
        @mtkbuild sys = NonlinearSystem([x^2 ~ 2y^2 + 1, sin(z) ~ y, z^3 + 4z + 1 ~ 0])
        mapping = map_variables_to_equations(sys)
        @test mapping[x] == (0 ~ 2y^2 + 1 - x^2)
        @test mapping[y] == (y ~ sin(z))
        @test mapping[z] == (0 ~ -1 - 4z - z^3)
        @test length(mapping) == 3
    end
end

@testset "Issue#3480: Derivatives of time-dependent parameters" begin
    @component function FilteredInput(; name, x0 = 0, T = 0.1)
        params = @parameters begin
            k(t) = x0
            T = T
        end
        vars = @variables begin
            x(t) = k
            dx(t) = 0
            ddx(t)
        end
        systems = []
        eqs = [D(x) ~ dx
               D(dx) ~ ddx
               dx ~ (k - x) / T]
        return ODESystem(eqs, t, vars, params; systems, name)
    end

    @component function FilteredInputFix(; name, x0 = 0, T = 0.1)
        params = @parameters begin
            k(t) = x0
            T = T
        end
        vars = @variables begin
            x(t) = k
            dx(t) = 0
            ddx(t)
        end
        systems = []
        eqs = [D(x) ~ dx
               D(dx) ~ ddx
               dx ~ (k - x) / T
               D(k) ~ 0]
        return ODESystem(eqs, t, vars, params; systems, name)
    end

    @named sys = FilteredInput()
    @test_throws ["derivative of discrete variable", "k(t)"] structural_simplify(sys)

    @mtkbuild sys = FilteredInputFix()
    vs = Set()
    for eq in equations(sys)
        ModelingToolkit.vars!(vs, eq)
    end
    for eq in observed(sys)
        ModelingToolkit.vars!(vs, eq)
    end

    @test !(D(sys.k) in vs)

    @testset "Called parameter still has derivative" begin
        @component function FilteredInput2(; name, x0 = 0, T = 0.1)
            ts = collect(0.0:0.1:10.0)
            spline = LinearInterpolation(ts .^ 2, ts)
            params = @parameters begin
                (k::LinearInterpolation)(..) = spline
                T = T
            end
            vars = @variables begin
                x(t) = k(t)
                dx(t) = 0
                ddx(t)
            end
            systems = []
            eqs = [D(x) ~ dx
                   D(dx) ~ ddx
                   dx ~ (k(t) - x) / T]
            return ODESystem(eqs, t, vars, params; systems, name)
        end

        @mtkbuild sys = FilteredInput2()
        vs = Set()
        for eq in equations(sys)
            ModelingToolkit.vars!(vs, eq)
        end
        for eq in observed(sys)
            ModelingToolkit.vars!(vs, eq)
        end

        @test D(sys.k(t)) in vs
    end
end
