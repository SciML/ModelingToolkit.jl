using Test
using ModelingToolkit
using Graphs
using SparseArrays
using UnPack
using ModelingToolkit: t_nounits as t, D_nounits as D, default_toterm
using Symbolics: unwrap
using DataInterpolations
using OrdinaryDiffEq, NonlinearSolve, StochasticDiffEq
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
pendulum = System(eqs, t, [x, y, w, z, T], [L, g], name = :pendulum)
state = TearingState(pendulum)
StructuralTransformations.find_solvables!(state)
sss = state.structure
@unpack graph, solvable_graph, var_to_diff = sss
@test sort(graph.fadjlist) == [[1, 7], [2, 8], [3, 5, 9], [4, 6, 9], [5, 6]]
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
    @mtkcompile sys = System(
        [D(x) ~ z[1] + z[2] + foo(z)[1], y[1] ~ 2t, y[2] ~ 3t, z ~ foo(y)], t)
    @test length(equations(sys)) == 1
    @test length(observed(sys)) == 6
    @test any(obs -> isequal(obs, y), observables(sys))
    @test any(obs -> isequal(obs, z), observables(sys))
    prob = ODEProblem(sys, [x => 1.0, foo => _tmp_fn], (0.0, 1.0))
    @test_nowarn prob.f(prob.u0, prob.p, 0.0)

    isys = ModelingToolkit.generate_initializesystem(sys)
    @test length(unknowns(isys)) == 5
    @test length(equations(isys)) == 4
    @test !any(equations(isys)) do eq
        iscall(eq.rhs) && operation(eq.rhs) in [StructuralTransformations.change_origin]
    end
end

@testset "array hack can be disabled" begin
    @testset "fully_determined = true" begin
        @variables x(t) y(t)[1:2] z(t)[1:2]
        @parameters foo(::AbstractVector)[1:2]
        _tmp_fn(x) = 2x
        @named sys = System(
            [D(x) ~ z[1] + z[2] + foo(z)[1], y[1] ~ 2t, y[2] ~ 3t, z ~ foo(y)], t)

        sys2 = mtkcompile(sys; array_hack = false)
        @test length(observed(sys2)) == 4
        @test !any(observed(sys2)) do eq
            iscall(eq.rhs) && operation(eq.rhs) == StructuralTransformations.change_origin
        end
    end

    @testset "fully_determined = false" begin
        @variables x(t) y(t)[1:2] z(t)[1:2] w(t)
        @parameters foo(::AbstractVector)[1:2]
        _tmp_fn(x) = 2x
        @named sys = System(
            [D(x) ~ z[1] + z[2] + foo(z)[1] + w, y[1] ~ 2t, y[2] ~ 3t, z ~ foo(y)], t)

        sys2 = mtkcompile(sys; array_hack = false, fully_determined = false)
        @test length(observed(sys2)) == 4
        @test !any(observed(sys2)) do eq
            iscall(eq.rhs) && operation(eq.rhs) == StructuralTransformations.change_origin
        end
    end
end

@testset "additional passes" begin
    @variables x(t) y(t)
    @named sys = System([D(x) ~ x, y ~ x + t], t)
    value = Ref(0)
    pass(sys; kwargs...) = (value[] += 1; return sys)
    mtkcompile(sys; additional_passes = [pass])
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
    @testset "Requires simplified system" begin
        @variables x(t) y(t)
        @named sys = System([D(x) ~ x, y ~ 2x], t)
        sys = complete(sys)
        @test_throws ArgumentError map_variables_to_equations(sys)
    end
    @testset "`ODESystem`" begin
        @variables x(t) y(t) z(t)
        @mtkcompile sys = System([D(x) ~ 2x + y, y ~ x + z, z^3 + x^3 ~ 12], t)
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
            @mtkcompile sys = System(eqs, t)
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
            @mtkcompile sys = compose(coupledOsc, systems)
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
        @mtkcompile sys = System([x^2 ~ 2y^2 + 1, sin(z) ~ y, z^3 + 4z + 1 ~ 0])
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
        return System(eqs, t, vars, params; systems, name)
    end

    @component function FilteredInputExplicit(; name, x0 = 0, T = 0.1)
        params = @parameters begin
            k(t)[1:1] = [x0]
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
               D(k[1]) ~ 1.0
               dx ~ (k[1] - x) / T]
        return System(eqs, t, vars, params; systems, name)
    end

    @component function FilteredInputErr(; name, x0 = 0, T = 0.1)
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
               D(k) ~ missing]
        return System(eqs, t, vars, params; systems, name)
    end

    @named sys = FilteredInputErr()
    @test_throws ["derivative of discrete variable", "k(t)"] mtkcompile(sys)

    @mtkcompile sys = FilteredInput()
    vs = Set()
    for eq in equations(sys)
        ModelingToolkit.vars!(vs, eq)
    end
    for eq in observed(sys)
        ModelingToolkit.vars!(vs, eq)
    end

    @test !(D(sys.k) in vs)

    @mtkcompile sys = FilteredInputExplicit()
    obsfn1 = ModelingToolkit.build_explicit_observed_function(sys, sys.ddx)
    obsfn2 = ModelingToolkit.build_explicit_observed_function(sys, sys.dx)
    u = [1.0]
    p = MTKParameters(sys, [sys.k => [2.0], sys.T => 3.0])
    @test obsfn1(u, p, 0.0) ≈ (1 - obsfn2(u, p, 0.0)) / 3.0

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
            return System(eqs, t, vars, params; systems, name)
        end

        @mtkcompile sys = FilteredInput2()
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

@testset "Don't rely on metadata" begin
    @testset "ODESystem" begin
        @variables x(t) p
        @parameters y(t) q
        @mtkcompile sys = System([D(x) ~ x * q, x^2 + y^2 ~ p], t, [x, y],
            [p, q]; initialization_eqs = [p + q ~ 3],
            defaults = [p => missing], guesses = [p => 1.0, y => 1.0])
        @test length(equations(sys)) == 2
        @test length(parameters(sys)) == 2
        prob = ODEProblem(sys, [x => 1.0, q => 2.0], (0.0, 1.0))
        integ = init(prob, Rodas5P(); abstol = 1e-10, reltol = 1e-8)
        @test integ.ps[p]≈1.0 atol=1e-6
        @test integ[y]≈0.0 atol=1e-5
    end

    @testset "NonlinearSystem" begin
        @variables x p
        @parameters y q
        @mtkcompile sys = System([0 ~ p * x + y, x^3 + y^3 ~ q], [x, y],
            [p, q]; initialization_eqs = [p ~ q + 1],
            guesses = [p => 1.0], defaults = [p => missing])
        @test length(equations(sys)) == length(unknowns(sys)) == 1
        @test length(observed(sys)) == 1
        @test observed(sys)[1].lhs in Set([x, y])
        @test length(parameters(sys)) == 2
        prob = NonlinearProblem(sys, [x => 1.0, y => 1.0, q => 1.0])
        integ = init(prob, NewtonRaphson())
        @test prob.ps[p] ≈ 2.0
    end

    @testset "SDESystem" begin
        @variables x(t) p a
        @parameters y(t) q b
        @brownians c
        @mtkcompile sys = System([D(x) ~ x + q * a, D(y) ~ y + p * b + c], t, [x, y],
            [p, q], [a, b, c]; initialization_eqs = [p + q ~ 4],
            guesses = [p => 1.0], defaults = [p => missing])
        @test length(equations(sys)) == 2
        @test issetequal(unknowns(sys), [x, y])
        @test issetequal(parameters(sys), [p, q])
        @test isempty(brownians(sys))
        neqs = ModelingToolkit.get_noise_eqs(sys)
        @test issetequal(sum.(eachrow(neqs)), [q, 1 + p])
        prob = SDEProblem(sys, [x => 1.0, y => 1.0, q => 1.0], (0.0, 1.0))
        integ = init(prob, ImplicitEM())
        @test integ.ps[p] ≈ 3.0
    end
end

@testset "Deprecated `mtkcompile` and `@mtkcompile`" begin
    @variables x(t)
    @test_deprecated @mtkbuild sys = System([D(x) ~ x], t)
    @named sys = System([D(x) ~ x], t)
    @test_deprecated structural_simplify(sys)
end
