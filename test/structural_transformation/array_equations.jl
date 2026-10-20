using ModelingToolkit, Test, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEqBDF, OrdinaryDiffEqRosenbrock, SciMLBase
using DiffEqBase: BrownFullBasicInit
import ModelingToolkitTearing as MTKT
import SymbolicUtils as SU
using SymbolicUtils: unwrap

is_array_eq(eq) = SU.is_array_shape(SU.shape(unwrap(eq.lhs)))
countnodes(ex) = ex isa Expr ? 1 + sum(countnodes, ex.args; init = 0) : 1

# Finite-difference heat equation in the residual form emitted by MethodOfLines.
function heat_array_system(n)
    @variables u(t)[1:n]
    dx = 1 / (n - 1)
    lap = (u[1:(n - 2)] .- 2 .* u[2:(n - 1)] .+ u[3:n]) ./ dx^2
    interior = broadcast(-, D(u[2:(n - 1)]), lap) ~ zeros(n - 2)
    @named heat = System([interior, u[1] ~ 0.0, u[n] ~ 0.0], t, collect(u), [])
    return heat, u
end

@testset "Pure array differential equation" begin
    @variables x(t)[1:4]
    @named sys = System([D(x) ~ -x], t)
    ssys = mtkcompile(sys; preserve_array_equations = true)
    @test length(equations(ssys)) == 1
    eq = only(equations(ssys))
    @test isequal(eq.lhs, D(x))
    @test issetequal(unknowns(ssys), collect(x))
    ts = ModelingToolkit.get_tearing_state(ssys)
    @test length(ts.array_groups) == 1
    @test !only(ts.array_groups).dirty

    # Default behaviour is unchanged: everything is scalarized.
    ssys0 = mtkcompile(sys)
    @test length(equations(ssys0)) == 4
    @test !any(is_array_eq, equations(ssys0))
end

@testset "Array DE with observed and torn algebraic equations" begin
    @variables x(t)[1:5] y(t) z(t)
    @parameters k
    eqs = [D(x) ~ -k .* x .+ y, y ~ sum(x), 0 ~ z^3 + z - y]
    @named sys = System(eqs, t, [collect(x); y; z], [k])
    ssys = mtkcompile(sys; preserve_array_equations = true)
    eqs = equations(ssys)
    @test length(eqs) == 2
    @test count(is_array_eq, eqs) == 1
    @test isequal(eqs[findfirst(is_array_eq, eqs)].lhs, D(x))
    @test issetequal(unknowns(ssys), [collect(x); z])
    @test isequal(only(observed(ssys)).lhs, y)

    # `ODEProblem` from a system with array equations needs the array-aware explicit
    # codegen (SciML/ModelingToolkit.jl#5101); the implicit DAE path already supports it.
    op = [x => ones(5), z => 1.5, k => 1.0, D(x) => zeros(5), D(z) => 0.0]
    prob = DAEProblem(ssys, op, (0.0, 1.0); build_initializeprob = false)
    sol = solve(prob, DFBDF(); initializealg = BrownFullBasicInit(), reltol = 1.0e-8, abstol = 1.0e-8)
    @test SciMLBase.successful_retcode(sol)
    prob0 = DAEProblem(mtkcompile(sys), op, (0.0, 1.0); build_initializeprob = false)
    sol0 = solve(prob0, DFBDF(); initializealg = BrownFullBasicInit(), reltol = 1.0e-8, abstol = 1.0e-8)
    @test sol[x][end] ≈ sol0[x][end] rtol = 1.0e-6
    @test sol[z][end] ≈ sol0[z][end] rtol = 1.0e-6
end

@testset "Alias elimination leaves the array equation intact" begin
    @variables x(t)[1:3] a(t) b(t)
    eqs = [D(x) ~ -x .+ a, a ~ 2b, b ~ sin(t)]
    @named sys = System(eqs, t)
    ssys = mtkcompile(sys; preserve_array_equations = true)
    @test length(equations(ssys)) == 1
    @test isequal(only(equations(ssys)).lhs, D(x))
    @test issetequal(unknowns(ssys), collect(x))
    @test length(observed(ssys)) == 2
end

@testset "High index systems scalarize the array equation" begin
    # Pendulum in Cartesian coordinates with array positions/velocities.
    @variables q(t)[1:2] v(t)[1:2] T(t)
    @parameters g L
    eqs = [
        D(q) ~ v,
        D(v) ~ T .* q .+ [0, -g],
        0 ~ sum(q .^ 2) - L^2,
    ]
    @named pend = System(eqs, t)
    ssys = mtkcompile(pend; preserve_array_equations = true)
    ts = ModelingToolkit.get_tearing_state(ssys)
    # `D(v) ~ ...` is differentiated by Pantelides / dummy derivatives and cannot survive.
    @test any(g -> g.dirty, ts.array_groups)
    @test !any(is_array_eq, equations(ssys))
    ssys0 = mtkcompile(pend)
    @test length(equations(ssys)) == length(equations(ssys0))
    op = [q => [1.0, 0.0], v => [0.0, 0.0], g => 9.81, L => 1.0]
    prob = ODEProblem(ssys, op, (0.0, 1.0); guesses = [T => 0.0])
    sol = solve(prob, FBDF())
    @test SciMLBase.successful_retcode(sol)
    sol0 = solve(ODEProblem(ssys0, op, (0.0, 1.0); guesses = [T => 0.0]), FBDF())
    @test sol[q][end] ≈ sol0[q][end]
end

@testset "MethodOfLines-style heat equation" begin
    n = 21
    heat, u = heat_array_system(n)
    ssys = mtkcompile(heat; preserve_array_equations = true)
    eqs = equations(ssys)
    @test length(eqs) == 1
    @test isequal(only(eqs).lhs, D(u[2:(n - 1)]))
    @test issetequal(unknowns(ssys), collect(u)[2:(n - 1)])
    @test any(eq -> isequal(eq.lhs, u[1]), observed(ssys))
    @test any(eq -> isequal(eq.lhs, u[n]), observed(ssys))

    xs = range(0.0, 1.0, length = n)
    op = [[u[i] => sinpi(xs[i]) for i in 2:(n - 1)]; [D(u[i]) => 0.0 for i in 2:(n - 1)]]
    prob = DAEProblem(ssys, op, (0.0, 0.1); build_initializeprob = false)
    sol = solve(prob, DFBDF(); initializealg = BrownFullBasicInit(), reltol = 1.0e-8, abstol = 1.0e-8)
    @test SciMLBase.successful_retcode(sol)
    exact = [exp(-pi^2 * 0.1) * sinpi(xi) for xi in xs]
    @test maximum(abs.(sol[u][end] .- exact)) < 1.0e-3

    prob0 = DAEProblem(mtkcompile(heat), op, (0.0, 0.1); build_initializeprob = false)
    sol0 = solve(prob0, DFBDF(); initializealg = BrownFullBasicInit(), reltol = 1.0e-8, abstol = 1.0e-8)
    @test sol[u][end] ≈ sol0[u][end] atol = 1.0e-7

    # Generated code does not grow with `n`.
    sizes = map((24, 48, 96)) do nn
        s, _ = heat_array_system(nn)
        cs = mtkcompile(s; preserve_array_equations = true)
        @test length(equations(cs)) == 1
        ex = ModelingToolkit.generate_rhs(cs; expression = Val{true}, implicit_dae = true)
        countnodes.(ex)
    end
    @test allequal(sizes)
    ex = ModelingToolkit.generate_rhs(ssys; expression = Val{true}, implicit_dae = true)
    @test !occursin("array_literal", string(ex))
end

@testset "`map_variables_to_equations` with preserved array equations" begin
    @variables x(t)[1:3] y(t) z(t)
    @named sys = System([D(x) ~ -x .+ y, y ~ sum(x), 0 ~ z^3 + z - y], t)
    ssys = mtkcompile(sys; preserve_array_equations = true)
    mapping = map_variables_to_equations(ssys)
    arr_eq = only(filter(is_array_eq, equations(ssys)))
    for i in 1:3
        @test isequal(mapping[x[i]], arr_eq)
    end
    @test isequal(mapping[y].lhs, y)
    @test isequal(mapping[z], only(filter(!is_array_eq, equations(ssys))))
end
