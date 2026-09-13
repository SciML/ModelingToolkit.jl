using ModelingToolkitBase, SciMLBase, SymbolicIndexingInterface, Test
using OrdinaryDiffEqBDF: DFBDF

@testset "DAE array unknowns" begin
    @parameters t
    D = Differential(t)
    for dims in ((3,), (2, 3))
        @testset "shape $dims" begin
            ax = map(Base.OneTo, dims)
            @variables x(t)[ax...] z(t) y(t)[1:2]
            values = reshape(collect(1.0:prod(dims)), dims)
            yvalues = [2.0, 3.0]
            eqs = [D(x) ~ -x, z ~ sum(x), D(y) ~ -2y]
            @named raw = System(eqs, t, [x, z, y], [])
            sys = complete(raw)
            op = [x => values, z => sum(values), y => yvalues, D(z) => 0.0]
            prob = DAEProblem(sys, op, (0.0, 0.1))
            expected_u0 = vcat(vec(values), sum(values), yvalues)
            for system in (raw, sys)
                @test is_variable(system, 1)
                @test is_variable(system, length(expected_u0))
                @test !is_variable(system, 0)
                @test !is_variable(system, length(expected_u0) + 1)
            end
            @test length(unknowns(sys)) == 3
            @test prob.u0 == expected_u0
            @test prob.du0 == vcat(-vec(values), 0.0, -2yvalues)
            @test prob.differential_vars == vcat(trues(length(values)), false, trues(2))
            @test getu(sys, x)(prob) == values
            @test getu(sys, z)(prob) == sum(values)
            @test getu(sys, y)(prob) == yvalues
            residual = similar(prob.u0)
            prob.f(residual, prob.du0, prob.u0, prob.p, 0.0)
            @test residual == zeros(length(expected_u0))
            sol = solve(prob, DFBDF(); abstol = 1.0e-9, reltol = 1.0e-9)
            @test SciMLBase.successful_retcode(sol)
            expected_final = vcat(exp(-0.1) .* vec(values), exp(-0.1) * sum(values), exp(-0.2) .* yvalues)
            @test sol.u[end] ≈ expected_final rtol = 1.0e-6
            @test sol[x][end] ≈ exp(-0.1) .* values rtol = 1.0e-6
            remade = remake(prob; u0 = [x => 2values, z => 2sum(values), y => 2yvalues])
            @test remade.u0 == 2expected_u0
        end
    end
end

@testset "Initialize algebraic array unknowns" begin
    @parameters t
    D = Differential(t)
    for dims in ((3,), (2, 2))
        @testset "shape $dims" begin
            ax = map(Base.OneTo, dims)
            @variables x(t)[ax...] z(t)[ax...]
            values = reshape(collect(1.0:prod(dims)), dims)
            @named raw = System([D(x) ~ -x, z ~ x .^ 2], t, [x, z], [])
            sys = complete(raw)
            prob = DAEProblem(sys, [x => values, D(z) => zeros(dims)], (0.0, 0.1); guesses = [z => zeros(dims)])
            sol = solve(prob, DFBDF(); abstol = 1.0e-9, reltol = 1.0e-9)
            @test SciMLBase.successful_retcode(sol)
            @test sol[z][1] ≈ values .^ 2
            @test sol[z][end] ≈ exp(-0.2) .* values .^ 2 rtol = 1.0e-6
        end
    end
end

@testset "Initialize scalar equations with array unknowns" begin
    @parameters t
    @variables u(t)[1:6]
    D = Differential(t)
    dx = 0.2
    eqs = [D(u[i]) ~ (u[i - 1] - 2u[i] + u[i + 1]) / dx^2 for i in 2:5]
    append!(eqs, [u[1] ~ 0.0, u[6] ~ 0.0])
    @named raw = System(eqs, t, [u], [])
    sys = complete(raw)
    values = vcat(0.0, sinpi.(collect(1:4) .* dx), 0.0)
    eigenvalue = -4sinpi(dx / 2)^2 / dx^2
    prob = DAEProblem(sys, [u => values, D(u[1]) => 0.0, D(u[6]) => 0.0], (0.0, 0.1))
    @test length(unknowns(sys)) == 1
    @test prob.u0 == values
    @test prob.du0 ≈ eigenvalue .* values
    @test prob.differential_vars == [false, true, true, true, true, false]
    sol = solve(prob, DFBDF(); abstol = 1.0e-9, reltol = 1.0e-9)
    @test SciMLBase.successful_retcode(sol)
    @test sol[u][end] ≈ exp(0.1eigenvalue) .* values rtol = 1.0e-6
end
