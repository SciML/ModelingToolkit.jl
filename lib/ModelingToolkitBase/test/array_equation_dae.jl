using ModelingToolkitBase, Test
using ModelingToolkitBase: unwrap, complete, unknowns
using Symbolics
using SciMLBase
using OrdinaryDiffEqBDF: DFBDF
using DiffEqBase: BrownFullBasicInit

# A system whose interior is written as one array equation over slices, as produced by a
# finite-difference PDE discretization that does not scalarize.
function heat_array_system(n)
    @independent_variables t
    @variables u(t)[1:n]
    D = Differential(t)
    dx = 1 / (n - 1)
    # Residual (cardinalized) form, as a finite-difference discretization emits it: the
    # derivative sits inside the expression rather than being the equation's whole LHS.
    lap = (u[1:(n - 2)] .- 2 .* u[2:(n - 1)] .+ u[3:n]) ./ dx^2
    interior = broadcast(-, D(u[2:(n - 1)]), lap) ~ zeros(n - 2)
    eqs = [interior, u[1] ~ 0.0, u[n] ~ 0.0]
    @named sys = System(eqs, t, collect(u), [])
    return complete(sys), u, t, D
end

@testset "array equations reach DAEProblem" begin
    n = 11
    sys, u, t, D = heat_array_system(n)
    xs = range(0.0, 1.0, length = n)
    op = vcat(
        [u[i] => sinpi(xs[i]) for i in 1:n],
        [D(u[i]) => 0.0 for i in 1:n]
    )

    prob = DAEProblem(sys, op, (0.0, 0.1); build_initializeprob = false)

    # one output row per element of the array equation, not one per equation
    @test length(prob.u0) == n
    @test prob.u0 isa Vector{Float64}

    # the interior points are differential, the two boundary points algebraic
    @test prob.differential_vars !== nothing
    @test count(prob.differential_vars) == n - 2

    # the residual evaluates: no `Differential` survives into the generated code
    out = zeros(n)
    du = zeros(n)
    prob.f(out, du, prob.u0, prob.p, 0.0)
    @test all(isfinite, out)
    # with du = 0 the interior residual is minus the Laplacian, which is nonzero here
    @test any(!iszero, out)
end

@testset "other problem types still require scalarized equations" begin
    n = 11
    sys, u, t, D = heat_array_system(n)
    op = [u[i] => 0.0 for i in 1:n]
    # ODEProblem cannot consume array equations; the guard must remain in place
    @test_throws Exception ODEProblem(sys, op, (0.0, 0.1); build_initializeprob = false)
end

@testset "array-equation DAE solves to the analytic solution" begin
    n = 21
    sys, u, t, D = heat_array_system(n)
    xs = range(0.0, 1.0, length = n)
    op = vcat(
        [u[i] => sinpi(xs[i]) for i in 1:n],
        [D(u[i]) => 0.0 for i in 1:n]
    )
    tend = 0.1
    prob = DAEProblem(sys, op, (0.0, tend); build_initializeprob = false)
    # `du0` above is not consistent; the solver's own DAE initialization supplies it.
    sol = solve(
        prob, DFBDF(); initializealg = BrownFullBasicInit(),
        reltol = 1.0e-8, abstol = 1.0e-8, saveat = [tend]
    )
    @test SciMLBase.successful_retcode(sol)
    exact = [exp(-pi^2 * tend) * sinpi(xi) for xi in xs]
    # second-order spatial discretization on 21 points
    @test maximum(abs.(sol.u[end] .- exact)) < 5.0e-3
end

@testset "array equations over a 2D slice keep their shape" begin
    # A derivative of a 2D slice must substitute a 2D array of scalar derivatives; a
    # flattened one does not broadcast against the surrounding slices and codegen fails
    # with a DimensionMismatch.
    n = 6
    @independent_variables t
    @variables w(t)[1:n, 1:n]
    D = Differential(t)
    dx = 1 / (n - 1)
    inner = 2:(n - 1)
    lap = (
        w[1:(n - 2), inner] .+ w[3:n, inner] .+ w[inner, 1:(n - 2)] .+
            w[inner, 3:n] .- 4 .* w[inner, inner]
    ) ./ dx^2
    eqs = Equation[broadcast(-, D(w[inner, inner]), lap) ~ zeros(n - 2, n - 2)]
    for i in 1:n
        push!(eqs, w[i, 1] ~ 0.0)
        push!(eqs, w[i, n] ~ 0.0)
    end
    for j in inner
        push!(eqs, w[1, j] ~ 0.0)
        push!(eqs, w[n, j] ~ 0.0)
    end
    @named sys2d = System(eqs, t, vec(collect(w)), [])
    sys2d = complete(sys2d)

    op = vcat(
        [w[i, j] => 0.25 for i in 1:n, j in 1:n] |> vec,
        [D(w[i, j]) => 0.0 for i in 1:n, j in 1:n] |> vec
    )
    prob = DAEProblem(sys2d, op, (0.0, 0.01); build_initializeprob = false)
    @test length(prob.u0) == n * n
    out = zeros(n * n)
    prob.f(out, zeros(n * n), prob.u0, prob.p, 0.0)
    @test all(isfinite, out)
end

@testset "array equations written as `D(slice) ~ rhs`" begin
    # The residual form above puts the derivative inside the expression. The equivalent
    # `D(u[2:n-1]) ~ rhs` form has no scalar `toterm` name for its LHS, which the
    # derivative-substitution machinery has to skip rather than trip over.
    n = 11
    @independent_variables t
    @variables u(t)[1:n]
    D = Differential(t)
    dx = 1 / (n - 1)
    lap = (u[1:(n - 2)] .- 2 .* u[2:(n - 1)] .+ u[3:n]) ./ dx^2
    eqs = [D(u[2:(n - 1)]) ~ lap, u[1] ~ 0.0, u[n] ~ 0.0]
    @named sys = System(eqs, t, collect(u), [])
    sys = complete(sys)

    xs = range(0.0, 1.0, length = n)
    op = vcat([u[i] => sinpi(xs[i]) for i in 1:n], [D(u[i]) => 0.0 for i in 1:n])
    prob = DAEProblem(sys, op, (0.0, 0.1); build_initializeprob = false)
    @test length(prob.u0) == n

    # the residual matches the analytic derivative of the initial condition
    out = zeros(n)
    du = zeros(n)
    du[2:(n - 1)] .= [-pi^2 * sinpi(x) for x in xs[2:(n - 1)]]
    prob.f(out, du, prob.u0, prob.p, 0.0)
    @test maximum(abs, out) < 1.0e-1

    sol = solve(
        prob, DFBDF(); initializealg = BrownFullBasicInit(), reltol = 1.0e-8,
        abstol = 1.0e-8, saveat = [0.1]
    )
    @test SciMLBase.successful_retcode(sol)
    @test maximum(abs, sol.u[end] .- [exp(-pi^2 * 0.1) * sinpi(x) for x in xs]) < 1.0e-2
end
