using ModelingToolkitBase, Test
using ModelingToolkitBase: unwrap, complete, unknowns
using ModelingToolkitBase: has_array_equations, accepts_array_equations
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

@testset "DAE initialization accepts derivative guesses" begin
    @independent_variables t
    @variables y(t)[1:3] z(t)
    D = Differential(t)
    @named sys = System(
        [zeros(3) ~ D(y[1:3]) .+ z .* y[1:3], 0 ~ z - sum(y)],
        t, [collect(y); z], []
    )
    sys = complete(sys)
    y0 = [1.0, 2.0, 3.0]
    op = [y[i] => y0[i] for i in 1:3]

    for provide_derivative_guesses in (true, false)
        guesses = [z => 0.0]
        if provide_derivative_guesses
            append!(guesses, [D(y[i]) => 0.5 for i in 1:3])
        end
        prob = DAEProblem(sys, op, (0.0, 0.1); guesses)
        @test prob.u0[4] == 0.0
        @test prob.du0[1:3] == fill(provide_derivative_guesses ? 0.5 : 0.0, 3)
        initprob = prob.f.initialization_data.initializeprob
        init_sol = solve(initprob)

        @test SciMLBase.successful_retcode(init_sol)
        init_unknowns = unknowns(initprob.f.sys)
        z_index = findfirst(isequal(z), init_unknowns)
        derivative_values = init_sol.u[setdiff(eachindex(init_unknowns), [z_index])]
        @test sort(derivative_values) ≈ sort(-sum(y0) .* y0)
        @test init_sol.u[z_index] ≈ sum(y0)

        integ = init(prob, DFBDF(); initializealg = BrownFullBasicInit())
        @test integ.du[1:3] ≈ -sum(y0) .* y0
        residual = zeros(4)
        prob.f(residual, integ.du, integ.u, integ.p, 0.0)
        @test maximum(abs, residual) < 1.0e-10
    end

    fixed_du0 = [D(y[i]) => 0.25 for i in 1:3]
    prob = DAEProblem(
        sys, [op; [z => 0.0]; fixed_du0], (0.0, 0.1); build_initializeprob = false
    )
    @test prob.du0[1:3] == fill(0.25, 3)
end

@testset "array unknowns flatten under an array equation" begin
    n = 11
    @independent_variables t
    @variables u(t)[1:n]
    D = Differential(t)
    @named sys = System([D(u) ~ -u], t, [u], [])
    sys = complete(sys)
    op = [u => zeros(n), D(u) => zeros(n)]
    # `flat_unknowns` flattens `u`, so the array equation over the array unknown
    # constructs directly
    prob = ODEProblem(sys, op, (0.0, 0.1); build_initializeprob = false)
    @test length(prob.u0) == n
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

@testset "has_array_equations detects every array-equation form" begin
    n = 5
    @independent_variables t
    @variables u(t)[1:n]
    D = Differential(t)
    lap = u[1:(n - 2)] .- 2 .* u[2:(n - 1)] .+ u[3:n]
    @test has_array_equations([zeros(n - 2) ~ broadcast(-, lap)])
    @test has_array_equations([broadcast(-, lap) ~ zeros(n - 2)])
    @test has_array_equations([D(u[2:(n - 1)]) ~ lap])
    @test !has_array_equations([u[1] ~ 0.0, u[n] ~ 0.0])
    @test !has_array_equations(Equation[])
end

@testset "accepting array equations is a per-constructor capability" begin
    @test accepts_array_equations(DAEFunction)
    @test accepts_array_equations(NonlinearFunction)
    @test accepts_array_equations(ODEFunction)
    @test !accepts_array_equations(SDEFunction)
    @test !accepts_array_equations(ImplicitDiscreteFunction)
    # Optimization vectorizes `costs` and `constraints`, not `equations`:
    # `check_no_equations` rejects them before this gate is reached.
    @test !accepts_array_equations(OptimizationFunction)
    @test !accepts_array_equations(MultiObjectiveOptimizationFunction)
end

@testset "symbolic jacobian from array residuals" begin
    n = 11
    sys, u, t, D = heat_array_system(n)
    xs = range(0.0, 1.0, length = n)
    op = vcat(
        [u[i] => sinpi(xs[i]) for i in 1:n],
        [D(u[i]) => 0.0 for i in 1:n]
    )
    # the jacobian is built from the scalarized `full_equations`, one row per residual;
    # `sparse = true` is not tested because `W_sparsity` requires a semi-explicit mass
    # matrix, which no residual-form DAE has
    prob = DAEProblem(sys, op, (0.0, 0.1); build_initializeprob = false, jac = true)
    J = zeros(n, n)
    γ = 2.0
    prob.f.jac(J, prob.du0, prob.u0, prob.p, γ, 0.0)
    dx = 1 / (n - 1)
    # residual row 1 is `lap[1] - D(u[2])`
    @test J[1, 1:3] ≈ [1, -2 - γ * dx^2, 1] ./ dx^2
    # residual row `n - 1` is `0 - u[1]`
    @test J[n - 1, 1] ≈ -1
    @test count(!iszero, J) == 3 * (n - 2) + 2
    sol = solve(
        prob, DFBDF(); initializealg = BrownFullBasicInit(), reltol = 1.0e-8,
        abstol = 1.0e-8, saveat = [0.1]
    )
    @test SciMLBase.successful_retcode(sol)
    @test maximum(abs, sol.u[end] .- [exp(-pi^2 * 0.1) * sinpi(x) for x in xs]) < 1.0e-2
end
