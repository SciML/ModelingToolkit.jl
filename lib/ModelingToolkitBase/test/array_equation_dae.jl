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

@testset "array-keyed DAE operating points preserve element order" begin
    n = 4
    @independent_variables t
    @variables u(t)[1:n, 1:n]
    D = Differential(t)
    @named sys = System([D(u) ~ zeros(n, n)], t, [u], [])
    sys = complete(sys)
    uval = reshape(Float64.(1:(n * n)), n, n)
    duval = reshape(Float64.((n * n + 1):(2 * n * n)), n, n)

    prob = DAEProblem(
        sys, [u => uval, D(u) => duval], (0.0, 1.0); build_initializeprob = false
    )

    @test prob.u0 == vec(uval)
    @test prob.du0 == vec(duval)
end

@testset "array operating points substitute specified symbolic elements" begin
    @independent_variables t
    @variables x(t)[1:3] y(t)
    @parameters p q
    D = Differential(t)
    @named sys = System([D(x) ~ -x, D(y) ~ -y], t, [x, y], [p, q])
    sys = complete(sys)
    ir = ModelingToolkitBase.get_irstructure(sys)
    dvs = Symbolics.unwrap.([x[1], x[2], x[3], y])

    constant_filled = ModelingToolkitBase.varmap_to_vars(
        Dict(x[1] => 2p, p => 3.0, y => 1.0), dvs; ir,
        missing_values = ModelingToolkitBase.MissingGuessValue.Constant(0.0)
    )
    @test constant_filled == [6.0, 0.0, 0.0, 1.0]

    hashed_filled = ModelingToolkitBase.varmap_to_vars(
        Dict(x[1] => 2p, x[3] => y, p => 3.0, y => 1.0), dvs; ir,
        missing_values = ModelingToolkitBase.MissingGuessValue.HashedRandom()
    )
    @test hashed_filled[[1, 3, 4]] == [6.0, 1.0, 1.0]

    fully_specified = ModelingToolkitBase.varmap_to_vars(
        Dict(x[1] => 2p, x[2] => y + 1, x[3] => x[1] + q, p => 3.0, q => 10.0, y => 1.0),
        dvs; ir
    )
    @test fully_specified == [6.0, 2.0, 16.0, 1.0]

    # `vars` may mix a whole-array key with element keys of the same array
    for (whole_idx, mixed_vars) in ((1, [x, x[1]]), (2, [x[1], x]))
        mixed = ModelingToolkitBase.varmap_to_vars(
            Dict(x => [p, 2p, 3p], p => 3.0), Symbolics.unwrap.(mixed_vars);
            ir, allow_symbolic = true
        )
        @test mixed[whole_idx] == [3.0, 6.0, 9.0]
        @test mixed[3 - whole_idx] == 3.0
    end

    # the non-IR path substitutes specified elements of partially filled arrays too
    non_ir = ModelingToolkitBase.varmap_to_vars(
        Dict(x[1] => 2p, p => 3.0, y => 1.0), dvs;
        missing_values = ModelingToolkitBase.MissingGuessValue.Constant(0.0)
    )
    @test non_ir == [6.0, 0.0, 0.0, 1.0]

    @variables u(t)[1:2, 1:2]
    @parameters r
    @named sys2 = System([D(u) ~ -u], t, [u], [r])
    sys2 = complete(sys2)
    ir2 = ModelingToolkitBase.get_irstructure(sys2)
    dvs2 = Symbolics.unwrap.([u[1, 1], u[2, 1], u[1, 2], u[2, 2]])
    twod = ModelingToolkitBase.varmap_to_vars(
        Dict(u[1, 1] => r, u[2, 2] => 4r, r => 1.5), dvs2; ir = ir2,
        missing_values = ModelingToolkitBase.MissingGuessValue.Constant(0.0)
    )
    @test twod == [1.5, 0.0, 0.0, 6.0]
end

@testset "acyclic cross-element array operating points resolve incrementally" begin
    # Chains longer than `substitution_limit` must still resolve when elements are
    # visited in index order with immediate write-back (master semantics). A single
    # whole-array fixpoint wrongly leaves the tail missing / filled with the default.
    @independent_variables t
    @variables x(t)[1:110]
    D = Differential(t)
    @named chain_sys = System([D(x) ~ -x], t, [x], [])
    chain_sys = complete(chain_sys)
    op = Dict(x[i] => x[i - 1] + 1 for i in 2:110)
    op[x[1]] = 1.0
    ir = ModelingToolkitBase.get_irstructure(chain_sys)
    expected = collect(1.0:110.0)
    for options in ((; ir), (;))
        values = ModelingToolkitBase.varmap_to_vars(
            op, Symbolics.unwrap.(collect(x)); options...,
            missing_values = ModelingToolkitBase.MissingGuessValue.Constant(-999.0)
        )
        @test values == expected
    end
    prob = ODEProblem(chain_sys, op, (0.0, 1.0); build_initializeprob = false)
    @test prob.u0 == expected
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
