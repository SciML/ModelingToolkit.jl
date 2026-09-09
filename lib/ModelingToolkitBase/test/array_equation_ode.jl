using ModelingToolkitBase, Test
using ModelingToolkitBase: complete, unknowns, generate_rhs, GeneratedFunctionOptions,
    calculate_massmatrix, full_equations, scalarize_array_equations
using Symbolics
using SciMLBase
using LinearAlgebra
using SparseArrays: nnz
using OrdinaryDiffEqTsit5: Tsit5
using OrdinaryDiffEqRosenbrock: Rodas5P

# Array differential equations reach `ODEProblem` from `complete(sys)` directly, without
# `mtkcompile` scalarizing them: the equation stays one equation, and its right-hand side
# is written to a contiguous block of `du`.

# Interior of the 1D heat equation as one array equation over slices, as a
# finite-difference discretization emits it. `:explicit` is `D(u[2:n-1]) ~ lap`,
# `:residual` is the cardinalized `D(u[2:n-1]) .- lap ~ 0`.
function heat_array_system(n, form)
    @independent_variables t
    @variables u(t)[1:n]
    D = Differential(t)
    dx = 1 / (n - 1)
    lap = (u[1:(n - 2)] .- 2 .* u[2:(n - 1)] .+ u[3:n]) ./ dx^2
    interior = if form === :explicit
        D(u[2:(n - 1)]) ~ lap
    else
        broadcast(-, D(u[2:(n - 1)]), lap) ~ zeros(n - 2)
    end
    eqs = [0 ~ u[1], interior, 0 ~ u[n]]
    @named sys = System(eqs, t, collect(u), [])
    return sys, u, t, D
end

count_expr_nodes(ex) = ex isa Expr ? 1 + sum(count_expr_nodes, ex.args; init = 0) : 0
function contains_symbol(ex, s::Symbol)
    ex === s && return true
    ex isa Expr || return false
    return any(arg -> contains_symbol(arg, s), ex.args)
end

@testset "`D(u) ~ f` over a whole array" begin
    n = 5
    @independent_variables t
    @variables u(t)[1:n]
    D = Differential(t)
    @named sys = System([D(u) ~ -u], t, collect(u), [])
    sys = complete(sys)
    @test length(equations(sys)) == 1
    @test length(full_equations(sys)) == n

    prob = ODEProblem(sys, [u => ones(n)], (0.0, 1.0); build_initializeprob = false)
    @test prob.f.mass_matrix === I
    du = zeros(n)
    prob.f(du, prob.u0, prob.p, 0.0)
    @test du ≈ -ones(n)
    @test prob.f(prob.u0, prob.p, 0.0) ≈ -ones(n)

    sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-8)
    @test SciMLBase.successful_retcode(sol)
    @test sol.u[end] ≈ fill(exp(-1.0), n) rtol = 1.0e-6

    # the default initialization problem is built from the scalarized equations
    prob2 = ODEProblem(sys, [u => ones(n)], (0.0, 1.0))
    @test prob2.f.initialization_data !== nothing
    sol2 = solve(prob2, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-8)
    @test SciMLBase.successful_retcode(sol2)
    @test sol2.u[end] ≈ fill(exp(-1.0), n) rtol = 1.0e-6
end

@testset "heat equation over a slice, $form form" for form in (:explicit, :residual)
    n = 21
    sys, u, t, D = heat_array_system(n, form)
    sys = complete(sys)
    xs = range(0.0, 1.0, length = n)
    tend = 0.1
    u0 = sinpi.(xs)

    M = calculate_massmatrix(sys)
    @test M isa Diagonal
    @test diag(M) == [0; ones(n - 2); 0]

    prob = ODEProblem(sys, [u => u0], (0.0, tend); build_initializeprob = false)
    @test length(prob.u0) == n
    @test prob.f.mass_matrix isa Diagonal
    du = zeros(n)
    prob.f(du, prob.u0, prob.p, 0.0)
    dx = step(xs)
    lap = (u0[1:(n - 2)] .- 2 .* u0[2:(n - 1)] .+ u0[3:n]) ./ dx^2
    @test du[2:(n - 1)] ≈ lap
    @test du[1] == u0[1] == 0.0
    @test du[n] == u0[n] == 0.0
    @test prob.f(prob.u0, prob.p, 0.0) ≈ du

    sol = solve(prob, Rodas5P(); reltol = 1.0e-8, abstol = 1.0e-8, saveat = [tend])
    @test SciMLBase.successful_retcode(sol)
    exact = [exp(-pi^2 * tend) * sinpi(xi) for xi in xs]
    # second-order spatial discretization on 21 points
    @test maximum(abs.(sol.u[end] .- exact)) < 1.0e-3

    # parity with the scalarized equations, through `mtkcompile` (which, without
    # ModelingToolkit loaded, needs the explicit form) and through `complete`
    explicit_sys, = heat_array_system(n, :explicit)
    csys = mtkcompile(explicit_sys)
    @test all(eq -> !Symbolics.isarraysymbolic(eq.lhs), equations(csys))
    cprob = ODEProblem(csys, [u => u0], (0.0, tend); build_initializeprob = false)
    csol = solve(cprob, Rodas5P(); reltol = 1.0e-8, abstol = 1.0e-8, saveat = [tend])
    @test sol[u][end] ≈ csol[u][end] atol = 1.0e-10

    scalar_sys, = heat_array_system(n, form)
    @named ssys = System(scalarize_array_equations(equations(scalar_sys)), t, collect(u), [])
    ssys = complete(ssys)
    @test length(equations(ssys)) == n
    sprob = ODEProblem(ssys, [u => u0], (0.0, tend); build_initializeprob = false)
    @test sprob.f.mass_matrix == prob.f.mass_matrix
    sdu = zeros(n)
    sprob.f(sdu, sprob.u0, sprob.p, 0.0)
    @test sdu ≈ du

    # the default initialization problem handles the array equation
    iprob = ODEProblem(sys, [u => u0], (0.0, tend); warn_initialize_determined = false)
    isol = solve(iprob, Rodas5P(); reltol = 1.0e-8, abstol = 1.0e-8, saveat = [tend])
    @test SciMLBase.successful_retcode(isol)
    @test isol[u][end] ≈ csol[u][end] atol = 1.0e-10
end

@testset "generated code is independent of the array length" begin
    sizes = map((24, 48, 96)) do n
        sys, = heat_array_system(n, :explicit)
        sys = complete(sys)
        f_oop, f_iip = generate_rhs(sys, GeneratedFunctionOptions(; expression = Val{true}))
        @test !contains_symbol(f_iip, :array_literal)
        (count_expr_nodes(f_oop), count_expr_nodes(f_iip))
    end
    @test allequal(sizes)
end

@testset "array equations over a 2D slice" begin
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
    eqs = Equation[D(w[inner, inner]) ~ lap]
    for i in 1:n
        push!(eqs, 0 ~ w[i, 1])
        push!(eqs, 0 ~ w[i, n])
    end
    for j in inner
        push!(eqs, 0 ~ w[1, j])
        push!(eqs, 0 ~ w[n, j])
    end
    @named sys2d = System(eqs, t, vec(collect(w)), [])
    sys2d = complete(sys2d)
    w0 = [sinpi(x) * sinpi(y) for x in range(0, 1, length = n), y in range(0, 1, length = n)]

    prob = ODEProblem(sys2d, [w => w0], (0.0, 0.01); build_initializeprob = false)
    @test length(prob.u0) == n * n
    @test count(isone, prob.f.mass_matrix) == (n - 2)^2
    du = zeros(n * n)
    prob.f(du, prob.u0, prob.p, 0.0)
    @test all(isfinite, du)

    # same rows as the scalarized equations
    @named ssys = System(scalarize_array_equations(eqs), t, vec(collect(w)), [])
    ssys = complete(ssys)
    sprob = ODEProblem(ssys, [w => w0], (0.0, 0.01); build_initializeprob = false)
    @test sprob.f.mass_matrix == prob.f.mass_matrix
    sdu = zeros(n * n)
    sprob.f(sdu, sprob.u0, sprob.p, 0.0)
    @test sdu ≈ du

    sol = solve(prob, Rodas5P())
    @test SciMLBase.successful_retcode(sol)
end

@testset "symbolic jacobian and sparsity from array equations" begin
    n = 11
    sys, u, t, D = heat_array_system(n, :explicit)
    sys = complete(sys)
    u0 = sinpi.(range(0.0, 1.0, length = n))
    prob = ODEProblem(
        sys, [u => u0], (0.0, 0.1);
        jac = true, sparse = true, build_initializeprob = false
    )
    J = similar(prob.f.jac_prototype)
    prob.f.jac(J, prob.u0, prob.p, 0.0)
    # tridiagonal interior rows plus the two boundary rows
    @test size(J) == (n, n)
    @test nnz(J) == 3 * (n - 2) + 2
    dx = 1 / (n - 1)
    @test J[2, 1] ≈ 1 / dx^2
    @test J[2, 2] ≈ -2 / dx^2
    @test J[1, 1] ≈ 1
    sol = solve(prob, Rodas5P(); reltol = 1.0e-8, abstol = 1.0e-8, saveat = [0.1])
    @test SciMLBase.successful_retcode(sol)
end

@testset "`full_equations` expands array equations into rows" begin
    n = 7
    for form in (:explicit, :residual)
        sys, u, t, D = heat_array_system(n, form)
        sys = complete(sys)
        eqs = full_equations(sys)
        @test length(eqs) == n
        @test all(eq -> !Symbolics.isarraysymbolic(eq.lhs), eqs)
        # the residual form is rewritten to `D(u[i]) ~ f_i`
        @test all(isequal(D(u[i]), eqs[i].lhs) for i in 2:(n - 1))
    end
end

@testset "invalid array differential equations are rejected" begin
    @independent_variables t
    @variables u(t)[1:4]
    D = Differential(t)
    opts = GeneratedFunctionOptions(; expression = Val{true})

    # a derivative that appears both in a slice and on its own
    @named dup = System([D(u[1:3]) ~ -u[1:3], D(u[3]) ~ 0, 0 ~ u[4]], t, collect(u), [])
    @test_throws ["LHS operator must be unique"] generate_rhs(complete(dup), opts)

    # a slice whose elements are not all unknowns
    @named notunknown = System([D(u[1:3]) ~ -u[1:3]], t, collect(u[1:2]), [])
    @test_throws ["not a valid LHS"] generate_rhs(complete(notunknown), opts)

    # array unknowns still need `mtkcompile`
    @named arrunknown = System([D(u) ~ -u], t, [u], [])
    @test_throws ["array unknowns"] ODEProblem(
        complete(arrunknown), [u => ones(4)], (0.0, 1.0); build_initializeprob = false
    )
end
