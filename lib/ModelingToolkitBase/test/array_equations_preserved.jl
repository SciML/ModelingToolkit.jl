using ModelingToolkitBase, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using ModelingToolkitBase: unwrap, arrays_scalarized, default_toterm, AtomicArrayDict,
    write_possibly_indexed_array!, generate_rhs, GeneratedFunctionOptions,
    explicit_array_differential_equation
using Symbolics, SymbolicUtils, LinearAlgebra
using SciMLBase
using OrdinaryDiffEqBDF: DFBDF
using OrdinaryDiffEqRosenbrock: Rodas5P

# Method-of-lines discretization of the heat equation with Dirichlet boundaries, written
# as one array equation over slices of the unknowns plus scalar boundary equations.
function heat_slice_system(n)
    @variables u(t)[1:n]
    @parameters α = 1.0
    dx = 1 / (n - 1)
    lap = (u[1:(n - 2)] .- 2 .* u[2:(n - 1)] .+ u[3:n]) ./ dx^2
    eqs = [D(u[2:(n - 1)]) ~ α .* lap, u[1] ~ 0.0, u[n] ~ 0.0]
    sys = System(eqs, t, [u], [α]; name = :heat)
    return sys, u, α
end

# Number of nodes in the generated code, as a proxy for its size.
count_expr_nodes(ex) = ex isa Expr ? 1 + sum(count_expr_nodes, ex.args; init = 0) : 1

@testset "`scalarize_arrays = false` keeps array equations" begin
    for n in (6, 24, 96)
        sys, u, α = heat_slice_system(n)
        ssys = mtkcompile(sys; scalarize_arrays = false)
        @test !arrays_scalarized(ssys)
        @test arrays_scalarized(sys)
        @test length(equations(ssys)) == 1
        eq = only(equations(ssys))
        @test isequal(eq.lhs, unwrap(D(u[2:(n - 1)])))
        @test length(unknowns(ssys)) == n - 2
        @test isequal(unknowns(ssys), [unwrap(u[i]) for i in 2:(n - 1)])
        @test length(observed(ssys)) == 2
        @test issetequal([eq.lhs for eq in observed(ssys)], [unwrap(u[1]), unwrap(u[n])])
    end

    @testset "default path still scalarizes" begin
        sys, u, α = heat_slice_system(6)
        ssys = mtkcompile(sys)
        @test arrays_scalarized(ssys)
        @test length(equations(ssys)) == 4
    end
end

@testset "`toterm` of a derivative of a slice" begin
    @variables u(t)[1:6]
    ttk = default_toterm(unwrap(D(u[2:4])))
    @test SymbolicUtils.is_array_shape(SymbolicUtils.shape(ttk))
    @test isequal(ttk[1], default_toterm(unwrap(D(u[2]))))
    @test isequal(ttk[3], default_toterm(unwrap(D(u[4]))))
end

@testset "scalar value written to an array key is broadcast" begin
    @variables x(t)[1:3]
    dd = AtomicArrayDict{SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}}(
        Dict{SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}, SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal}}()
    )
    zero = SymbolicUtils.Const{SymbolicUtils.SymReal}(false)
    nothing_sym = SymbolicUtils.Const{SymbolicUtils.SymReal}(nothing)
    write_possibly_indexed_array!(dd, unwrap(x), zero, nothing_sym)
    @test SymbolicUtils.unwrap_const(dd[unwrap(x)]) == [false, false, false]
    write_possibly_indexed_array!(dd, unwrap(x[2:3]), SymbolicUtils.Const{SymbolicUtils.SymReal}(true), nothing_sym)
    @test SymbolicUtils.unwrap_const(dd[unwrap(x)]) == [false, true, true]
end

@testset "generated DAE code does not grow with the array length" begin
    sizes = Dict{Int, Int}()
    for n in (24, 48, 96)
        sys, u, α = heat_slice_system(n)
        ssys = mtkcompile(sys; scalarize_arrays = false)
        oop, iip = generate_rhs(
            ssys, GeneratedFunctionOptions(; expression = Val{true}); implicit_dae = true
        )
        sizes[n] = count_expr_nodes(iip)
        # no per-element reads of `du` or the unknowns
        @test !occursin("array_literal", string(iip))
    end
    @test sizes[24] == sizes[48] == sizes[96]
end

@testset "DAEProblem after `mtkcompile(sys; scalarize_arrays = false)`" begin
    n = 24
    sys, u, α = heat_slice_system(n)
    ssys = mtkcompile(sys; scalarize_arrays = false)
    xs = range(0.0, 1.0, length = n)
    u0 = sinpi.(xs)
    tend = 0.1
    exact = [exp(-pi^2 * tend) * sinpi(x) for x in xs]

    # consistent initial derivatives: the discrete Laplacian of `u0` in the interior
    dx = 1 / (n - 1)
    du0 = zeros(n)
    du0[2:(n - 1)] .= (u0[1:(n - 2)] .- 2 .* u0[2:(n - 1)] .+ u0[3:n]) ./ dx^2
    prob = DAEProblem(ssys, [u => u0, D(u) => du0], (0.0, tend); build_initializeprob = false)
    @test length(prob.u0) == n - 2
    @test prob.u0 == u0[2:(n - 1)]
    @test prob.du0 == du0[2:(n - 1)]
    resid = zeros(n - 2)
    du = zeros(n - 2)
    prob.f(resid, du, prob.u0, prob.p, 0.0)
    # with `du = 0` the residual is the discrete Laplacian of `sinpi`, `-π² sinpi`
    @test maximum(abs, resid .- (-pi^2 .* u0[2:(n - 1)])) < 0.05
    # the out-of-place form agrees
    @test prob.f(du, prob.u0, prob.p, 0.0) ≈ resid
    # and the consistent `du0` gives a zero residual
    prob.f(resid, prob.du0, prob.u0, prob.p, 0.0)
    @test maximum(abs, resid) < 1.0e-10

    sol = solve(
        prob, DFBDF(); initializealg = SciMLBase.NoInit(), reltol = 1.0e-8,
        abstol = 1.0e-8, saveat = [tend]
    )
    @test SciMLBase.successful_retcode(sol)
    # observed boundary values and the whole array are accessible
    @test sol[u[1]][end] == 0.0
    @test sol[u[n]][end] == 0.0
    @test maximum(abs, sol[u][end] .- exact) < 5.0e-3

    @testset "matches the scalarized compilation" begin
        ssys2 = mtkcompile(sys)
        prob2 = ODEProblem(ssys2, [u => u0], (0.0, tend))
        sol2 = solve(prob2, Rodas5P(); reltol = 1.0e-8, abstol = 1.0e-8, saveat = [tend])
        @test maximum(abs, sol[u][end] .- sol2[u][end]) < 1.0e-5
    end

    @testset "with initialization" begin
        prob = DAEProblem(ssys, [u => u0], (0.0, tend); warn_initialize_determined = false)
        @test prob.f.initialization_data !== nothing
        # the initialization determines the derivatives of the interior
        @test maximum(abs, prob.du0 .- (-pi^2 .* u0[2:(n - 1)])) < 0.05
        sol = solve(prob, DFBDF(); reltol = 1.0e-8, abstol = 1.0e-8, saveat = [tend])
        @test SciMLBase.successful_retcode(sol)
        @test maximum(abs, sol[u][end] .- exact) < 5.0e-3
    end
end

@testset "2D slice" begin
    n = 8
    @variables w(t)[1:n, 1:n]
    dx = 1 / (n - 1)
    inner = 2:(n - 1)
    lap = (
        w[1:(n - 2), inner] .+ w[3:n, inner] .+ w[inner, 1:(n - 2)] .+
            w[inner, 3:n] .- 4 .* w[inner, inner]
    ) ./ dx^2
    eqs = Equation[D(w[inner, inner]) ~ lap]
    for i in 1:n
        push!(eqs, w[i, 1] ~ 0.0)
        push!(eqs, w[i, n] ~ 0.0)
    end
    for j in inner
        push!(eqs, w[1, j] ~ 0.0)
        push!(eqs, w[n, j] ~ 0.0)
    end
    @named sys2d = System(eqs, t, [w], [])
    ssys = mtkcompile(sys2d; scalarize_arrays = false)
    @test length(equations(ssys)) == 1
    @test length(unknowns(ssys)) == (n - 2)^2
    @test length(observed(ssys)) == 4n - 4

    xs = range(0.0, 1.0, length = n)
    w0 = [sinpi(x) * sinpi(y) for x in xs, y in xs]
    prob = DAEProblem(ssys, [w => w0, D(w) => zeros(n, n)], (0.0, 0.01); build_initializeprob = false)
    resid = zeros((n - 2)^2)
    prob.f(resid, zeros((n - 2)^2), prob.u0, prob.p, 0.0)
    # `-2π² sinpi(x) sinpi(y)` up to discretization error
    @test maximum(abs, resid .- vec(-2pi^2 .* w0[inner, inner])) < 1.0
end

@testset "residual form `D(x) .- f ~ 0` is accepted" begin
    n = 24
    @variables u(t)[1:n]
    dx = 1 / (n - 1)
    lap = (u[1:(n - 2)] .- 2 .* u[2:(n - 1)] .+ u[3:n]) ./ dx^2
    der = unwrap(D(u[2:(n - 1)]))
    @test isequal(
        explicit_array_differential_equation(D(u[2:(n - 1)]) .- lap ~ zeros(n - 2)),
        der ~ unwrap(lap)
    )
    @test isequal(
        explicit_array_differential_equation(zeros(n - 2) ~ D(u[2:(n - 1)]) .- lap),
        der ~ unwrap(lap)
    )
    @test isequal(
        explicit_array_differential_equation(lap .- D(u[2:(n - 1)]) ~ zeros(n - 2)),
        der ~ unwrap(lap)
    )
    @test isequal(
        explicit_array_differential_equation(D(u[2:(n - 1)]) .+ lap ~ zeros(n - 2)),
        der ~ unwrap(-1 .* lap)
    )
    @test isequal(
        explicit_array_differential_equation(D(u[2:(n - 1)]) .- lap ~ u[2:(n - 1)]),
        der ~ unwrap(u[2:(n - 1)] .+ lap)
    )
    # the derivative must not appear elsewhere; other equations are left alone
    eq = D(u[2:(n - 1)]) .- lap ~ D(u[2:(n - 1)])
    @test explicit_array_differential_equation(eq) === eq
    eq = u[2:(n - 1)] ~ lap
    @test explicit_array_differential_equation(eq) === eq

    # as MethodOfLines emits it: a residual array equation and boundary residuals
    eqs = [D(u[2:(n - 1)]) .- lap ~ zeros(n - 2), u[1] - 0.0 ~ 0.0, u[n] - 0.0 ~ 0.0]
    @named heat = System(eqs, t, [u], [])
    ssys = mtkcompile(heat; scalarize_arrays = false)
    @test length(equations(ssys)) == 1
    @test isequal(only(equations(ssys)).lhs, der)
    @test length(unknowns(ssys)) == n - 2
    @test length(observed(ssys)) == 2

    xs = range(0.0, 1.0, length = n)
    u0 = sinpi.(xs)
    prob = DAEProblem(ssys, [u => u0], (0.0, 0.1); warn_initialize_determined = false)
    sol = solve(prob, DFBDF(); reltol = 1.0e-8, abstol = 1.0e-8, saveat = [0.1])
    @test SciMLBase.successful_retcode(sol)
    @test maximum(abs, sol[u][end] .- exp(-pi^2 * 0.1) .* u0) < 5.0e-3
end

@testset "ODEProblem still requires scalar equations" begin
    sys, u, α = heat_slice_system(12)
    ssys = mtkcompile(sys; scalarize_arrays = false)
    @test_throws ArgumentError ODEProblem(ssys, [u => zeros(12)], (0.0, 0.1))
end

@testset "unsupported systems" begin
    @variables x(t)[1:3]
    @brownians a
    @test_throws ArgumentError mtkcompile(
        System([D(x) ~ -x .+ a], t; name = :noisy); scalarize_arrays = false
    )
    # a higher order array equation is not preserved, but scalarized and order-reduced
    ssys = mtkcompile(System([D(D(x)) ~ -x], t; name = :secondorder); scalarize_arrays = false)
    @test length(equations(ssys)) == 6
    @test length(unknowns(ssys)) == 6
end
