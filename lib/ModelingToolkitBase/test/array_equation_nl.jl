using ModelingToolkitBase, Test
using ModelingToolkitBase: complete, mtkcompile, NonlinearSystem, is_time_dependent
using ModelingToolkitBase: array_residual_maker, count_equation_rows,
    equation_row_count, generate_rhs, generate_diffusion_function, eval_or_rgf
using ModelingToolkitBase: _iszero, SymbolicT
using ModelingToolkitBase: array_derivative_is_atomic, array_derivative_expansion,
    collect_applied_operators, has_array_equations
using Symbolics
using SciMLBase
using SciMLStructures: replace, Tunable
using NonlinearSolve
using ForwardDiff
using StaticArrays
import SymbolicUtils as SU

# 1D Laplace interior as one array equation.
function laplace_array_system(n; pdebase_lhs_zeros = true, compile = complete)
    @variables u[1:n]
    dx = 1 / (n - 1)
    lap = (u[1:(n - 2)] .- 2 .* u[2:(n - 1)] .+ u[3:n]) ./ dx^2
    interior = if pdebase_lhs_zeros
        zeros(n - 2) ~ broadcast(-, lap)
    else
        lap ~ zeros(n - 2)
    end
    eqs = [interior, u[1] ~ 0.0, u[n] ~ 1.0]
    @named sys = System(eqs, collect(u), [])
    return compile(sys), u
end

function residual_vector(eqs)
    return SymbolicT[_iszero(eq.lhs) ? eq.rhs : eq.rhs - eq.lhs for eq in eqs]
end

# Time-dependent heat residual for steady-state conversion.
function heat_array_nl_system(n; lhs_derivative = false, compile = complete)
    @independent_variables t
    @variables u(t)[1:n]
    D = Differential(t)
    dx = 1 / (n - 1)
    lap = (u[1:(n - 2)] .- 2 .* u[2:(n - 1)] .+ u[3:n]) ./ dx^2
    interior = if lhs_derivative
        D(u[2:(n - 1)]) ~ lap
    else
        broadcast(-, D(u[2:(n - 1)]), lap) ~ zeros(n - 2)
    end
    eqs = [interior, u[1] ~ 0.0, u[n] ~ 1.0]
    @named sys = System(eqs, t, collect(u), [])
    return compile(sys), u, t, D
end

function unscalarized_unknown_system(n)
    @variables u[1:n]
    @named sys = System([zeros(n) ~ u .- 1], [u], [])
    return complete(sys), u
end

function nonsquare_array_system()
    @variables x y
    eqs = Equation[zeros(3) ~ [x - 1, y - 2, x + y - 3]]
    @named sys = System(eqs, [x, y], [])
    return complete(sys), x, y
end

@testset "equation_row_count accepts literal zero arrays" begin
    @variables x[1:5]
    @test equation_row_count(zeros(5) ~ x) == 5
    @variables y[1:3, 1:3]
    @test equation_row_count(zeros(3, 3) ~ y) == 9
    @variables z
    @test equation_row_count(0 ~ z) == 1
    @test count_equation_rows([zeros(5) ~ x, z ~ 0.0]) == 6
end

@testset "array equations reach NonlinearProblem" begin
    n = 11
    sys, u = laplace_array_system(n)
    op = [u[i] => 0.5 for i in 1:n]

    prob = NonlinearProblem(sys, op)

    @test length(prob.u0) == n
    @test prob.u0 isa Vector{Float64}
    @test prob.f.resid_prototype === nothing

    out = zeros(n)
    prob.f(out, prob.u0, prob.p)
    @test all(isfinite, out)
    @test any(!iszero, out)
end

@testset "both residual equation forms construct" begin
    n = 11
    for pdebase_lhs_zeros in (true, false)
        sys, u = laplace_array_system(n; pdebase_lhs_zeros)
        op = [u[i] => 0.5 for i in 1:n]
        prob = NonlinearProblem(sys, op)
        @test length(prob.u0) == n
        out = zeros(n)
        prob.f(out, prob.u0, prob.p)
        @test all(isfinite, out)
    end
end

@testset "symbolic residual assembly stays one ArrayMaker as n grows" begin
    n21 = 21
    n101 = 101
    sys21, _ = laplace_array_system(n21)
    sys101, _ = laplace_array_system(n101)

    @test length(equations(sys21)) == length(equations(sys101)) == 3
    @test count_equation_rows(equations(sys21)) == n21
    @test count_equation_rows(equations(sys101)) == n101

    made21 = array_residual_maker(residual_vector(equations(sys21)))
    made101 = array_residual_maker(residual_vector(equations(sys101)))
    @test !(made21 isa AbstractVector)
    @test !(made101 isa AbstractVector)
    @test prod(length, SU.shape(made21)) == n21
    @test prod(length, SU.shape(made101)) == n101

    regs21 = SU.unwrap_const(first(arguments(made21)))
    regs101 = SU.unwrap_const(first(arguments(made101)))
    @test length(regs21) == length(regs101) == 3
    @test only.(regs21) == [1:(n21 - 2), (n21 - 1):(n21 - 1), n21:n21]
    @test only.(regs101) == [1:(n101 - 2), (n101 - 1):(n101 - 1), n101:n101]

    iip21 = generate_rhs(sys21; expression = Val{true})[2]
    iip101 = generate_rhs(sys101; expression = Val{true})[2]
    s21 = string(iip21)
    s101 = string(iip101)
    @test count(==('\n'), s21) == count(==('\n'), s101)
end

@testset "array-equation Nonlinear solves to the analytic solution" begin
    n = 21
    sys, u = laplace_array_system(n)
    op = [u[i] => 0.5 for i in 1:n]
    prob = NonlinearProblem(sys, op)
    sol = solve(prob, NewtonRaphson())
    @test SciMLBase.successful_retcode(sol)
    exact = range(0.0, 1.0, length = n)
    @test maximum(abs.(sol.u .- exact)) < 1.0e-10
end

@testset "array equations over a 2D slice keep their shape" begin
    n = 6
    @variables w[1:n, 1:n]
    dx = 1 / (n - 1)
    inner = 2:(n - 1)
    lap = (
        w[1:(n - 2), inner] .+ w[3:n, inner] .+ w[inner, 1:(n - 2)] .+
            w[inner, 3:n] .- 4 .* w[inner, inner]
    ) ./ dx^2
    eqs = Equation[zeros(n - 2, n - 2) ~ broadcast(-, lap)]
    for i in 1:n
        x = (i - 1) * dx
        push!(eqs, w[i, 1] ~ x)
        push!(eqs, w[i, n] ~ x)
    end
    for j in inner
        push!(eqs, w[1, j] ~ 0.0)
        push!(eqs, w[n, j] ~ 1.0)
    end
    @named sys2d = System(eqs, vec(collect(w)), [])
    sys2d = complete(sys2d)

    op = [w[i, j] => 0.5 for i in 1:n, j in 1:n] |> vec
    prob = NonlinearProblem(sys2d, op)
    @test length(prob.u0) == n * n
    @test prob.f.resid_prototype === nothing
    out = zeros(n * n)
    prob.f(out, prob.u0, prob.p)
    @test all(isfinite, out)

    sol = solve(prob, NewtonRaphson())
    @test SciMLBase.successful_retcode(sol)
    exact = [(i - 1) * dx for i in 1:n, j in 1:n]
    @test maximum(abs.(reshape(sol.u, n, n) .- exact)) < 1.0e-8
end

@testset "other problem types still require scalarized equations" begin
    n = 11
    sys, u, _, _ = heat_array_nl_system(n)
    op = [u[i] => 0.0 for i in 1:n]
    err = @test_throws ArgumentError ODEProblem(
        sys, op, (0.0, 0.1); build_initializeprob = false
    )
    @test occursin("array equations", err.value.msg)
end

@testset "scalar 0 ~ f systems do not assemble an ArrayMaker" begin
    @variables x y z
    @parameters σ ρ β
    eqs = [
        0 ~ σ * (y - x),
        0 ~ x * (ρ - z) - y,
        0 ~ x * y - β * z,
    ]
    @named ns = System(eqs, [x, y, z], [σ, ρ, β])
    ns = complete(ns)
    made = array_residual_maker(residual_vector(equations(ns)))
    @test made isa AbstractVector
    prob = NonlinearProblem(
        ns, [x => 1.0, y => 1.0, z => 1.0, σ => 1.0, ρ => 1.0, β => 1.0]
    )
    @test prob.f.resid_prototype === nothing
    sol = solve(prob, NewtonRaphson())
    @test SciMLBase.successful_retcode(sol)
    @test sol[x] ≈ sol[y]
end

@testset "array residual matches the mtkcompile solution" begin
    n = 21
    sys, u = laplace_array_system(n)
    sys_sc, _ = laplace_array_system(n; compile = mtkcompile)
    op = [u[i] => 0.5 for i in 1:n]
    sol_arr = solve(NonlinearProblem(sys, op), NewtonRaphson())
    sol_sc = solve(NonlinearProblem(sys_sc, op), NewtonRaphson())
    @test SciMLBase.successful_retcode(sol_arr)
    @test SciMLBase.successful_retcode(sol_sc)
    @test sol_arr[u] ≈ sol_sc[u] rtol = 1.0e-8
end

@testset "array derivative expansion preserves shape" begin
    @independent_variables t
    @variables u(t)[1:5] w(t)[1:4, 1:4]
    D = Differential(t)
    term1 = unwrap(D(u[2:4]))
    @test array_derivative_is_atomic(term1)
    exp1 = array_derivative_expansion(term1)
    @test SU.shape(exp1) == SU.shape(unwrap(u[2:4]))
    @test length(SU.shape(exp1)) == 1

    term2 = unwrap(D(w[2:3, 2:3]))
    @test array_derivative_is_atomic(term2)
    exp2 = array_derivative_expansion(term2)
    @test SU.shape(exp2) == SU.shape(unwrap(w[2:3, 2:3]))
    @test length(SU.shape(exp2)) == 2
    @test !array_derivative_is_atomic(unwrap(D(u[1])))
end

@testset "time-dependent D(slice) becomes a finite steady-state residual" begin
    n = 21
    exact = range(0.0, 1.0, length = n)
    for lhs_derivative in (false, true)
        sys, u, _, _ = heat_array_nl_system(n; lhs_derivative)
        nlsys = NonlinearSystem(sys)
        @test !is_time_dependent(nlsys)
        for eq in equations(nlsys)
            @test isempty(collect_applied_operators(eq, Differential))
        end
        @test has_array_equations(equations(nlsys))

        op = [u[i] => 0.5 for i in 1:n]
        prob = NonlinearProblem(sys, op)
        out = zeros(n)
        prob.f(out, prob.u0, prob.p)
        @test all(isfinite, out)
        @test any(!iszero, out)

        sol = solve(prob, NewtonRaphson())
        @test SciMLBase.successful_retcode(sol)
        @test maximum(abs.(sol.u .- exact)) < 1.0e-10
    end
end

@testset "2D D(matrix_slice) steady-state keeps its shape" begin
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
        x = (i - 1) * dx
        push!(eqs, w[i, 1] ~ x)
        push!(eqs, w[i, n] ~ x)
    end
    for j in inner
        push!(eqs, w[1, j] ~ 0.0)
        push!(eqs, w[n, j] ~ 1.0)
    end
    @named sys2d = System(eqs, t, vec(collect(w)), [])
    sys2d = complete(sys2d)
    nlsys = NonlinearSystem(sys2d)
    for eq in equations(nlsys)
        @test isempty(collect_applied_operators(eq, Differential))
    end
    interior = first(equations(nlsys))
    @test SU.is_array_shape(SU.shape(interior.lhs))
    @test length(SU.shape(interior.lhs)) == 2

    op = [w[i, j] => 0.5 for i in 1:n, j in 1:n] |> vec
    prob = NonlinearProblem(sys2d, op)
    out = zeros(n * n)
    prob.f(out, prob.u0, prob.p)
    @test all(isfinite, out)

    sol = solve(prob, NewtonRaphson())
    @test SciMLBase.successful_retcode(sol)
    exact = [(i - 1) * dx for i in 1:n, j in 1:n]
    @test maximum(abs.(reshape(sol.u, n, n) .- exact)) < 1.0e-8
end

@testset "array unknown guard covers every nonlinear entry" begin
    sys, u = unscalarized_unknown_system(4)
    @test any(Symbolics.isarraysymbolic, unknowns(sys))
    op = [u => ones(4)]

    err = @test_throws ArgumentError NonlinearProblem(sys, op)
    @test occursin("array unknowns", err.value.msg)
    @test occursin("collect(u)", err.value.msg)

    err = @test_throws ArgumentError NonlinearLeastSquaresProblem(sys, op)
    @test occursin("array unknowns", err.value.msg)

    err = @test_throws ArgumentError NonlinearFunction(sys)
    @test occursin("array unknowns", err.value.msg)

    # scalarized unknowns remain the happy path
    sys_sc, u_sc = laplace_array_system(5)
    @test_nowarn NonlinearProblem(sys_sc, [u_sc[i] => 0.5 for i in 1:5])
end

@testset "jac and sparse fail fast on array residuals" begin
    n = 11
    sys, u = laplace_array_system(n)
    op = [u[i] => 0.5 for i in 1:n]
    for kwargs in ((; jac = true), (; sparse = true), (; jac = true, sparse = true))
        err = @test_throws ArgumentError NonlinearProblem(sys, op; kwargs...)
        @test occursin("jac = true", err.value.msg)
        @test occursin("sparse = true", err.value.msg)
        @test occursin("mtkcompile", err.value.msg)
    end
    err = @test_throws ArgumentError NonlinearProblem{false}(sys, op; jac = true)
    @test occursin("array residuals", lowercase(err.value.msg))
    err = @test_throws ArgumentError NonlinearFunction(sys; jac = true)
    @test occursin("array residuals", lowercase(err.value.msg))
    err = @test_throws ArgumentError NonlinearFunction(sys; sparse = true, expression = Val{true})
    @test occursin("array residuals", lowercase(err.value.msg))

    sys_sc, _ = laplace_array_system(n; compile = mtkcompile)
    prob_jac = NonlinearProblem(sys_sc, op; jac = true)
    J = prob_jac.f.jac(prob_jac.u0, prob_jac.p)
    @test size(J) == (length(prob_jac.u0), length(prob_jac.u0))
    @test all(isfinite, J)
    prob_sparse = NonlinearProblem(sys_sc, op; sparse = true)
    @test prob_sparse.f.jac_prototype !== nothing
end

@testset "non-square array residual contract" begin
    sys, x, y = nonsquare_array_system()
    @test count_equation_rows(equations(sys)) == 3
    @test length(unknowns(sys)) == 2
    op = [x => 0.0, y => 0.0]

    prob = NonlinearLeastSquaresProblem(sys, op)
    @test length(prob.u0) == 2
    @test prob.f.resid_prototype !== nothing
    @test length(prob.f.resid_prototype) == 3

    out = zeros(3)
    prob.f(out, prob.u0, prob.p)
    @test length(out) == 3
    @test all(isfinite, out)

    proboop = NonlinearLeastSquaresProblem{false}(sys, op)
    resid = proboop.f(proboop.u0, proboop.p)
    @test length(resid) == 3
    J = ForwardDiff.jacobian(u -> collect(proboop.f(u, proboop.p)), proboop.u0)
    @test size(J) == (3, 2)
    @test all(isfinite, J)

    sol = solve(prob)
    @test SciMLBase.successful_retcode(sol)
    @test sol[x] ≈ 1
    @test sol[y] ≈ 2

    probsa = NonlinearLeastSquaresProblem{false}(sys, op; u0_constructor = splat(SVector))
    @test probsa.u0 isa SVector
    @test length(probsa.u0) == 2
    @test length(probsa.f.resid_prototype) == 3
    @test length(probsa.f.resid_prototype) != length(probsa.u0)
    resid_sa = probsa.f(probsa.u0, probsa.p)
    @test length(resid_sa) == 3
    @test all(isfinite, collect(resid_sa))
    Jsa = ForwardDiff.jacobian(u -> collect(probsa.f(SVector{2}(u), probsa.p)), collect(probsa.u0))
    @test size(Jsa) == (3, 2)
    @test all(isfinite, Jsa)

    probfw = NonlinearLeastSquaresProblem{true, SciMLBase.FunctionWrapperSpecialize}(sys, op)
    @test SciMLBase.specialization(probfw.f) === SciMLBase.FunctionWrapperSpecialize
    @test length(probfw.f.resid_prototype) == 3
    @test length(probfw.f.resid_prototype) != length(probfw.u0)
    outfw = zeros(3)
    probfw.f(outfw, probfw.u0, probfw.p)
    @test length(outfw) == 3
    @test all(isfinite, outfw)
end

@testset "OOP and expression entry points" begin
    n = 11
    sys, u = laplace_array_system(n)
    op = [u[i] => 0.5 for i in 1:n]

    proboop = NonlinearProblem{false}(sys, op)
    resid = proboop.f(proboop.u0, proboop.p)
    @test length(resid) == n
    @test all(isfinite, resid)

    oop_expr, iip_expr = generate_rhs(sys; expression = Val{true})
    @test occursin("similar_for_residual", string(oop_expr))
    @test !occursin("similar_for_residual", string(iip_expr))
    f_oop = eval_or_rgf(oop_expr)
    f_iip = eval_or_rgf(iip_expr)
    resid_expr = f_oop(proboop.u0, proboop.p)
    @test resid_expr ≈ resid
    out_expr = zeros(n)
    f_iip(out_expr, proboop.u0, proboop.p)
    @test out_expr ≈ resid

    J = ForwardDiff.jacobian(u -> collect(proboop.f(u, proboop.p)), proboop.u0)
    @test size(J) == (n, n)
    @test all(isfinite, J)

    sol = solve(proboop, NewtonRaphson())
    @test SciMLBase.successful_retcode(sol)
    @test maximum(abs.(sol.u .- range(0.0, 1.0, length = n))) < 1.0e-10

    probeval = eval(NonlinearProblem(sys, op; expression = Val{true}))
    @test probeval isa NonlinearProblem
    out_eval = zeros(n)
    probeval.f(out_eval, probeval.u0, probeval.p)
    @test all(isfinite, out_eval)
    feval = eval(NonlinearFunction(sys; expression = Val{true}))
    @test feval isa NonlinearFunction
    feval(out_eval, probeval.u0, probeval.p)
    @test all(isfinite, out_eval)
end

@testset "parameterized stencil residual" begin
    n = 11
    @variables u[1:n]
    @parameters dx
    lap = (u[1:(n - 2)] .- 2 .* u[2:(n - 1)] .+ u[3:n]) ./ dx^2
    eqs = [zeros(n - 2) ~ broadcast(-, lap), u[1] ~ 0.0, u[n] ~ 1.0]
    @named sys = System(eqs, collect(u), [dx])
    sys = complete(sys)
    op = vcat([u[i] => 0.5 for i in 1:n], [dx => 1 / (n - 1)])
    prob = NonlinearProblem(sys, op)
    out = zeros(n)
    prob.f(out, prob.u0, prob.p)
    @test all(isfinite, out)
    sol = solve(prob, NewtonRaphson())
    @test SciMLBase.successful_retcode(sol)
    @test maximum(abs.(sol.u .- range(0.0, 1.0, length = n))) < 1.0e-10
end

@testset "assertions apply to each row of an array residual" begin
    n = 11
    @variables u[1:n]
    dx = 1 / (n - 1)
    lap = (u[1:(n - 2)] .- 2 .* u[2:(n - 1)] .+ u[3:n]) ./ dx^2
    # The array equation is last so the assertion is broadcast onto its rows.
    eqs = [u[1] ~ 0.0, u[n] ~ 1.0, zeros(n - 2) ~ broadcast(-, lap)]
    op = [u[i] => 0.5 for i in 1:n]

    @named sys_ok = System(eqs, collect(u), []; assertions = [(u[2] > 0) => "u[2] positive"])
    sys_ok = complete(sys_ok)
    prob_ok = NonlinearProblem(sys_ok, op)
    out_ok = zeros(n)
    prob_ok.f(out_ok, prob_ok.u0, prob_ok.p)
    @test all(isfinite, out_ok)

    @named sys_bad = System(eqs, collect(u), []; assertions = [(u[2] < 0) => "u[2] negative"])
    sys_bad = complete(sys_bad)
    prob_bad = NonlinearProblem(sys_bad, op)
    out_bad = zeros(n)
    prob_bad.f(out_bad, prob_bad.u0, prob_bad.p)
    @test any(isnan, out_bad)
end

@testset "n=101 evaluates a real IIP residual" begin
    n = 101
    sys, u = laplace_array_system(n)
    @test count_equation_rows(equations(sys)) == n
    op = [u[i] => 0.5 for i in 1:n]
    prob = NonlinearProblem(sys, op)
    out = zeros(n)
    prob.f(out, prob.u0, prob.p)
    @test length(out) == n
    @test all(isfinite, out)
    @test any(!iszero, out)
end

@testset "SDEProblem rejects array equations" begin
    n = 5
    @independent_variables t
    @variables u(t)[1:n]
    D = Differential(t)
    eqs = [D(u[2:(n - 1)]) ~ -u[2:(n - 1)], u[1] ~ 0.0, u[n] ~ 0.0]
    @named sys = System(eqs, t, collect(u), []; noise_eqs = [0.01, 0.01, 0.01])
    sys = complete(sys)
    op = [u[i] => 0.0 for i in 1:n]
    err = @test_throws ArgumentError SDEProblem(
        sys, op, (0.0, 0.1); build_initializeprob = false
    )
    @test occursin("array equations", err.value.msg)
end

@testset "OOP nonlinear residual promotes Dual inputs" begin
    n = 11
    @variables u[1:n]
    @parameters dx
    lap = (u[1:(n - 2)] .- 2 .* u[2:(n - 1)] .+ u[3:n]) ./ dx^2
    eqs = [zeros(n - 2) ~ broadcast(-, lap), u[1] ~ 0.0, u[n] ~ 1.0]
    @named sys = System(eqs, collect(u), [dx])
    sys = complete(sys)
    dxval = 1 / (n - 1)
    op = vcat(
        [u[i] => 0.5 + 0.1 * sinpi((i - 1) / (n - 1)) for i in 1:n],
        [dx => dxval]
    )
    proboop = NonlinearProblem{false}(sys, op)
    @test eltype(proboop.u0) === Float64

    p_dual = replace(Tunable(), proboop.p, [ForwardDiff.Dual(dxval, 1.0)])
    resid_p = proboop.f(proboop.u0, p_dual)
    @test eltype(resid_p) <: ForwardDiff.Dual
    @test length(resid_p) == n
    @test all(isfinite, ForwardDiff.value.(resid_p))
    Jp = ForwardDiff.jacobian(
        θs -> collect(
            proboop.f(proboop.u0, replace(Tunable(), proboop.p, θs))
        ),
        [dxval]
    )
    @test size(Jp) == (n, 1)
    @test all(isfinite, Jp)
    @test any(!iszero, Jp)

    sys_ns = complete(
        System(eqs, collect(u), [dx]; name = :sys_ns);
        split = false
    )
    proboop_ns = NonlinearProblem{false}(sys_ns, op)
    resid_ns = proboop_ns.f(proboop_ns.u0, [ForwardDiff.Dual(dxval, 1.0)])
    @test eltype(resid_ns) <: ForwardDiff.Dual
    @test length(resid_ns) == n
    @test all(isfinite, ForwardDiff.value.(resid_ns))

    n21 = 21
    sys21, u21 = laplace_array_system(n21)
    op21 = [u21[i] => 0.5 for i in 1:n21]
    proboop21 = NonlinearProblem{false}(sys21, op21)
    J21 = ForwardDiff.jacobian(u -> collect(proboop21.f(u, proboop21.p)), proboop21.u0)
    @test size(J21) == (n21, n21)
    @test all(isfinite, J21)
end

@testset "scalar ODE/SDE/DAE codegen stays scalar" begin
    @independent_variables t
    @variables x(t) y(t)
    D = Differential(t)
    @named odesys = System([D(x) ~ -x, D(y) ~ -2y], t)
    odesys = complete(odesys)
    iip_ode = generate_rhs(odesys; expression = Val{true})[2]
    @test !occursin("ArrayMaker", string(iip_ode))
    @test !occursin("similar_for_residual", string(iip_ode))

    @named daesys = System([D(x) ~ -x, 0 ~ x + y], t)
    daesys = complete(daesys)
    iip_dae = generate_rhs(daesys; implicit_dae = true, expression = Val{true})[2]
    @test !occursin("ArrayMaker", string(iip_dae))
    @test !occursin("similar_for_residual", string(iip_dae))

    @named sdesys = System([D(x) ~ -x, D(y) ~ -2y], t; noise_eqs = [0.01, 0.01])
    sdesys = complete(sdesys)
    iip_sde = generate_rhs(sdesys; expression = Val{true})[2]
    @test !occursin("ArrayMaker", string(iip_sde))
    @test !occursin(
        "ArrayMaker", string(generate_diffusion_function(sdesys; expression = Val{true}))
    )

    @variables a b
    @named nlsys = System([0 ~ a + b, 0 ~ a - b], [a, b], [])
    nlsys = complete(nlsys)
    iip_nl = generate_rhs(nlsys; expression = Val{true})[2]
    @test !occursin("ArrayMaker", string(iip_nl))
    @test !occursin("similar_for_residual", string(iip_nl))
end
