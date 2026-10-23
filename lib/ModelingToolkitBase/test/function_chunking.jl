using ModelingToolkitBase, SciMLBase, ForwardDiff
using ModelingToolkitBase: t_nounits as t, D_nounits as D, chunk_function_expr
using Test

@testset "chunk_function_expr" begin
    # `a_i` form a dependency chain so that every helper needs the previous one's results.
    stmts = Any[:(local k = 3), :(local u1 = x[1])]
    for i in 1:60
        push!(stmts, :(local $(Symbol(:a, i)) = $(i == 1 ? :u1 : Symbol(:a, i - 1)) * k + x[2]))
        i % 7 == 0 && push!(stmts, :(out[$(i ÷ 7)] = $(Symbol(:a, i)) + $(Symbol(:a, i - 3))))
    end
    push!(stmts, :(out))
    fex = Expr(:function, Expr(:tuple, :out, :x), Expr(:block, stmts...))
    expected = eval(fex)

    chunked = chunk_function_expr(fex, 10, 10)
    @test chunked != fex
    @test occursin("__mtk_chunk_", string(chunked))
    out_a, out_b = zeros(8), zeros(8)
    x = [0.5, 0.25]
    f = eval(chunked)
    @test Base.invokelatest(f, out_a, x) === out_a
    @test Base.invokelatest(expected, out_b, x) === out_b
    @test out_a == out_b
    @test (@allocated Base.invokelatest(f, out_a, x)) <= (@allocated Base.invokelatest(expected, out_b, x)) + 64

    # Below the threshold the expression is unchanged.
    @test chunk_function_expr(fex, 10, 10_000) == fex
    # Statements that rebind variables stay in the calling function.
    odd = :(
        function (x, y)
            z = x + y
            w = z * 2
            for i in 1:3
                z += i
            end
            z + w
        end
    )
    @test Base.invokelatest(eval(chunk_function_expr(odd, 1, 1)), 1.0, 2.0) ==
        Base.invokelatest(eval(odd), 1.0, 2.0)
    # No statement block to split.
    one_liner = :(
        function (x, y)
            x + y
        end
    )
    @test chunk_function_expr(one_liner, 1, 1) == one_liner
end

@testset "chunked generated functions match unchunked ones" begin
    n = 40
    @variables x(t)[1:n]
    @parameters c = 0.3
    xs = collect(x)
    eqs = [
        D(xs[i]) ~ c * (xs[mod1(i - 1, n)] - 2xs[i] + xs[mod1(i + 1, n)]) +
            sin(xs[i]) * exp(-xs[mod1(i + 2, n)]^2) / (1 + abs(xs[mod1(i - 2, n)]))
            for i in 1:n
    ]
    @named sys = System(eqs, t)
    sys = complete(sys)
    u0 = [xs[i] => 0.1i for i in 1:n]

    prob0 = ODEProblem(sys, u0, (0.0, 1.0); compiler_options = CompilerOptions(chunk_size = 0))
    for eval_expression in (false, true)
        prob = ODEProblem(
            sys, u0, (0.0, 1.0); eval_expression,
            compiler_options = CompilerOptions(chunk_size = 50)
        )
        du0, du = similar(prob0.u0), similar(prob0.u0)
        prob0.f(du0, prob0.u0, prob0.p, 0.3)
        Base.invokelatest(prob.f, du, prob.u0, prob.p, 0.3)
        @test du == du0
        @test Base.invokelatest(prob.f, prob.u0, prob.p, 0.3) == prob0.f(prob0.u0, prob0.p, 0.3)
        # Dual numbers flow through the helpers
        jac(p) = ForwardDiff.jacobian(u -> Base.invokelatest(p.f, u, p.p, 0.3), p.u0)
        @test jac(prob) == jac(prob0)
    end

    # The default threshold leaves small functions alone.
    expr = generate_rhs(sys; expression = Val(true))[2]
    @test chunk_function_expr(expr, 500, 2000) == expr
    @test occursin("__mtk_chunk_", string(chunk_function_expr(expr, 50, 50)))
end
