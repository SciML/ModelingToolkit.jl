using ModelingToolkitBase, SciMLBase, ForwardDiff
using ModelingToolkitBase: t_nounits as t, D_nounits as D, chunk_function_expr
using Test

# `getindex` with a side effect: reads must happen exactly once, where the original does them.
mutable struct CountingVector
    data::Vector{Float64}
    reads::Int
end
Base.getindex(v::CountingVector, i::Int) = (v.reads += 1; v.data[i])
Base.setindex!(v::CountingVector, x, i::Int) = (v.data[i] = x)

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

    # Helpers that exchange many values are not worth it.
    @test chunk_function_expr(fex, 10, 10, 0.0) == fex
    @test chunk_function_expr(fex, 10, 10, Inf) == chunked
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

@testset "chunk_function_expr is conservative" begin
    run(ex, args...) = Base.invokelatest(eval(ex), args...)
    # Constant-index reads must not be recomputed after the array was mutated.
    mutating = :(
        function (x)
            a = x[1]
            x[1] = 4
            b = a + 1
            b
        end
    )
    @test run(chunk_function_expr(mutating, 1, 1), [2.0]) == run(mutating, [2.0]) == 3.0
    # Aliases of an argument that are rebound before a store.
    alias_store = :(
        function (x, y)
            s = x[1]
            a = x
            a[1] = 4
            a = y
            b = s + 1
            b
        end
    )
    @test run(chunk_function_expr(alias_store, 1, 1), [2.0], [0.0]) ==
        run(alias_store, [2.0], [0.0]) == 3.0
    # Compound assignment to an argument.
    compound = :(
        function (x)
            a = x
            x += 1
            b = a + 1
            b
        end
    )
    @test run(chunk_function_expr(compound, 1, 1), 2.0) == run(compound, 2.0) == 3.0
    # An alias of a read must not be re-read after a later store.
    alias_read = :(
        function (x)
            a = x[1]
            b = a
            x[1] = 4
            c = b + 1
            c
        end
    )
    @test run(chunk_function_expr(alias_read, 1, 1), [2.0]) == run(alias_read, [2.0]) == 3.0
    # Broadcast calls obey the same callee rules.
    broadcast = :(
        function (x, y, fn)
            a = x[1]
            fn.(y)
            b = a + 1
            b
        end
    )
    bmut = v -> (x0[1] = 4.0; v)
    x0 = [2.0]
    @test run(chunk_function_expr(broadcast, 1, 1), x0, [1.0], bmut) == 3.0
    x0 = [2.0]
    @test run(broadcast, x0, [1.0], bmut) == 3.0
    # Reads of an argument happen exactly once.
    counted = :(
        function (x)
            a = x[1]
            b = a + 1
            c = b * 2
            d = a + c
            d
        end
    )
    cv1, cv2 = CountingVector([2.0], 0), CountingVector([2.0], 0)
    @test run(chunk_function_expr(counted, 1, 1), cv1) == run(counted, cv2)
    @test cv1.reads == cv2.reads == 1
    # Distinct names are not independent storage: the same array passed for both arguments.
    same_array = :(
        function (x, y)
            a = x[1]
            y[1] = 4
            b = a + 1
            b
        end
    )
    v = [2.0]
    w = [2.0]
    @test run(chunk_function_expr(same_array, 1, 1), v, v) == run(same_array, w, w) == 3.0
    # An embedded `Fix1` mutates through a Base function that does not end in `!`.
    for chunk in (false, true)
        arr = [2.0]
        fix = Expr(
            :function, Expr(:tuple, :x),
            Expr(:block, :(a = x[1]), Expr(:call, Base.Fix1(pushfirst!, arr), 5.0), :(b = a + 1), :b)
        )
        fix = chunk ? chunk_function_expr(fix, 1, 1) : fix
        @test run(fix, arr) == 3.0
    end
    # A parameter named like a whitelisted predicate is not that predicate.
    shadowed = :(
        function (x, isequal)
            a = x[1]
            isequal(x, 1)
            b = a + 1
            b
        end
    )
    mutator = (v, k) -> (v[1] = 4.0; false)
    @test run(chunk_function_expr(shadowed, 1, 1), [2.0], mutator) == run(shadowed, [2.0], mutator) == 3.0
    # Function objects embedded in the expression are not recognized as pure callees.
    mutate!(v) = (v[1] = 4.0; nothing)
    touch(v) = nothing
    for fn in (mutate!, touch)
        callee = Expr(
            :function, Expr(:tuple, :x),
            Expr(
                :block, :(a = x[1]), Expr(:call, fn, :x), :(b = a + 1), :b
            )
        )
        @test run(chunk_function_expr(callee, 1, 1), [2.0]) == run(callee, [2.0]) == 3.0
    end
    # A mutating callee given by name.
    named = :(
        function (x)
            a = x[1]
            mutate!(x)
            b = a + 1
            b
        end
    )
    @test chunk_function_expr(named, 1, 1) == named
    # Helper names never collide with existing bindings.
    clash = :(
        function (__mtk_chunk_0)
            a = __mtk_chunk_0 + 1
            b = a + 1
            b
        end
    )
    @test run(chunk_function_expr(clash, 1, 1), 2) == run(clash, 2) == 4
    # Rebinding a name.
    rebinding = :(
        function (x)
            a = x + 1
            a = a + 1
            a
        end
    )
    @test run(chunk_function_expr(rebinding, 1, 1), 1.0) == 3.0
    # Closures capture variables that the liveness analysis cannot see.
    closure = :(
        function (x)
            f = () -> a
            a = x + 2
            f()
        end
    )
    @test chunk_function_expr(closure, 1, 1) == closure
    # Empty blocks never throw.
    for body in (Expr(:block), Expr(:block, Expr(:block)), Expr(:block, Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(0), Expr(:block))))
        empty_fn = Expr(:function, Expr(:tuple, :x), body)
        @test chunk_function_expr(empty_fn, 1, 1) == empty_fn
    end
    # Bounds checks are only disabled in helpers if the original function disabled them.
    body = Any[:(local $(Symbol(:a, i)) = x[$i] + 1) for i in 1:20]
    push!(body, :(a20))
    checked = Expr(:function, Expr(:tuple, :x), Expr(:block, body...))
    chunked = chunk_function_expr(checked, 5, 5)
    @test chunked != checked
    @test !occursin("@inbounds", string(chunked))
    @test_throws BoundsError run(chunked, zeros(3))
    unchecked = Expr(
        :function, Expr(:tuple, :x),
        Expr(:macrocall, Symbol("@inbounds"), LineNumberNode(0), Expr(:block, body...))
    )
    @test occursin("@inbounds", string(chunk_function_expr(unchecked, 5, 5)))
    # The 3-argument constructor still exists.
    @test CompilerOptions(1, 2, 3).chunk_size == -1
    @test CompilerOptions(-1, -1, -1) == CompilerOptions()
    @test CompilerOptions(1, 2, true).infer == 1
    @test CompilerOptions(Int32(1), Int32(2), Int32(0)) == CompilerOptions(1, 2, 0)
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
