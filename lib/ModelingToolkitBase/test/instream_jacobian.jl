function stream_finite_difference(f, x; step = 1.0e-5)
    return reduce(
        hcat, map(eachindex(x)) do i
            xp, xm = copy(x), copy(x)
            xp[i] += step
            xm[i] -= step
            (f(xp) - f(xm)) / (2 * step)
        end
    )
end

@testset "instream runtime derivatives" begin
    for (ni, no) in [(2, 0), (3, 0), (1, 1), (2, 1)]
        n = ni + no
        @variables v[1:(2 * n)]
        vars = collect(v)
        flows = [1:ni; (2 * ni + 1):(2 * ni + no)]
        streams = [(ni + 1):(2 * ni); (2 * ni + no + 1):(2 * n)]
        term = SU.term(
            ModelingToolkitBase.instream_rt, Val(ni), Val(no), vars...; type = Real
        )
        grad = Symbolics.gradient(Num(term), vars)
        gradfun = Symbolics.build_function(grad, vars; expression = Val(false))[1]
        runtime(x) = ModelingToolkitBase.instream_rt(Val(ni), Val(no), x...)
        for flow in [-1.0, -5.0e-5, -1.0e-4 / n, -1.0e-5, 0.0, 1.0]
            x = zeros(2 * n)
            x[flows] .= flow
            x[flows[(ni + 1):end]] .*= -1
            x[streams] .= 2:(n + 1)
            actual = gradfun(x)
            @test all(isfinite, actual)
            # Central differences converge linearly at the C¹ smoothing cutoff.
            expected = vec(stream_finite_difference(runtime, BigFloat.(x); step = big"1e-20"))
            @test actual ≈ expected atol = 1.0e-8 rtol = 1.0e-10
        end
        x = zeros(2 * n)
        x[flows] .= -1.0
        x[streams] .= 2:(n + 1)
        x[first(flows)] = 0.0
        @test all(isfinite, gradfun(x))
        @test iszero(gradfun(x)[first(flows)])
    end
end

if @isdefined(ModelingToolkit)
    @connector function JacobianStreamPort(; name)
        @variables p(t) m(t) [connect = Flow] h(t) [connect = Stream]
        return System(Equation[], t, [p, m, h], []; name)
    end

    function JacobianStreamSource(; name, rate, start, variable_flow = false)
        @named port = JacobianStreamPort()
        @variables x(t) = start
        eqs = [port.h ~ x, D(x) ~ -x]
        if variable_flow
            @variables q(t) = rate
            append!(eqs, [port.m ~ -q, D(q) ~ -1])
        else
            push!(eqs, port.m ~ -rate)
        end
        return System(eqs, t; name, systems = [port])
    end

    function JacobianStreamReceiver(; name)
        @named port = JacobianStreamPort()
        @variables z(t) = 3.0
        return System(
            [port.p ~ 1.0, port.h ~ z, D(z) ~ instream(port.h) - z], t;
            name, systems = [port]
        )
    end

    @testset "Multi-source stream Jacobian" begin
        for variable_flow in (false, true)
            @named a = JacobianStreamSource(rate = 1.0, start = 2.0, variable_flow = variable_flow)
            @named b = JacobianStreamSource(rate = 2.0, start = 4.0)
            @named c = JacobianStreamReceiver()
            sys = mtkcompile(
                System(
                    [connect(a.port, b.port, c.port)], t;
                    name = :mixer, systems = [a, b, c]
                )
            )
            prob = ODEProblem(sys, [], (0.0, 2.0); jac = true)
            idx(v) = findfirst(isequal(v), unknowns(sys))
            for rate in (variable_flow ? (1.0, -1.0) : (1.0,))
                u = copy(prob.u0)
                if variable_flow
                    u[idx(sys.a.q)] = rate
                end
                J = zeros(length(u), length(u))
                prob.f.jac(J, u, prob.p, 1 - rate)
                fd = stream_finite_difference(x -> prob.f(x, prob.p, 1 - rate), u)
                @test J ≈ fd atol = 1.0e-8 rtol = 0
                if !variable_flow
                    @test J[idx(sys.c.z), idx(sys.a.x)] ≈ 1 / 3
                    @test J[idx(sys.c.z), idx(sys.b.x)] ≈ 2 / 3
                end
            end
        end
    end
end
