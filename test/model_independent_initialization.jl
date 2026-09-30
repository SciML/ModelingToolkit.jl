using Test, ModelingToolkit, SciMLBase, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: ModelingToolkitBase

function chain_ode(n, spec = SciMLBase.AutoDespecialize)
    @variables x(t)[1:n]
    @parameters a[1:n]
    eqs = [D(x[i]) ~ -a[i] * x[i] + (i == 1 ? 0 : x[i - 1]) for i in 1:n]
    sys = mtkcompile(System(eqs, t; name = :chain_ode))
    return ODEProblem{true, spec}(
        sys, [
            [sys.x[i] => 1.0 / i for i in 1:n];
            [sys.a[i] => 0.5 + i for i in 1:n]
        ], (0.0, 1.0)
    )
end

function chain_dae(n, spec = SciMLBase.AutoDespecialize)
    @variables x(t)[1:n] y(t)[1:n]
    @parameters a[1:n]
    eqs = Equation[]
    for i in 1:n
        push!(eqs, D(x[i]) ~ -a[i] * x[i] + y[i])
        push!(eqs, 0 ~ y[i]^3 + y[i] - x[i] - (i == 1 ? 0 : y[i - 1]))
    end
    sys = mtkcompile(System(eqs, t; name = :chain_dae))
    return ODEProblem{true, spec}(
        sys, [
            [sys.x[i] => 1.0 / i for i in 1:n];
            [sys.a[i] => 0.5 + i for i in 1:n]
        ], (0.0, 1.0);
        guesses = [sys.y[i] => 0.5 for i in 1:n]
    )
end

@testset "Model-independent initialization data" begin
    # DAE sizes both exceed 5 SCC blocks so type equality does not depend on
    # the ≤5-block tuple container (left to ModelingToolkit.jl#5202).
    for (makeprob, alg, ns) in (
            (chain_ode, Tsit5(), (6, 9)),
            (chain_dae, FBDF(), (6, 9)),
        )
        probs = [makeprob(n) for n in ns]
        @test typeof(probs[1].f.initialization_data) === typeof(probs[2].f.initialization_data)
        for (n, prob) in zip(ns, probs)
            fullprob = makeprob(n, SciMLBase.FullSpecialize)
            sol = solve(prob, alg)
            fullsol = solve(fullprob, alg)
            @test SciMLBase.successful_retcode(sol)
            @test SciMLBase.successful_retcode(fullsol)
            @test sol.u[1] ≈ fullsol.u[1]
            @test sol.u[end] ≈ fullsol.u[end]
            remade = remake(prob; u0 = copy(sol.u[1]))
            remadesol = solve(remade, alg)
            @test SciMLBase.successful_retcode(remadesol)
            @test remadesol.u[end] ≈ sol.u[end]
        end
        autospec = makeprob(ns[1], SciMLBase.AutoSpecialize)
        @test !(
            autospec.f.initialization_data.initializeprobpmap isa
                ModelingToolkitBase.InitializeprobParameterMap
        )
    end
end

@testset "CopyParamsByTemplate records non-tunable array-parameter indices" begin
    @parameters a[1:3]
    @parameters c[1:2, 1:2] [tunable = false]
    @variables x(t)[1:3]
    eqs = [
        D(x[1]) ~ -a[1] * x[1] + c[1, 1] + c[2, 2],
        D(x[2]) ~ -a[2] * x[2] + c[1, 2],
        D(x[3]) ~ -a[3] * x[3] + c[2, 1],
    ]
    sys = mtkcompile(System(eqs, t; name = :nontunable_template))
    syms = [
        sys.a[1], sys.c[1, 1], sys.a[2], sys.c[1, 2],
        sys.a[3], sys.c[2, 1], sys.c[2, 2],
    ]
    getter = ModelingToolkitBase.CopyParamsByTemplate(
        sys,
        ModelingToolkitBase.SymbolicT[ModelingToolkitBase.unwrap(s) for s in syms];
        eval_expression = false, eval_module = @__MODULE__
    )
    @test any(t -> t isa ModelingToolkitBase.ParameterIndex && t.idx isa Tuple, getter.template)
    op = vcat(
        [sys.x[i] => 1.0 for i in 1:3],
        [sys.a[i] => Float64(i) for i in 1:3],
        [sys.c[i, j] => 10.0 * i + j for i in 1:2 for j in 1:2],
    )
    prob = ODEProblem(sys, op, (0.0, 1.0))
    @test getter(prob) ≈ [1.0, 11.0, 2.0, 12.0, 3.0, 21.0, 22.0]
end
