using Test, ModelingToolkit, SciMLBase, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

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
    for (makeprob, alg) in ((chain_ode, Tsit5()), (chain_dae, FBDF()))
        probs = [makeprob(n) for n in (5, 9)]
        @test typeof(probs[1].f.initialization_data) === typeof(probs[2].f.initialization_data)
        for (n, prob) in zip((5, 9), probs)
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
    end
end
