using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Zygote
using SymbolicIndexingInterface
using SciMLStructures
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqNonlinearSolve
using NonlinearSolve
using SciMLSensitivity
using ForwardDiff
using StableRNGs
using ChainRulesCore
using ChainRulesCore: NoTangent
using ChainRulesTestUtils: test_rrule, rand_tangent

@variables x(t)[1:3] y(t)
@parameters p[1:3, 1:3] q
eqs = [D(x) ~ p * x
       D(y) ~ sum(p) + q * y]
u0 = [x => zeros(3),
    y => 1.0]
ps = [p => zeros(3, 3),
    q => 1.0]
tspan = (0.0, 10.0)
@mtkcompile sys = System(eqs, t)
prob = ODEProblem(sys, [u0; ps], tspan)
sol = solve(prob, Tsit5())

mtkparams = parameter_values(prob)
new_p = rand(14)
gs = gradient(new_p) do new_p
    new_params = SciMLStructures.replace(SciMLStructures.Tunable(), mtkparams, new_p)
    new_prob = remake(prob, p = new_params)
    new_sol = solve(new_prob, Tsit5())
    sum(new_sol)
end

@testset "Issue#2997" begin
    pars = @parameters y0 mh Tγ0 Th0 h ργ0
    vars = @variables x(t)
    @named sys = System([D(x) ~ y0],
        t,
        vars,
        pars;
        defaults = [
            y0 => mh * 3.1 / (2.3 * Th0),
            mh => 123.4,
            Th0 => (4 / 11)^(1 / 3) * Tγ0,
            Tγ0 => (15 / π^2 * ργ0 * (2 * h)^2 / 7)^(1 / 4) / 5
        ])
    sys = mtkcompile(sys)

    function x_at_0(θ)
        prob = ODEProblem(sys, [sys.x => 1.0, sys.ργ0 => θ[1], sys.h => θ[2]], (0.0, 1.0))
        return prob.u0[1]
    end

    @test ForwardDiff.gradient(x_at_0, [0.3, 0.7]) == zeros(2)
end

@parameters a b[1:3] c(t) d::Integer e[1:3] f[1:3, 1:3]::Int g::Vector{AbstractFloat} h::String
@named sys = System(
    Equation[], t, [], [a, b, c, d, e, f, g, h],
    continuous_events = [ModelingToolkit.SymbolicContinuousCallback(
        [a ~ 0] => [c ~ 0], discrete_parameters = c, iv = t)])
sys = complete(sys)

ivs = Dict(c => 3a, b => ones(3), a => 1.0, d => 4, e => [5.0, 6.0, 7.0],
    f => ones(Int, 3, 3), g => [0.1, 0.2, 0.3], h => "foo")

ps = MTKParameters(sys, ivs)

varmap = Dict(a => 1.0f0, b => 3ones(Float32, 3), c => 2.0,
    e => Float32[0.4, 0.5, 0.6], g => ones(Float32, 4))
get_values = getp(sys, [a, b..., c, e...])
get_g = getp(sys, g)
for (_idxs, vals) in [
    # all portions
    (collect(keys(varmap)), collect(values(varmap))),
    # non-arrays
    (keys(varmap), values(varmap)),
    # tunable only
    ([a], [varmap[a]]),
    ([a, b], (varmap[a], varmap[b])),
    ([a, b[2]], (varmap[a], varmap[b][2]))
]
    for idxs in [_idxs, map(i -> parameter_index(sys, i), collect(_idxs))]
        loss = function (p)
            newps = remake_buffer(sys, ps, idxs, p)
            return sum(get_values(newps)) + sum(get_g(newps))
        end

        grad = Zygote.gradient(loss, vals)[1]
        for (val, g) in zip(vals, grad)
            @test eltype(val) == eltype(g)
            if val isa Number
                @test isone(g)
            else
                @test all(isone, g)
            end
        end
    end
end

idxs = (parameter_index(sys, a), parameter_index(sys, b))
vals = (1.0f0, 3ones(Float32, 3))
tangent = rand_tangent(ps)
fwd, back = ChainRulesCore.rrule(remake_buffer, sys, ps, idxs, vals)
@inferred back(tangent)

@testset "Dual type promotion in remake with dummy derivatives" begin # https://github.com/SciML/ModelingToolkit.jl/issues/3336
    # Throw ball straight up into the air
    @variables y(t)
    eqs = [D(D(y)) ~ -9.81]
    initialization_eqs = [y^2 ~ 0] # initialize y = 0 in a way that builds an initialization problem
    @named sys = System(eqs, t; initialization_eqs)
    sys = mtkcompile(sys)

    # Find initial throw velocity that reaches exactly 10 m after 1 s
    dprob0 = ODEProblem(sys, [D(y) => NaN], (0.0, 1.0); guesses = [y => 0.0])
    function f(ics, _)
        dprob = remake(dprob0, u0 = Dict(D(y) => ics[1]))
        dsol = solve(dprob, Tsit5())
        return [dsol[y][end] - 10.0]
    end
    nprob = NonlinearProblem(f, [1.0])
    nsol = solve(nprob, NewtonRaphson())
    @test nsol[1] ≈ 10.0 / 1.0 + 9.81 * 1.0 / 2 # anal free fall solution is y = v0*t - g*t^2/2 -> v0 = y/t + g*t/2
end

@testset "`sys.var` is non-differentiable" begin
    @variables x(t)
    @mtkcompile sys = System(D(x) ~ x, t)
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0))

    grad = Zygote.gradient(prob) do prob
        prob[sys.x]
    end
end

@testset "`p` provided to `solve` is respected" begin
    @mtkmodel Linear begin
        @variables begin
            x(t) = 1.0, [description = "Prey"]
        end
        @parameters begin
            α = 1.5
        end
        @equations begin
            D(x) ~ -α * x
        end
    end

    @mtkcompile linear = Linear()
    problem = ODEProblem(linear, [], (0.0, 1.0))
    solution = solve(problem, Tsit5(), saveat = 0.1)
    rng = StableRNG(42)
    data = (;
        t = solution.t,
        # [[y, x], :]
        measurements = Array(solution)
    )
    data.measurements .+= 0.05 * randn(rng, size(data.measurements))

    p0, repack, _ = SciMLStructures.canonicalize(SciMLStructures.Tunable(), problem.p)

    objective = let repack = repack, problem = problem
        (p, data) -> begin
            pnew = repack(p)
            sol = solve(problem, Tsit5(), p = pnew, saveat = data.t)
            sum(abs2, sol .- data.measurements) / size(data.t, 1)
        end
    end

    # Check 0.0031677344878386607 
    @test_nowarn objective(p0, data)

    fd = ForwardDiff.gradient(Base.Fix2(objective, data), p0)
    zg = Zygote.gradient(Base.Fix2(objective, data), p0)

    @test fd≈zg[1] atol=1e-6
end
