using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Zygote
using SymbolicIndexingInterface
using SciMLStructures
using OrdinaryDiffEq
using NonlinearSolve
using SciMLSensitivity
using ForwardDiff
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
@mtkbuild sys = ODESystem(eqs, t)
prob = ODEProblem(sys, u0, tspan, ps)
sol = solve(prob, Tsit5())

mtkparams = parameter_values(prob)
new_p = rand(10)
gs = gradient(new_p) do new_p
    new_params = SciMLStructures.replace(SciMLStructures.Tunable(), mtkparams, new_p)
    new_prob = remake(prob, p = new_params)
    new_sol = solve(new_prob, Tsit5())
    sum(new_sol)
end

@testset "Issue#2997" begin
    pars = @parameters y0 mh Tγ0 Th0 h ργ0
    vars = @variables x(t)
    @named sys = ODESystem([D(x) ~ y0],
        t,
        vars,
        pars;
        defaults = [
            y0 => mh * 3.1 / (2.3 * Th0),
            mh => 123.4,
            Th0 => (4 / 11)^(1 / 3) * Tγ0,
            Tγ0 => (15 / π^2 * ργ0 * (2 * h)^2 / 7)^(1 / 4) / 5
        ])
    sys = structural_simplify(sys)

    function x_at_0(θ)
        prob = ODEProblem(sys, [sys.x => 1.0], (0.0, 1.0), [sys.ργ0 => θ[1], sys.h => θ[2]])
        return prob.u0[1]
    end

    @test ForwardDiff.gradient(x_at_0, [0.3, 0.7]) == zeros(2)
end

@parameters a b[1:3] c(t) d::Integer e[1:3] f[1:3, 1:3]::Int g::Vector{AbstractFloat} h::String
@named sys = ODESystem(
    Equation[], t, [], [a, b, c, d, e, f, g, h],
    continuous_events = [[a ~ 0] => [c ~ 0]])
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
    @named sys = ODESystem(eqs, t; initialization_eqs)
    sys = structural_simplify(sys)

    # Find initial throw velocity that reaches exactly 10 m after 1 s
    dprob0 = ODEProblem(sys, [D(y) => NaN], (0.0, 1.0), []; guesses = [y => 0.0])
    function f(ics, _)
        dprob = remake(dprob0, u0 = Dict(D(y) => ics[1]))
        dsol = solve(dprob, Tsit5())
        return [dsol[y][end] - 10.0]
    end
    nprob = NonlinearProblem(f, [1.0])
    nsol = solve(nprob, NewtonRaphson())
    @test nsol[1] ≈ 10.0 / 1.0 + 9.81 * 1.0 / 2 # anal free fall solution is y = v0*t - g*t^2/2 -> v0 = y/t + g*t/2
end
