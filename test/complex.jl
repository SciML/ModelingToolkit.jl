using ModelingToolkit
using ModelingToolkit: t_nounits as t
using Test

@mtkmodel ComplexModel begin
    @variables begin
        x(t)
        y(t)
        z(t)::Complex
    end
    @equations begin
        z ~ x + im * y
    end
end
@named mixed = ComplexModel()
@test length(equations(mixed)) == 2

@testset "Complex ODEProblem" begin
    using ModelingToolkit: t_nounits as t, D_nounits as D

    vars = @variables x(t) y(t) z(t)
    pars = @parameters a b

    eqs = [
        D(x) ~ y - x,
        D(y) ~ -x * z + b * abs(z),
        D(z) ~ x * y - a
    ]
    @named modlorenz = System(eqs, t)
    sys = mtkcompile(modlorenz)

    ic = ModelingToolkit.get_index_cache(sys)
    @test ic.tunable_buffer_size.type == Number

    u0 = ComplexF64[-4.0, 5.0, 0.0] .+ randn(ComplexF64, 3)
    p = ComplexF64[5.0, 0.1]
    dict = merge(Dict(unknowns(sys) .=> u0), Dict(parameters(sys) .=> p))
    prob = ODEProblem(sys, dict, (0.0, 1.0))

    using OrdinaryDiffEq
    sol = solve(prob, Tsit5(), saveat = 0.1)

    @test sol.u[1] isa Vector{ComplexF64}
end
