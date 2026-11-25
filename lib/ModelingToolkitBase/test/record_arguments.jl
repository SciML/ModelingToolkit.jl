using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using Symbolics
using Test

struct RecordArgumentPair
    x::Float64
    y::Float64
end
struct RecordArgumentOuter
    pair::RecordArgumentPair
    unused::Float64
end
struct RecordArgumentArray
    values::Vector{Float64}
end
struct RecordArgumentEmpty end

@symstruct RecordArgumentPair
@symstruct RecordArgumentOuter
@symstruct RecordArgumentArray
@symstruct RecordArgumentEmpty

record_argument_sum(r::RecordArgumentPair) = r.x + r.y
record_array_argument_sum(r::RecordArgumentArray, x) = sum(r.values) + x
record_empty_argument_sum(::RecordArgumentEmpty, x) = 3.0 + x
@register_symbolic record_argument_sum(r::RecordArgumentPair)
@register_symbolic record_array_argument_sum(r::RecordArgumentArray, x)
@register_symbolic record_empty_argument_sum(r::RecordArgumentEmpty, x)

@testset "Whole record with differential leaves" begin
    @variables s(t)::RecordArgumentPair
    @named sys = System([D(s.x) ~ -s.x, D(s.y) ~ record_argument_sum(s)], t, [s], [])
    sys = mtkcompile(sys)
    @test issetequal(unknowns(sys), Symbolics.unwrap.([s.x, s.y]))
    prob = ODEProblem(sys, [s.x => 1.0, s.y => 0.0], (0.0, 1.0))
    sol = solve(prob, Tsit5())
    @test sol[s.y][end] ≈ sinh(1) rtol = 1.0e-5
end

@testset "Whole record with an eliminated leaf" begin
    @variables s(t)::RecordArgumentPair z(t)
    @named sys = System([D(s.x) ~ -s.x, s.y ~ 2s.x, D(z) ~ record_argument_sum(s)], t, [s, z], [])
    sys = mtkcompile(sys)
    @test issetequal(unknowns(sys), Symbolics.unwrap.([s.x, z]))
    prob = ODEProblem(sys, [s.x => 1.0, z => 0.0], (0.0, 1.0))
    sol = solve(prob, Tsit5())
    @test sol[z][end] ≈ 3(1 - exp(-1)) rtol = 1.0e-5
    @test sol[record_argument_sum(s)][end] ≈ 3 * sol[s.x][end]
end

@testset "Record-valued field argument" begin
    @variables s(t)::RecordArgumentOuter z(t)
    @named sys = System(
        [D(s.pair.x) ~ -s.pair.x, D(s.pair.y) ~ -s.pair.y, D(z) ~ record_argument_sum(s.pair)],
        t, [s, z], []
    )
    prob = ODEProblem(mtkcompile(sys), [s.pair.x => 1.0, s.pair.y => 2.0, z => 0.0], (0.0, 1.0))
    sol = solve(prob, Tsit5())
    @test sol[z][end] ≈ 3(1 - exp(-1)) rtol = 1.0e-5
end

@testset "Concrete record arguments" begin
    @variables z(t)
    @named sys = System(
        [
            D(z) ~ record_array_argument_sum(RecordArgumentArray([1.0, 2.0]), z) +
                record_empty_argument_sum(RecordArgumentEmpty(), z),
        ], t, [z], []
    )
    prob = ODEProblem(mtkcompile(sys), [z => 0.0], (0.0, 1.0))
    @test prob.f(prob.u0, prob.p, 0.0) == [6.0]
end
