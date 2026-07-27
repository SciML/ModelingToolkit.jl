using ModelingToolkitBase, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D, as_any_buffer
using SymbolicIndexingInterface: getp, setp, remake_buffer
import SciMLStructures

# Two distinct nonnumeric parameter types. Under the opaque nonnumeric buffer they
# must share one `MTKParameters` type while every read/write path keeps seeing the
# concrete payload.

struct FA
    c::Vector{Float64}
end
(f::FA)(τ) = f.c[1] * τ

struct FB
    c::Vector{Float64}
end
(f::FB)(τ) = f.c[1] * τ^2

@parameters (f::FA)(..) d
@variables z(t)
sys = complete(System([D(z) ~ -d * z + f(t)], t; name = :R))
prob = ODEProblem(sys, Dict(z => 0.0, d => 1.0, f => FA([2.0])), (0.0, 1.0))
p = prob.p

@testset "opaque nonnumeric buffer" begin
    @testset "reads" begin
        @test p.nonnumeric isa Tuple
        @test p.nonnumeric[1][1] isa FA
        @test getp(sys, d)(p) == 1.0
        @test getp(sys, f)(p) isa FA
    end

    @testset "setp" begin
        setp(sys, f)(p, FA([7.0]))
        @test getp(sys, f)(p).c == [7.0]
        setp(sys, d)(p, 3.0)
        @test getp(sys, d)(p) == 3.0
    end

    # `remake_buffer` computes promoted buffer element types from the
    # `MTKParameters` type parameters. The nonnumeric element types are no longer
    # in the type, so that buffer is rebuilt dynamically — cover both a numeric
    # and a nonnumeric edit.
    @testset "remake_buffer" begin
        p2 = remake_buffer(sys, p, (d,), (9.0,))
        @test getp(sys, d)(p2) == 9.0
        @test getp(sys, f)(p2).c == [7.0]

        p3 = remake_buffer(sys, p, (f,), (FA([1.5]),))
        @test getp(sys, f)(p3).c == [1.5]
    end

    @testset "copy / equality / as_any_buffer" begin
        pc = copy(p)
        @test pc == p
        @test pc.nonnumeric[1][1].c == p.nonnumeric[1][1].c
        @test as_any_buffer(p).nonnumeric isa Tuple
    end

    @testset "SciMLStructures" begin
        buf, _, _ = SciMLStructures.canonicalize(SciMLStructures.Tunable(), p)
        @test buf isa AbstractVector
        p4 = SciMLStructures.replace(SciMLStructures.Tunable(), p, buf .* 2)
        @test getp(sys, d)(p4) == 2 * getp(sys, d)(p)
        @test getp(sys, f)(p4) isa FA
    end

    @testset "type is uniform across nonnumeric types" begin
        @parameters (g::FB)(..) e
        @variables w(t)
        sys2 = complete(System([D(w) ~ -e * w + g(t)], t; name = :S))
        p2 = ODEProblem(sys2, Dict(w => 0.0, e => 1.0, g => FB([3.0])), (0.0, 1.0)).p
        @test typeof(p2) === typeof(p)
    end
end
