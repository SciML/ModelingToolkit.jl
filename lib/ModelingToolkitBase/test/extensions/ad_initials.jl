# https://github.com/SciML/SciMLSensitivity.jl/issues/1582
# Lives under the parent ModelingToolkit Extensions group (`@mtktestset`): the
# mass-matrix RHS needs full MTK `mtkcompile`, which ModelingToolkitBase refuses.
# Imports mirror the context `ad.jl` provided when this testset lived there.
#
# Buffer is taken from `p` inside the Zygote closure (as in #4905's verification
# table: `replace(Initials, p, p.initials)`), not closed over from `p0` outside.
# That is what #4905 meant to test — the Initials cotangent survives an
# `MTKParameters` rebuild (see that PR's body / commit that added this testset).
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Zygote
using SymbolicIndexingInterface
using SciMLStructures
using OrdinaryDiffEqTsit5
using SciMLSensitivity
using ForwardDiff
using ChainRulesCore
using Test

@testset "`Initials` cotangent survives an `MTKParameters` rebuild" begin
    @parameters m = 1.5 d = 9.0
    @variables s(t) = 1.0 v(t) = 1.0
    sys = mtkcompile(System([D(s) ~ v, m * D(v) ~ 1 - d * v], t; name = :model))
    prob = ODEProblem(sys, [], (0.0, 1.0))
    p0 = parameter_values(prob)

    # `setproperties`/`@set!` on an `MTKParameters` goes through its positional
    # constructor. `promote_u0_p` does exactly this inside `remake`, so a rebuild
    # that drops the non-tunable cotangents silently zeroes every `Initial(x)`
    # gradient downstream.
    for portion in (SciMLStructures.Tunable(), SciMLStructures.Initials())
        g = Zygote.gradient(p0) do p
            buf = SciMLStructures.canonicalize(portion, p)[1]
            return sum(SciMLStructures.replace(portion, p, buf).initials)
        end[1]
        @test g.initials == ones(length(p0.initials))
    end

    # End-to-end: `u0` is a derived quantity here, so `Initial(v)` must reach it.
    set_v0 = setp_oop(prob, [Initial(sys.v)])
    f = x -> sum(state_values(remake(prob; p = set_v0(prob, x))))
    @test Zygote.gradient(f, [1.1])[1] ≈ ForwardDiff.gradient(f, [1.1])
end
