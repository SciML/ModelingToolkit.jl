using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, MTKParameters
using SciMLBase
using SciMLStructures: SciMLStructures, canonicalize, Tunable, Discrete, Constants
using JET
using Test

# scalar parameters only
function level1()
    @parameters p1 = 0.5 [tunable = true] p2 = 1 [tunable = true] p3 = 3 [tunable = false] p4 = 3 [tunable = true] y0
    @variables x(t) = 2 y(t) = y0
    D = Differential(t)

    eqs = [
        D(x) ~ p1 * x - p2 * x * y
        D(y) ~ -p3 * y + p4 * x * y
    ]

    sys = mtkcompile(
        complete(
            System(
                eqs, t, name = :sys, bindings = [y0 => 2p4]
            )
        )
    )
    return prob = ODEProblem{true, SciMLBase.FullSpecialize}(sys, [], (0.0, 3.0))
end

# scalar and vector parameters
function level2()
    @parameters p1 = 0.5 [tunable = true] (p23[1:2] = [1, 3.0]) [tunable = true] p4 = 3 [tunable = false] y0
    @variables x(t) = 2 y(t) = y0
    D = Differential(t)

    eqs = [
        D(x) ~ p1 * x - p23[1] * x * y
        D(y) ~ -p23[2] * y + p4 * x * y
    ]

    sys = mtkcompile(
        complete(
            System(
                eqs, t, name = :sys, bindings = [y0 => 2p4]
            )
        )
    )
    return prob = ODEProblem{true, SciMLBase.FullSpecialize}(sys, [], (0.0, 3.0))
end

# scalar and vector parameters with different scalar types
function level3()
    @parameters p1 = 0.5 [tunable = true] (p23[1:2] = [1, 3.0]) [tunable = true] p4::Int = 3 [tunable = true] y0::Int
    @variables x(t) = 2 y(t) = y0
    D = Differential(t)

    eqs = [
        D(x) ~ p1 * x - p23[1] * x * y
        D(y) ~ -p23[2] * y + p4 * x * y
    ]

    sys = mtkcompile(
        complete(
            System(
                eqs, t, name = :sys, bindings = [y0 => 2p4]
            )
        )
    )
    return prob = ODEProblem{true, SciMLBase.FullSpecialize}(sys, [], (0.0, 3.0))
end

@testset "level$i" for (i, prob) in enumerate([level1(), level2(), level3()])
    ps = prob.p
    @testset "Type stability of $portion" for portion in [
            Tunable(), Discrete(), Constants(),
        ]
        @test_call canonicalize(portion, ps)

        # broken because the size of a vector of vectors can't be determined at compile time
        @test_opt target_modules = (ModelingToolkitBase,) canonicalize(
            portion, ps
        )

        buffer, repack, alias = canonicalize(portion, ps)

        # broken because dependent update functions break inference
        @test_call target_modules = (ModelingToolkitBase,) SciMLStructures.replace(
            portion, ps, ones(length(buffer))
        )
        @test_opt target_modules = (ModelingToolkitBase,) SciMLStructures.replace(
            portion, ps, ones(length(buffer))
        )

        @test_call target_modules = (ModelingToolkitBase,) SciMLStructures.replace!(
            portion, ps, ones(length(buffer))
        )
        @test_opt target_modules = (ModelingToolkitBase,) SciMLStructures.replace!(
            portion, ps, ones(length(buffer))
        )
    end
end
