using Test
using ModelingToolkit, Graphs, JumpProcesses, RecursiveArrayTools
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit: value

#################################
#  testing for Jumps / all dgs
#################################
@testset "JumpSystem" begin
    @parameters k1 k2
    @variables S(t) I(t) R(t)
    j₁ = MassActionJump(k1, [0 => 1], [S => 1])
    j₂ = MassActionJump(k1, [S => 1], [S => -1])
    j₃ = MassActionJump(k2, [S => 1, I => 1], [S => -1, I => 1])
    j₄ = MassActionJump(k2, [S => 2, R => 1], [R => -1])
    j₅ = ConstantRateJump(k1 * I, [R ~ R + 1])
    j₆ = VariableRateJump(k1 * k2 / (1 + t) * S, [S ~ S - 1, R ~ R + 1])
    alleqs = [j₁, j₂, j₃, j₄, j₅, j₆]
    @named js = JumpSystem(alleqs, t, [S, I, R], [k1, k2])
    S = value(S)
    I = value(I)
    R = value(R)
    k1 = value(k1)
    k2 = value(k2)

    test_case_1 = (;
        eqs = jumps(js),
        # eq to vars they depend on
        eq_sdeps = [Variable[], [S], [S, I], [S, R], [I], [S]],
        eq_sidepsf = [Int[], [1], [1, 2], [1, 3], [2], [1]],
        eq_sidepsb = [[2, 3, 4, 6], [3, 5], [4]],
        # eq to params they depend on
        eq_pdeps = [[k1], [k1], [k2], [k2], [k1], [k1, k2]],
        eq_pidepsf = [[1], [1], [2], [2], [1], [1, 2]],
        eq_pidepsb = [[1, 2, 5, 6], [3, 4, 6]],
        # var to eqs that modify them
        s_eqdepsf = [[1, 2, 3, 6], [3], [4, 5, 6]],
        s_eqdepsb = [[1], [1], [1, 2], [3], [3], [1, 3]],
        var_eq_ne = 8,
        # eq to eqs that depend on them
        eq_eqdeps = [[2, 3, 4, 6], [2, 3, 4, 6], [2, 3, 4, 5, 6], [4], [4], [2, 3, 4, 6]],
        eq_eq_ne = 6,
        # var to vars that depend on them
        var_vardeps = [[1, 2, 3], [1, 2, 3], [3]],
        var_var_ne = 3
    )
    # testing when ignoring VariableRateJumps
    test_case_2 = (;
        # filter out vrjs in making graphs
        eqs = filter(x -> !(x isa VariableRateJump), jumps(js)),
        # eq to vars they depend on
        eq_sdeps = [Variable[], [S], [S, I], [S, R], [I]],
        eq_sidepsf = [Int[], [1], [1, 2], [1, 3], [2]],
        eq_sidepsb = [[2, 3, 4], [3, 5], [4]],
        # eq to params they depend on
        eq_pdeps = [[k1], [k1], [k2], [k2], [k1]],
        eq_pidepsf = [[1], [1], [2], [2], [1]],
        eq_pidepsb = [[1, 2, 5], [3, 4]],
        # var to eqs that modify them
        s_eqdepsf = [[1, 2, 3], [3], [4, 5]],
        s_eqdepsb = [[1], [1], [1, 2], [3], [3]],
        var_eq_ne = 6,
        # eq to eqs that depend on them
        eq_eqdeps = [[2, 3, 4], [2, 3, 4], [2, 3, 4, 5], [4], [4], [2, 3, 4]],
        eq_eq_ne = 5,
        # var to vars that depend on them
        var_vardeps = [[1, 2, 3], [1, 2, 3], [3]],
        var_var_ne = 3
    )

    @testset "Case $i" for (i, test_case) in enumerate([test_case_1, test_case_2])
        (;         # filter out vrjs in making graphs
            eqs,         # eq to vars they depend on
            eq_sdeps,
            eq_sidepsf,
            eq_sidepsb,         # eq to params they depend on
            eq_pdeps,
            eq_pidepsf,
            eq_pidepsb,         # var to eqs that modify them
            s_eqdepsf,
            s_eqdepsb,
            var_eq_ne,         # eq to eqs that depend on them
            eq_eqdeps,
            eq_eq_ne,         # var to vars that depend on them
            var_vardeps,
            var_var_ne
        ) = test_case
        deps = equation_dependencies(js; eqs)
        @test length(deps) == length(eq_sdeps)
        @test all([issetequal(a, b) for (a, b) in zip(eq_sdeps, deps)])
        # @test all(i -> )
        # @test all(i -> isequal(Set(eq_sdeps[i]), Set(deps[i])), 1:length(alleqs))
        depsbg = asgraph(js; eqs)
        @test depsbg.fadjlist == eq_sidepsf
        @test depsbg.badjlist == eq_sidepsb

        deps = equation_dependencies(js; variables = parameters(js), eqs)
        @test length(deps) == length(eq_pdeps)
        @test all([issetequal(a, b) for (a, b) in zip(eq_pdeps, deps)])
        depsbg2 = asgraph(js; variables = parameters(js), eqs)
        @test depsbg2.fadjlist == eq_pidepsf
        @test depsbg2.badjlist == eq_pidepsb

        bg = BipartiteGraph(var_eq_ne, s_eqdepsf, s_eqdepsb)
        deps2 = variable_dependencies(js; eqs)
        @test isequal(bg, deps2)

        dg = SimpleDiGraph(eq_eq_ne)
        for (eqidx, eqdeps) in enumerate(eq_eqdeps)
            for eqdepidx in eqdeps
                add_edge!(dg, eqidx, eqdepidx)
            end
        end
        dg3 = eqeq_dependencies(depsbg, deps2)
        @test dg == dg3

        dg = SimpleDiGraph(var_var_ne)
        for (vidx, vdeps) in enumerate(var_vardeps)
            for vdepidx in vdeps
                add_edge!(dg, vidx, vdepidx)
            end
        end
        dg4 = varvar_dependencies(depsbg, deps2)
        @test dg == dg4
    end
end

#####################################
#       testing for ODE/SDEs
#####################################

@testset "ODEs, SDEs" begin
    @parameters k1 k2
    @variables S(t) I(t) R(t)
    eqs = [D(S) ~ k1 - k1 * S - k2 * S * I - k1 * k2 / (1 + t) * S
           D(I) ~ k2 * S * I
           D(R) ~ -k2 * S^2 * R / 2 + k1 * I + k1 * k2 * S / (1 + t)]
    @named os = System(eqs, t, [S, I, R], [k1, k2])
    deps = equation_dependencies(os)
    S = value(S)
    I = value(I)
    R = value(R)
    k1 = value(k1)
    k2 = value(k2)
    eq_sdeps = [[S, I], [S, I], [S, I, R]]
    @test all(i -> isequal(Set(eq_sdeps[i]), Set(deps[i])), 1:length(deps))

    noiseeqs = [S, I, R]
    @named sdes = SDESystem(eqs, noiseeqs, t, [S, I, R], [k1, k2])
    deps = equation_dependencies(sdes)
    @test all(i -> isequal(Set(eq_sdeps[i]), Set(deps[i])), 1:length(deps))

    deps = variable_dependencies(os)
    s_eqdeps = [[1], [2], [3]]
    @test deps.fadjlist == s_eqdeps
end

#####################################
#       testing for nonlin sys
#####################################
@testset "Nonlinear" begin
    @variables x y z
    @parameters σ ρ β

    eqs = [0 ~ σ * (y - x),
        0 ~ ρ - y,
        0 ~ y - β * z]
    @named ns = System(eqs, [x, y, z], [σ, ρ, β])
    deps = equation_dependencies(ns)
    eq_sdeps = [[x, y], [y], [y, z]]
    @test all(i -> isequal(Set(deps[i]), Set(value.(eq_sdeps[i]))), 1:length(deps))
end
