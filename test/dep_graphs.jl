using Test
using ModelingToolkit, Graphs, JumpProcesses, RecursiveArrayTools
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit: value

#################################
#  testing for Jumps / all dgs
#################################
@parameters k1 k2
@variables S(t) I(t) R(t)
j₁ = MassActionJump(k1, [0 => 1], [S => 1])
j₂ = MassActionJump(k1, [S => 1], [S => -1])
j₃ = MassActionJump(k2, [S => 1, I => 1], [S => -1, I => 1])
j₄ = MassActionJump(k2, [S => 2, R => 1], [R => -1])
j₅ = ConstantRateJump(k1 * I, [R ~ R + 1])
j₆ = VariableRateJump(k1 * k2 / (1 + t) * S, [S ~ S - 1, R ~ R + 1])
eqs = [j₁, j₂, j₃, j₄, j₅, j₆]
@named js = JumpSystem(eqs, t, [S, I, R], [k1, k2])
S = value(S)
I = value(I)
R = value(R)
k1 = value(k1)
k2 = value(k2)
# eq to vars they depend on
eq_sdeps = [Variable[], [S], [S, I], [S, R], [I], [S]]
eq_sidepsf = [Int[], [1], [1, 2], [1, 3], [2], [1]]
eq_sidepsb = [[2, 3, 4, 6], [3, 5], [4]]
deps = equation_dependencies(js)
@test all(i -> isequal(Set(eq_sdeps[i]), Set(deps[i])), 1:length(eqs))
depsbg = asgraph(js)
@test depsbg.fadjlist == eq_sidepsf
@test depsbg.badjlist == eq_sidepsb

# eq to params they depend on
eq_pdeps = [[k1], [k1], [k2], [k2], [k1], [k1, k2]]
eq_pidepsf = [[1], [1], [2], [2], [1], [1, 2]]
eq_pidepsb = [[1, 2, 5, 6], [3, 4, 6]]
deps = equation_dependencies(js, variables = parameters(js))
@test all(i -> isequal(Set(eq_pdeps[i]), Set(deps[i])), 1:length(eqs))
depsbg2 = asgraph(js, variables = parameters(js))
@test depsbg2.fadjlist == eq_pidepsf
@test depsbg2.badjlist == eq_pidepsb

# var to eqs that modify them
s_eqdepsf = [[1, 2, 3, 6], [3], [4, 5, 6]]
s_eqdepsb = [[1], [1], [1, 2], [3], [3], [1, 3]]
ne = 8
bg = BipartiteGraph(ne, s_eqdepsf, s_eqdepsb)
deps2 = variable_dependencies(js)
@test isequal(bg, deps2)

# eq to eqs that depend on them
eq_eqdeps = [[2, 3, 4, 6], [2, 3, 4, 6], [2, 3, 4, 5, 6], [4], [4], [2, 3, 4, 6]]
dg = SimpleDiGraph(6)
for (eqidx, eqdeps) in enumerate(eq_eqdeps)
    for eqdepidx in eqdeps
        add_edge!(dg, eqidx, eqdepidx)
    end
end
dg3 = eqeq_dependencies(depsbg, deps2)
@test dg == dg3

# var to vars that depend on them
var_vardeps = [[1, 2, 3], [1, 2, 3], [3]]
ne = 7
dg = SimpleDiGraph(3)
for (vidx, vdeps) in enumerate(var_vardeps)
    for vdepidx in vdeps
        add_edge!(dg, vidx, vdepidx)
    end
end
dg4 = varvar_dependencies(depsbg, deps2)
@test dg == dg4

# testing when ignoring VariableRateJumps
let
    @parameters k1 k2
    @variables S(t) I(t) R(t)
    j₁ = MassActionJump(k1, [0 => 1], [S => 1])
    j₂ = MassActionJump(k1, [S => 1], [S => -1])
    j₃ = MassActionJump(k2, [S => 1, I => 1], [S => -1, I => 1])
    j₄ = MassActionJump(k2, [S => 2, R => 1], [R => -1])
    j₅ = ConstantRateJump(k1 * I, [R ~ R + 1])
    j₆ = VariableRateJump(k1 * k2 / (1 + t) * S, [S ~ S - 1, R ~ R + 1])
    eqs = [j₁, j₂, j₃, j₄, j₅, j₆]
    @named js = JumpSystem(eqs, t, [S, I, R], [k1, k2])
    S = value(S)
    I = value(I)
    R = value(R)
    k1 = value(k1)
    k2 = value(k2)
    # eq to vars they depend on
    eq_sdeps = [Variable[], [S], [S, I], [S, R], [I]]
    eq_sidepsf = [Int[], [1], [1, 2], [1, 3], [2]]
    eq_sidepsb = [[2, 3, 4], [3, 5], [4]]

    # filter out vrjs in making graphs
    eqs = ArrayPartition(equations(js).x[1], equations(js).x[2])
    deps = equation_dependencies(js; eqs)
    @test length(deps) == length(eq_sdeps)
    @test all(i -> isequal(Set(eq_sdeps[i]), Set(deps[i])), 1:length(eqs))
    depsbg = asgraph(js; eqs)
    @test depsbg.fadjlist == eq_sidepsf
    @test depsbg.badjlist == eq_sidepsb

    # eq to params they depend on
    eq_pdeps = [[k1], [k1], [k2], [k2], [k1]]
    eq_pidepsf = [[1], [1], [2], [2], [1]]
    eq_pidepsb = [[1, 2, 5], [3, 4]]
    deps = equation_dependencies(js; variables = parameters(js), eqs)
    @test length(deps) == length(eq_pdeps)
    @test all(i -> isequal(Set(eq_pdeps[i]), Set(deps[i])), 1:length(eqs))
    depsbg2 = asgraph(js; variables = parameters(js), eqs)
    @test depsbg2.fadjlist == eq_pidepsf
    @test depsbg2.badjlist == eq_pidepsb

    # var to eqs that modify them
    s_eqdepsf = [[1, 2, 3], [3], [4, 5]]
    s_eqdepsb = [[1], [1], [1, 2], [3], [3]]
    ne = 6
    bg = BipartiteGraph(ne, s_eqdepsf, s_eqdepsb)
    deps2 = variable_dependencies(js; eqs)
    @test isequal(bg, deps2)

    # eq to eqs that depend on them
    eq_eqdeps = [[2, 3, 4], [2, 3, 4], [2, 3, 4, 5], [4], [4], [2, 3, 4]]
    dg = SimpleDiGraph(5)
    for (eqidx, eqdeps) in enumerate(eq_eqdeps)
        for eqdepidx in eqdeps
            add_edge!(dg, eqidx, eqdepidx)
        end
    end
    dg3 = eqeq_dependencies(depsbg, deps2)
    @test dg == dg3

    # var to vars that depend on them
    var_vardeps = [[1, 2, 3], [1, 2, 3], [3]]
    ne = 7
    dg = SimpleDiGraph(3)
    for (vidx, vdeps) in enumerate(var_vardeps)
        for vdepidx in vdeps
            add_edge!(dg, vidx, vdepidx)
        end
    end
    dg4 = varvar_dependencies(depsbg, deps2)
    @test dg == dg4
end

#####################################
#       testing for ODE/SDEs
#####################################
@parameters k1 k2
@variables S(t) I(t) R(t)
eqs = [D(S) ~ k1 - k1 * S - k2 * S * I - k1 * k2 / (1 + t) * S,
    D(I) ~ k2 * S * I,
    D(R) ~ -k2 * S^2 * R / 2 + k1 * I + k1 * k2 * S / (1 + t)]
noiseeqs = [S, I, R]
@named os = ODESystem(eqs, t, [S, I, R], [k1, k2])
deps = equation_dependencies(os)
S = value(S);
I = value(I);
R = value(R);
k1 = value(k1);
k2 = value(k2);
eq_sdeps = [[S, I], [S, I], [S, I, R]]
@test all(i -> isequal(Set(eq_sdeps[i]), Set(deps[i])), 1:length(deps))

@parameters k1 k2
@variables S(t) I(t) R(t)
@named sdes = SDESystem(eqs, noiseeqs, t, [S, I, R], [k1, k2])
deps = equation_dependencies(sdes)
@test all(i -> isequal(Set(eq_sdeps[i]), Set(deps[i])), 1:length(deps))

deps = variable_dependencies(os)
s_eqdeps = [[1], [2], [3]]
@test deps.fadjlist == s_eqdeps

#####################################
#       testing for nonlin sys
#####################################
@variables x y z
@parameters σ ρ β

eqs = [0 ~ σ * (y - x),
    0 ~ ρ - y,
    0 ~ y - β * z]
@named ns = NonlinearSystem(eqs, [x, y, z], [σ, ρ, β])
deps = equation_dependencies(ns)
eq_sdeps = [[x, y], [y], [y, z]]
@test all(i -> isequal(Set(deps[i]), Set(value.(eq_sdeps[i]))), 1:length(deps))
