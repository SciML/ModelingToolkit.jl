using ModelingToolkit, Test, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: AnalysisPoint, sampletime, is_continuous_partition, is_linearized,
                       input_group, output_group
using SciMLBase: PeriodicClock, ContinuousClock
using OrdinaryDiffEqRosenbrock: Rodas5

"""
The plant of the sampled-data test cases, `1 / (s + 1)`, discretized with a zero-order hold
over `dt`. This is the discretization the assembled closed loop must agree with, since `Hold`
is a zero-order hold and `Sample` is ideal sampling.
"""
function zoh_first_order(dt)
    ad = exp(-dt)
    return ad, 1 - ad
end

@testset "A model without clocks gives one partition equal to `linearize`" begin
    @variables x(t) u(t) y(t)
    @named sys = System([D(x) ~ -2x + 3u, y ~ x], t)
    op = Dict(x => 0.0, u => 0.0)

    partitions = linearize_clocked(sys, [u], [y]; op)
    @test length(partitions) == 1
    part = only(partitions)
    @test is_continuous_partition(part)
    @test sampletime(part) === nothing
    @test is_linearized(part)

    matrices, _, _ = linearize(sys, [u], [y]; op)
    @test part.A == matrices.A
    @test part.B == matrices.B
    @test part.C == matrices.C
    @test part.D == matrices.D

    @test isequal(only(input_group(part, 0)), ModelingToolkit.unwrap(u))
    @test isequal(only(output_group(part, 0)), ModelingToolkit.unwrap(y))
    @test part.input_groups == [0 => 1:1]
    @test part.output_groups == [0 => 1:1]
end

@testset "A sampled-data loop splits into a plant and a clocked controller" begin
    dt = 0.1
    k = ShiftIndex(Clock(dt))
    @variables x(t) y(t) u(t) yd(t) ud(t) r(t)
    @parameters kp
    kpval = 2.0

    # A clocked PI: ud(k) = ud(k-1) + kp * (r(k) - y(k))
    @named sys = System(
        [
            yd ~ Sample(dt)(y),
            ud ~ ud(k - 1) + kp * (r - yd),
            u ~ Hold(ud),
            D(x) ~ -x + u,
            y ~ x,
        ], t)

    op = Dict(x => 0.0, y => 0.0, u => 0.0, ud => 0.0, ud(k - 1) => 0.0,
        r => 0.0, kp => kpval)
    partitions = linearize_clocked(sys, [r], [y]; op)
    @test length(partitions) == 2
    ctrl, plant = partitions

    # The continuous partition is last, and carries the clock it runs on.
    @test is_continuous_partition(plant)
    @test !is_continuous_partition(ctrl)
    @test ctrl.clock == PeriodicClock(dt, 0.0)
    @test sampletime(ctrl) == dt
    @test sampletime(plant) === nothing

    # The signals crossing the boundary are named by the model's own variables, the same on
    # both sides, and are grouped by the partition they come from or go to.
    @test isequal(only(input_group(ctrl, 0)), ModelingToolkit.unwrap(r))
    @test isequal(only(input_group(ctrl, 2)), ModelingToolkit.unwrap(y))
    @test isequal(only(output_group(ctrl, 2)), ModelingToolkit.unwrap(ud))
    @test isequal(only(input_group(plant, 1)), ModelingToolkit.unwrap(ud))
    @test isequal(only(output_group(plant, 0)), ModelingToolkit.unwrap(y))
    @test isequal(only(output_group(plant, 1)), ModelingToolkit.unwrap(y))
    @test isempty(input_group(plant, 0))
    @test isempty(output_group(ctrl, 0))

    # The controller is the difference equation it was written as.
    @test ctrl.A ≈ [1.0;;]
    @test ctrl.B ≈ [kpval -kpval]
    @test ctrl.C ≈ [1.0;;]
    @test ctrl.D ≈ [kpval -kpval]

    # The plant is returned in continuous time.
    @test plant.A ≈ [-1.0;;]
    @test plant.B ≈ [1.0;;]
    @test plant.C ≈ [1.0; 1.0;;]
    @test plant.D ≈ zeros(2, 1)

    # Assembling the loop: the plant held over one sample interval, closed with the
    # controller. Each partition supplies one block of the closed-loop dynamics.
    ad, bd = zoh_first_order(dt)
    reference = [ad-bd*kpval bd; -kpval 1]
    assembled = [ad+bd*ctrl.D[1, 2] bd*ctrl.C[]; ctrl.B[1, 2] ctrl.A[]]
    @test assembled ≈ reference
end

@testset "A clocked partition without a state variable is a static map" begin
    dt = 0.1
    @variables x(t) y(t) u(t) yd(t) ud(t) r(t)
    @parameters kp
    @named sys = System(
        [
            yd ~ Sample(dt)(y),
            ud ~ kp * (r - yd),
            r ~ 1.0,
            u ~ Hold(ud),
            D(x) ~ -x + u,
            y ~ x,
        ], t)

    op = Dict(x => 0.0, y => 0.0, u => 0.0, ud => 0.0, kp => 2.0)
    ctrl, plant = linearize_clocked(sys, ModelingToolkit.SymbolicT[], [y]; op)
    @test size(ctrl.A) == (0, 0)
    @test size(ctrl.B) == (0, 1)
    @test size(ctrl.C) == (1, 0)
    @test ctrl.D ≈ [-2.0;;]
    @test sampletime(ctrl) == dt
    @test plant.A ≈ [-1.0;;]
end

@testset "Two rates give one partition each, plus the continuous one" begin
    dt1, dt2 = 0.1, 0.2
    k2 = ShiftIndex(Clock(dt2))
    @variables x(t) y(t) u(t) yd(t) ud(t) rd(t) rdc(t) ys(t) rs(t) ref(t)
    @parameters kp ki
    kpval, kival = 2.0, 0.5

    @named sys = System(
        [
            D(x) ~ -x + u,
            y ~ x,
            u ~ Hold(ud),
            yd ~ Sample(dt1)(y),
            rd ~ Sample(dt1)(rdc),
            ud ~ kp * (rd - yd),
            ys ~ Sample(dt2)(y),
            rs ~ rs(k2 - 1) + ki * (ref - ys),
            rdc ~ Hold(rs),
        ], t)

    op = Dict(x => 0.0, y => 0.0, u => 0.0, ud => 0.0, rdc => 0.0, rs => 0.0,
        rs(k2 - 1) => 0.0, ref => 0.0, kp => kpval, ki => kival)
    partitions = linearize_clocked(sys, [ref], [y]; op)
    @test length(partitions) == 3
    slow, fast, plant = partitions

    @test sampletime(slow) == dt2
    @test sampletime(fast) == dt1
    @test is_continuous_partition(plant)

    @test slow.A ≈ [1.0;;]
    @test slow.B ≈ [kival -kival]
    @test slow.D ≈ [kival -kival]
    @test size(fast.A) == (0, 0)
    @test fast.D ≈ [kpval -kpval]

    # Each clocked partition reads the plant output and writes one plant input.
    @test isequal(only(input_group(slow, 3)), ModelingToolkit.unwrap(y))
    @test isequal(only(output_group(slow, 3)), ModelingToolkit.unwrap(rs))
    @test issetequal(input_group(fast, 3),
        ModelingToolkit.unwrap.([rdc, y]))
    @test isequal(only(output_group(fast, 3)), ModelingToolkit.unwrap(ud))
    @test isequal(only(input_group(plant, 1)), ModelingToolkit.unwrap(rs))
    @test isequal(only(input_group(plant, 2)), ModelingToolkit.unwrap(ud))
    @test issetequal(output_group(plant, 2), ModelingToolkit.unwrap.([rdc]))

    # The reference passed to the fast loop is held, so it reaches the fast loop without
    # going through the plant dynamics.
    rdc_row = findfirst(isequal(ModelingToolkit.unwrap(rdc)), plant.outputs)
    rs_col = findfirst(isequal(ModelingToolkit.unwrap(rs)), plant.inputs)
    @test plant.D[rdc_row, rs_col] ≈ 1.0
end

@testset "Analysis points name the inputs and outputs" begin
    dt = 0.1
    @variables x(t) y(t) [output = true] yd(t) ud(t) uh(t) [output = true]
    @variables u(t) [input = true]
    @parameters kp
    kpval = 2.0

    @named sys = System(
        [
            yd ~ Sample(dt)(y),
            ud ~ -kp * yd,
            uh ~ Hold(ud),
            connect(uh, AnalysisPoint(:plant_input), u),
            D(x) ~ -x + u,
            y ~ x,
        ], t)

    op = Dict(x => 0.0, y => 0.0, u => 0.0, uh => 0.0, ud => 0.0, kp => kpval)
    partitions = linearize_clocked(sys, :plant_input, [y]; op)
    @test length(partitions) == 2
    ctrl, plant = partitions

    # The perturbation of the analysis point is an input of the continuous partition.
    @test length(input_group(plant, 0)) == 1
    @test isempty(input_group(ctrl, 0))
    @test plant.A ≈ [-1.0;;]
    @test plant.B[:, 1] ≈ [1.0]
    @test ctrl.D ≈ [-kpval;;]
end

@testset "An unsupported clock and an event are reported" begin
    dt = 0.1
    @variables x(t) y(t) u(t) yd(t) ud(t)
    @parameters kp

    reset = ModelingToolkit.ImperativeAffect(modified = (; x)) do m, o, _, i
        return (; x = 0.0)
    end
    @named evsys = System(
        [
            yd ~ Sample(dt)(y),
            ud ~ -kp * yd,
            u ~ Hold(ud),
            D(x) ~ -x + u,
            y ~ x,
        ], t; continuous_events = [[x ~ 0.5] => reset])
    op = Dict(x => 0.0, y => 0.0, u => 0.0, ud => 0.0, kp => 1.0)
    @test_warn ["continuous event", "not included"] linearize_clocked(
        evsys, ModelingToolkit.SymbolicT[], [y]; op)

    @named stepsys = System(
        [
            yd ~ Sample(SolverStepClock())(y),
            ud ~ -kp * yd,
            u ~ Hold(ud),
            D(x) ~ -x + u,
            y ~ x,
        ], t)
    partitions = @test_warn ["not a periodic clock", "not included"] linearize_clocked(
        stepsys, ModelingToolkit.SymbolicT[], [y]; op)
    stepped = partitions[findfirst(!is_continuous_partition, partitions)]
    @test !is_linearized(stepped)
    @test stepped.A === nothing
    # The signal lists survive, so a model of that partition can be supplied by other means.
    @test isequal(only(input_group(stepped, 2)), ModelingToolkit.unwrap(y))
    @test isequal(only(output_group(stepped, 2)), ModelingToolkit.unwrap(ud))
end

@testset "An operating-point entry that names nothing is reported" begin
    dt = 0.1
    @variables x(t) y(t) u(t) yd(t) ud(t) nothere(t)
    @parameters kp
    @named sys = System(
        [
            yd ~ Sample(dt)(y),
            ud ~ -kp * yd,
            u ~ Hold(ud),
            D(x) ~ -x + u,
            y ~ x,
        ], t)
    op = Dict(x => 0.0, y => 0.0, u => 0.0, ud => 0.0, kp => 1.0, nothere => 3.0)
    @test_warn ["nothere", "ignored"] linearize_clocked(
        sys, ModelingToolkit.SymbolicT[], [y]; op)
end

@testset "An operating point taken from a solution" begin
    dt = 0.1
    k = ShiftIndex(Clock(dt))
    @variables x(t) y(t) u(t) yd(t) ud(t)
    @parameters kp
    @named sys = System(
        [
            yd ~ Sample(dt)(y),
            ud ~ ud(k - 1) - kp * yd,
            u ~ Hold(ud),
            D(x) ~ -x + u,
            y ~ x,
        ], t)

    # The plant alone can be simulated; the solution then supplies the continuous states of
    # the operating point, and the rest is given explicitly.
    @variables xp(t) up(t)
    @named plant = System([D(xp) ~ -xp + up], t)
    psys = mtkcompile(plant; inputs = [up])
    prob = ODEProblem(psys, [xp => 1.0, up => 0.0], (0.0, 1.0))
    sol = solve(prob, Rodas5())

    extra = Dict(y => sol(0.5; idxs = xp), u => 0.0, ud => 0.0,
        ud(k - 1) => 0.0, kp => 2.0)
    partitions = linearize_clocked(
        sys, ModelingToolkit.SymbolicT[], [y];
        op = Dict(x => sol(0.5; idxs = xp), extra...))
    ctrl, cont = partitions
    @test cont.extras.x ≈ [sol(0.5; idxs = xp)]
    @test ctrl.A ≈ [1.0;;]
end

"A function of a whole vector, the way a controller calls a solver."
vector_gain(v) = sum(v)
@register_symbolic vector_gain(v::AbstractVector)

@testset "An array-valued clocked variable used at a previous tick" begin
    dt = 0.1
    k = ShiftIndex(Clock(dt))
    @variables x1(t) x2(t) y1(t) y2(t) u(t) ud(t)
    @variables (xd(t))[1:2]
    @parameters kp
    kpval = 3.0

    # The elements of `xd` are aliased to the sampled measurements and the control law is a
    # function of the whole vector one tick back, as a model-predictive controller calling its
    # solver with the whole state estimate is.
    loop_equations(spelling) = [
        xd[1] ~ Sample(dt)(y1),
        xd[2] ~ Sample(dt)(y2),
        ud ~ kp * spelling,
        u ~ Hold(ud),
        D(x1) ~ -x1 + u,
        D(x2) ~ -2x2 + u,
        y1 ~ x1,
        y2 ~ x2,
    ]
    @named whole = System(loop_equations(vector_gain(xd(k - 1))), t)
    @named elementwise = System(loop_equations(xd(k - 1)[1] + xd(k - 1)[2]), t)

    op = Dict(x1 => 0.0, x2 => 0.0, y1 => 0.0, y2 => 0.0, u => 0.0, ud => 0.0,
        xd(k - 1)[1] => 0.0, xd(k - 1)[2] => 0.0, kp => kpval)
    ctrl, plant = linearize_clocked(whole, ModelingToolkit.SymbolicT[], [y1]; op)

    @test sampletime(ctrl) == dt
    # The state variables are the two elements of `xd` at the previous tick, each updated by
    # the measurement it is aliased to, and the control law sums them. The order the elements
    # are given is an implementation detail, so the assertion is on the input-output map.
    @test ctrl.A ≈ zeros(2, 2)
    @test ctrl.D ≈ zeros(1, 2)
    @test ctrl.C * ctrl.B ≈ [kpval kpval]
    @test sort(real(eigvals(plant.A))) ≈ [-2.0, -1.0]

    # The two spellings of the same control law must give the same matrices.
    ctrl2, plant2 = linearize_clocked(elementwise, ModelingToolkit.SymbolicT[], [y1]; op)
    @test ctrl2.A ≈ ctrl.A
    @test ctrl2.C * ctrl2.B ≈ ctrl.C * ctrl.B
    @test ctrl2.D ≈ ctrl.D
    @test plant2.A ≈ plant.A
end
