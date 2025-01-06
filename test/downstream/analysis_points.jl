using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, ControlSystemsBase
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: connect, t_nounits as t, D_nounits as D
import ControlSystemsBase as CS

@testset "Complicated model" begin
    # Parameters
    m1 = 1
    m2 = 1
    k = 1000 # Spring stiffness
    c = 10   # Damping coefficient
    @named inertia1 = Inertia(; J = m1)
    @named inertia2 = Inertia(; J = m2)
    @named spring = Spring(; c = k)
    @named damper = Damper(; d = c)
    @named torque = Torque()

    function SystemModel(u = nothing; name = :model)
        eqs = [connect(torque.flange, inertia1.flange_a)
               connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
               connect(inertia2.flange_a, spring.flange_b, damper.flange_b)]
        if u !== nothing
            push!(eqs, connect(torque.tau, u.output))
            return ODESystem(eqs, t;
                systems = [
                    torque,
                    inertia1,
                    inertia2,
                    spring,
                    damper,
                    u
                ],
                name)
        end
        ODESystem(eqs, t; systems = [torque, inertia1, inertia2, spring, damper], name)
    end

    @named r = Step(start_time = 0)
    model = SystemModel()
    @named pid = PID(k = 100, Ti = 0.5, Td = 1)
    @named filt = SecondOrder(d = 0.9, w = 10)
    @named sensor = AngleSensor()
    @named er = Add(k2 = -1)

    connections = [connect(r.output, :r, filt.input)
                   connect(filt.output, er.input1)
                   connect(pid.ctr_output, :u, model.torque.tau)
                   connect(model.inertia2.flange_b, sensor.flange)
                   connect(sensor.phi, :y, er.input2)
                   connect(er.output, :e, pid.err_input)]

    closed_loop = ODESystem(connections, t, systems = [model, pid, filt, sensor, r, er],
        name = :closed_loop, defaults = [
            model.inertia1.phi => 0.0,
            model.inertia2.phi => 0.0,
            model.inertia1.w => 0.0,
            model.inertia2.w => 0.0,
            filt.x => 0.0,
            filt.xd => 0.0
        ])

    sys = structural_simplify(closed_loop)
    prob = ODEProblem(sys, unknowns(sys) .=> 0.0, (0.0, 4.0))
    sol = solve(prob, Rodas5P(), reltol = 1e-6, abstol = 1e-9)

    matrices, ssys = linearize(closed_loop, :r, :y)
    lsys = ss(matrices...) |> sminreal
    @test lsys.nx == 8

    stepres = ControlSystemsBase.step(c2d(lsys, 0.001), 4)
    @test Array(stepres.y[:])≈Array(sol(0:0.001:4, idxs = model.inertia2.phi)) rtol=1e-4

    matrices, ssys = get_sensitivity(closed_loop, :y)
    So = ss(matrices...)

    matrices, ssys = get_sensitivity(closed_loop, :u)
    Si = ss(matrices...)

    @test tf(So) ≈ tf(Si)
end

@testset "Analysis points with subsystems" begin
    @named P = FirstOrder(k = 1, T = 1)
    @named C = Gain(; k = 1)
    @named add = Blocks.Add(k2 = -1)

    eqs = [connect(P.output, :plant_output, add.input2)
           connect(add.output, C.input)
           connect(C.output, :plant_input, P.input)]

    # eqs = [connect(P.output, add.input2)
    #        connect(add.output, C.input)
    #        connect(C.output, P.input)]

    sys_inner = ODESystem(eqs, t, systems = [P, C, add], name = :inner)

    @named r = Constant(k = 1)
    @named F = FirstOrder(k = 1, T = 3)

    eqs = [connect(r.output, F.input)
           connect(F.output, sys_inner.add.input1)]
    sys_outer = ODESystem(eqs, t, systems = [F, sys_inner, r], name = :outer)

    # test first that the structural_simplify works correctly
    ssys = structural_simplify(sys_outer)
    prob = ODEProblem(ssys, Pair[], (0, 10))
    @test_nowarn solve(prob, Rodas5())

    matrices, _ = get_sensitivity(sys_outer, sys_outer.inner.plant_input)
    lsys = sminreal(ss(matrices...))
    @test lsys.A[] == -2
    @test lsys.B[] * lsys.C[] == -1 # either one negative
    @test lsys.D[] == 1

    matrices_So, _ = get_sensitivity(sys_outer, sys_outer.inner.plant_output)
    lsyso = sminreal(ss(matrices_So...))
    @test lsys == lsyso || lsys == -1 * lsyso * (-1) # Output and input sensitivities are equal for SISO systems
end

@testset "multilevel system with loop openings" begin
    @named P_inner = FirstOrder(k = 1, T = 1)
    @named feedback = Feedback()
    @named ref = Step()
    @named sys_inner = ODESystem(
        [connect(P_inner.output, :y, feedback.input2)
         connect(feedback.output, :u, P_inner.input)
         connect(ref.output, :r, feedback.input1)],
        t,
        systems = [P_inner, feedback, ref])

    P_not_broken, _ = linearize(sys_inner, :u, :y)
    @test P_not_broken.A[] == -2
    P_broken, _ = linearize(sys_inner, :u, :y, loop_openings = [:u])
    @test P_broken.A[] == -1
    P_broken, _ = linearize(sys_inner, :u, :y, loop_openings = [:y])
    @test P_broken.A[] == -1

    Sinner = sminreal(ss(get_sensitivity(sys_inner, :u)[1]...))

    @named sys_inner = ODESystem(
        [connect(P_inner.output, :y, feedback.input2)
         connect(feedback.output, :u, P_inner.input)],
        t,
        systems = [P_inner, feedback])

    @named P_outer = FirstOrder(k = rand(), T = rand())

    @named sys_outer = ODESystem(
        [connect(sys_inner.P_inner.output, :y2, P_outer.input)
         connect(P_outer.output, :u2, sys_inner.feedback.input1)],
        t,
        systems = [P_outer, sys_inner])

    Souter = sminreal(ss(get_sensitivity(sys_outer, sys_outer.sys_inner.u)[1]...))

    Sinner2 = sminreal(ss(get_sensitivity(
        sys_outer, sys_outer.sys_inner.u, loop_openings = [:y2])[1]...))

    @test Sinner.nx == 1
    @test Sinner == Sinner2
    @test Souter.nx == 2
end

@testset "sensitivities in multivariate signals" begin
    A = [-0.994 -0.0794; -0.006242 -0.0134]
    B = [-0.181 -0.389; 1.1 1.12]
    C = [1.74 0.72; -0.33 0.33]
    D = [0.0 0.0; 0.0 0.0]
    @named P = Blocks.StateSpace(A, B, C, D)
    Pss = CS.ss(A, B, C, D)

    A = [-0.097;;]
    B = [-0.138 -1.02]
    C = [-0.076; 0.09;;]
    D = [0.0 0.0; 0.0 0.0]
    @named K = Blocks.StateSpace(A, B, C, D)
    Kss = CS.ss(A, B, C, D)

    eqs = [connect(P.output, :plant_output, K.input)
           connect(K.output, :plant_input, P.input)]
    sys = ODESystem(eqs, t, systems = [P, K], name = :hej)

    matrices, _ = get_sensitivity(sys, :plant_input)
    S = CS.feedback(I(2), Kss * Pss, pos_feedback = true)

    @test CS.tf(CS.ss(matrices...)) ≈ CS.tf(S)

    matrices, _ = get_comp_sensitivity(sys, :plant_input)
    T = -CS.feedback(Kss * Pss, I(2), pos_feedback = true)

    # bodeplot([ss(matrices...), T])
    @test CS.tf(CS.ss(matrices...)) ≈ CS.tf(T)

    matrices, _ = get_looptransfer(
        sys, :plant_input)
    L = Kss * Pss
    @test CS.tf(CS.ss(matrices...)) ≈ CS.tf(L)

    matrices, _ = linearize(sys, AnalysisPoint(:plant_input), :plant_output)
    G = CS.feedback(Pss, Kss, pos_feedback = true)
    @test CS.tf(CS.ss(matrices...)) ≈ CS.tf(G)
end

@testset "multiple analysis points" begin
    @named P = FirstOrder(k = 1, T = 1)
    @named C = Gain(; k = 1)
    @named add = Blocks.Add(k2 = -1)

    eqs = [connect(P.output, :plant_output, add.input2)
           connect(add.output, C.input)
           connect(C.output, :plant_input, P.input)]

    sys_inner = ODESystem(eqs, t, systems = [P, C, add], name = :inner)

    @named r = Constant(k = 1)
    @named F = FirstOrder(k = 1, T = 3)

    eqs = [connect(r.output, F.input)
           connect(F.output, sys_inner.add.input1)]
    sys_outer = ODESystem(eqs, t, systems = [F, sys_inner, r], name = :outer)

    matrices, _ = get_sensitivity(
        sys_outer, [sys_outer.inner.plant_input, sys_outer.inner.plant_output])

    Ps = tf(1, [1, 1]) |> ss
    Cs = tf(1) |> ss

    G = CS.ss(matrices...) |> sminreal
    Si = CS.feedback(1, Cs * Ps)
    @test tf(G[1, 1]) ≈ tf(Si)

    So = CS.feedback(1, Ps * Cs)
    @test tf(G[2, 2]) ≈ tf(So)
    @test tf(G[1, 2]) ≈ tf(-CS.feedback(Cs, Ps))
    @test tf(G[2, 1]) ≈ tf(CS.feedback(Ps, Cs))

    matrices, _ = get_comp_sensitivity(
        sys_outer, [sys_outer.inner.plant_input, sys_outer.inner.plant_output])

    G = CS.ss(matrices...) |> sminreal
    Ti = CS.feedback(Cs * Ps)
    @test tf(G[1, 1]) ≈ tf(Ti)

    To = CS.feedback(Ps * Cs)
    @test tf(G[2, 2]) ≈ tf(To)
    @test tf(G[1, 2]) ≈ tf(CS.feedback(Cs, Ps)) # The negative sign appears in a confusing place due to negative feedback not happening through Ps
    @test tf(G[2, 1]) ≈ tf(-CS.feedback(Ps, Cs))

    # matrices, _ = get_looptransfer(sys_outer, [:inner_plant_input, :inner_plant_output])
    matrices, _ = get_looptransfer(sys_outer, sys_outer.inner.plant_input)
    L = CS.ss(matrices...) |> sminreal
    @test tf(L) ≈ -tf(Cs * Ps)

    matrices, _ = get_looptransfer(sys_outer, sys_outer.inner.plant_output)
    L = CS.ss(matrices...) |> sminreal
    @test tf(L[1, 1]) ≈ -tf(Ps * Cs)

    # Calling looptransfer like below is not the intended way, but we can work out what it should return if we did so it remains a valid test
    matrices, _ = get_looptransfer(
        sys_outer, [sys_outer.inner.plant_input, sys_outer.inner.plant_output])
    L = CS.ss(matrices...) |> sminreal
    @test tf(L[1, 1]) ≈ tf(0)
    @test tf(L[2, 2]) ≈ tf(0)
    @test sminreal(L[1, 2]) ≈ ss(-1)
    @test tf(L[2, 1]) ≈ tf(Ps)

    matrices, _ = linearize(
        sys_outer, [sys_outer.inner.plant_input], [nameof(sys_inner.plant_output)])
    G = CS.ss(matrices...) |> sminreal
    @test tf(G) ≈ tf(CS.feedback(Ps, Cs))
end
