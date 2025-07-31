using ModelingToolkit, ModelingToolkitStandardLibrary.Blocks, ControlSystemsBase
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq, LinearAlgebra
using Test
using ModelingToolkit: t_nounits as t, D_nounits as D, AnalysisPoint, AbstractSystem
import ModelingToolkit as MTK
import ControlSystemsBase as CS
using Symbolics: NAMESPACE_SEPARATOR

@testset "AnalysisPoint is lowered to `connect`" begin
    @named P = FirstOrder(k = 1, T = 1)
    @named C = Gain(; k = -1)

    ap = AnalysisPoint(:plant_input)
    eqs = [connect(P.output, C.input)
           connect(C.output, ap, P.input)]
    sys_ap = System(eqs, t, systems = [P, C], name = :hej)
    sys_ap2 = @test_nowarn expand_connections(sys_ap)

    @test all(eq -> !(eq.lhs isa AnalysisPoint), equations(sys_ap2))

    eqs = [connect(P.output, C.input)
           connect(C.output, P.input)]
    sys_normal = System(eqs, t, systems = [P, C], name = :hej)
    sys_normal2 = @test_nowarn expand_connections(sys_normal)

    @test issetequal(equations(sys_ap2), equations(sys_normal2))
end

@testset "Inverse causality throws a warning" begin
    @named P = FirstOrder(k = 1, T = 1)
    @named C = Gain(; k = -1)

    ap = AnalysisPoint(:plant_input)
    @test_warn ["1-th argument", "plant_input", "not a output"] connect(
        P.input, ap, C.output)
    @test_nowarn connect(P.input, ap, C.output; verbose = false)
end

# also tests `connect(input, name::Symbol, outputs...)` syntax
@testset "AnalysisPoint is accessible via `getproperty`" begin
    @named P = FirstOrder(k = 1, T = 1)
    @named C = Gain(; k = -1)

    eqs = [connect(P.output, C.input), connect(C.output, :plant_input, P.input)]
    sys_ap = System(eqs, t, systems = [P, C], name = :hej)
    ap2 = @test_nowarn sys_ap.plant_input
    @test nameof(ap2) == Symbol(join(["hej", "plant_input"], NAMESPACE_SEPARATOR))
    @named sys = System(Equation[], t; systems = [sys_ap])
    ap3 = @test_nowarn sys.hej.plant_input
    @test nameof(ap3) == Symbol(join(["sys", "hej", "plant_input"], NAMESPACE_SEPARATOR))
    csys = complete(sys)
    ap4 = csys.hej.plant_input
    @test nameof(ap4) == Symbol(join(["hej", "plant_input"], NAMESPACE_SEPARATOR))
    nsys = toggle_namespacing(sys, false)
    ap5 = nsys.hej.plant_input
    @test nameof(ap4) == Symbol(join(["hej", "plant_input"], NAMESPACE_SEPARATOR))
end

### Ported from MTKStdlib

@named P = FirstOrder(k = 1, T = 1)
@named C = Gain(; k = -1)

ap = AnalysisPoint(:plant_input)
eqs = [connect(P.output, C.input), connect(C.output, ap, P.input)]
sys = System(eqs, t, systems = [P, C], name = :hej)
@named nested_sys = System(Equation[], t; systems = [sys])
nonamespace_sys = toggle_namespacing(nested_sys, false)

@testset "simplifies and solves" begin
    ssys = mtkcompile(sys)
    prob = ODEProblem(ssys, [P.x => 1], (0, 10))
    sol = solve(prob, Rodas5())
    @test norm(sol.u[1]) >= 1
    @test norm(sol.u[end]) < 1e-6 # This fails without the feedback through C
end

test_cases = [
    ("inner", sys, sys.plant_input),
    ("nested", nested_sys, nested_sys.hej.plant_input),
    ("nonamespace", nonamespace_sys, nonamespace_sys.hej.plant_input),
    ("inner - Symbol", sys, :plant_input),
    ("nested - Symbol", nested_sys, nameof(sys.plant_input))
]

@testset "get_sensitivity - $name" for (name, sys, ap) in test_cases
    matrices, _ = get_sensitivity(sys, ap)
    @test matrices.A[] == -2
    @test matrices.B[] * matrices.C[] == -1 # either one negative
    @test matrices.D[] == 1
end

@testset "get_comp_sensitivity - $name" for (name, sys, ap) in test_cases
    matrices, _ = get_comp_sensitivity(sys, ap)
    @test matrices.A[] == -2
    @test matrices.B[] * matrices.C[] == 1 # both positive or negative
    @test matrices.D[] == 0
end

#=
# Equivalent code using ControlSystems. This can be used to verify the expected results tested for above.
using ControlSystemsBase
P = tf(1.0, [1, 1])
C = 1                      # Negative feedback assumed in ControlSystems
S = sensitivity(P, C)      # or feedback(1, P*C)
T = comp_sensitivity(P, C) # or feedback(P*C)
=#

@testset "get_looptransfer - $name" for (name, sys, ap) in test_cases
    matrices, _ = get_looptransfer(sys, ap)
    @test matrices.A[] == -1
    @test matrices.B[] * matrices.C[] == -1 # either one negative
    @test matrices.D[] == 0
end

#=
# Equivalent code using ControlSystems. This can be used to verify the expected results tested for above.
using ControlSystemsBase
P = tf(1.0, [1, 1])
C = -1
L = P*C
=#

@testset "open_loop - $name" for (name, sys, ap) in test_cases
    open_sys, (du, u) = open_loop(sys, ap)
    matrices, _ = linearize(open_sys, [du], [u])
    @test matrices.A[] == -1
    @test matrices.B[] * matrices.C[] == -1 # either one negative
    @test matrices.D[] == 0
end

# Multiple analysis points

eqs = [connect(P.output, :plant_output, C.input)
       connect(C.output, :plant_input, P.input)]
sys = System(eqs, t, systems = [P, C], name = :hej)
@named nested_sys = System(Equation[], t; systems = [sys])

test_cases = [
    ("inner", sys, sys.plant_input),
    ("nested", nested_sys, nested_sys.hej.plant_input),
    ("inner - Symbol", sys, :plant_input),
    ("nested - Symbol", nested_sys, nameof(sys.plant_input))
]

@testset "get_sensitivity - $name" for (name, sys, ap) in test_cases
    matrices, _ = get_sensitivity(sys, ap)
    @test matrices.A[] == -2
    @test matrices.B[] * matrices.C[] == -1 # either one negative
    @test matrices.D[] == 1
end

@testset "linearize - $name" for (name, sys, inputap, outputap) in [
    ("inner", sys, sys.plant_input, sys.plant_output),
    ("nested", nested_sys, nested_sys.hej.plant_input, nested_sys.hej.plant_output)
]
    inputname = Symbol(join(
        MTK.namespace_hierarchy(nameof(inputap))[2:end], NAMESPACE_SEPARATOR))
    outputname = Symbol(join(
        MTK.namespace_hierarchy(nameof(outputap))[2:end], NAMESPACE_SEPARATOR))
    @testset "input - $(typeof(input)), output - $(typeof(output))" for (input, output) in [
        (inputap, outputap),
        (inputname, outputap),
        (inputap, outputname),
        (inputname, outputname),
        (inputap, [outputap]),
        (inputname, [outputap]),
        (inputap, [outputname]),
        (inputname, [outputname])
    ]
        matrices, _ = linearize(sys, input, output)
        # Result should be the same as feedpack(P, 1), i.e., the closed-loop transfer function from plant input to plant output
        @test matrices.A[] == -2
        @test matrices.B[] * matrices.C[] == 1 # both positive
        @test matrices.D[] == 0
    end
end

@testset "linearize - variable output - $name" for (name, sys, input, output) in [
    ("inner", sys, sys.plant_input, P.output.u),
    ("nested", nested_sys, nested_sys.hej.plant_input, sys.P.output.u)
]
    matrices, _ = linearize(sys, input, [output])
    @test matrices.A[] == -2
    @test matrices.B[] * matrices.C[] == 1 # both positive
    @test matrices.D[] == 0
end

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
            return System(eqs, t;
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
        System(eqs, t; systems = [torque, inertia1, inertia2, spring, damper], name)
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

    closed_loop = System(connections, t, systems = [model, pid, filt, sensor, r, er],
        name = :closed_loop, defaults = [
            model.inertia1.phi => 0.0,
            model.inertia2.phi => 0.0,
            model.inertia1.w => 0.0,
            model.inertia2.w => 0.0,
            filt.x => 0.0,
            filt.xd => 0.0
        ])

    sys = mtkcompile(closed_loop)
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

@testset "Duplicate `connect` statements across subsystems with AP transforms - standard `connect`" begin
    @named P = FirstOrder(k = 1, T = 1)
    @named C = Gain(; k = 1)
    @named add = Blocks.Add(k2 = -1)

    eqs = [connect(P.output, :plant_output, add.input2)
           connect(add.output, C.input)
           connect(C.output, P.input)]

    sys_inner = System(eqs, t, systems = [P, C, add], name = :inner)

    @named r = Constant(k = 1)
    @named F = FirstOrder(k = 1, T = 3)

    eqs = [connect(r.output, F.input)
           connect(sys_inner.P.output, sys_inner.add.input2)
           connect(sys_inner.C.output, :plant_input, sys_inner.P.input)
           connect(F.output, sys_inner.add.input1)]
    sys_outer = System(eqs, t, systems = [F, sys_inner, r], name = :outer)

    # test first that the mtkcompile works correctly
    ssys = mtkcompile(sys_outer)
    prob = ODEProblem(ssys, Pair[], (0, 10))
    @test_nowarn solve(prob, Rodas5())

    matrices, _ = get_sensitivity(sys_outer, sys_outer.plant_input)
    lsys = sminreal(ss(matrices...))
    @test lsys.A[] == -2
    @test lsys.B[] * lsys.C[] == -1 # either one negative
    @test lsys.D[] == 1

    matrices_So, _ = get_sensitivity(sys_outer, sys_outer.inner.plant_output)
    lsyso = sminreal(ss(matrices_So...))
    @test lsys == lsyso || lsys == -1 * lsyso * (-1) # Output and input sensitivities are equal for SISO systems
end

@testset "Duplicate `connect` statements across subsystems with AP transforms - causal variable `connect`" begin
    @named P = FirstOrder(k = 1, T = 1)
    @named C = Gain(; k = 1)
    @named add = Blocks.Add(k2 = -1)

    eqs = [connect(P.output.u, :plant_output, add.input2.u)
           connect(add.output, C.input)
           connect(C.output.u, P.input.u)]

    sys_inner = System(eqs, t, systems = [P, C, add], name = :inner)

    @named r = Constant(k = 1)
    @named F = FirstOrder(k = 1, T = 3)

    eqs = [connect(r.output, F.input)
           connect(sys_inner.P.output.u, sys_inner.add.input2.u)
           connect(sys_inner.C.output.u, :plant_input, sys_inner.P.input.u)
           connect(F.output, sys_inner.add.input1)]
    sys_outer = System(eqs, t, systems = [F, sys_inner, r], name = :outer)

    # test first that the mtkcompile works correctly
    ssys = mtkcompile(sys_outer)
    prob = ODEProblem(ssys, Pair[], (0, 10))
    @test_nowarn solve(prob, Rodas5())

    matrices, _ = get_sensitivity(sys_outer, sys_outer.plant_input)
    lsys = sminreal(ss(matrices...))
    @test lsys.A[] == -2
    @test lsys.B[] * lsys.C[] == -1 # either one negative
    @test lsys.D[] == 1

    matrices_So, _ = get_sensitivity(sys_outer, sys_outer.inner.plant_output)
    lsyso = sminreal(ss(matrices_So...))
    @test lsys == lsyso || lsys == -1 * lsyso * (-1) # Output and input sensitivities are equal for SISO systems
end

@testset "Duplicate `connect` statements across subsystems with AP transforms - mixed `connect`" begin
    @named P = FirstOrder(k = 1, T = 1)
    @named C = Gain(; k = 1)
    @named add = Blocks.Add(k2 = -1)

    eqs = [connect(P.output.u, :plant_output, add.input2.u)
           connect(add.output, C.input)
           connect(C.output, P.input)]

    sys_inner = System(eqs, t, systems = [P, C, add], name = :inner)

    @named r = Constant(k = 1)
    @named F = FirstOrder(k = 1, T = 3)

    eqs = [connect(r.output, F.input)
           connect(sys_inner.P.output, sys_inner.add.input2)
           connect(sys_inner.C.output.u, :plant_input, sys_inner.P.input.u)
           connect(F.output, sys_inner.add.input1)]
    sys_outer = System(eqs, t, systems = [F, sys_inner, r], name = :outer)

    # test first that the mtkcompile works correctly
    ssys = mtkcompile(sys_outer)
    prob = ODEProblem(ssys, Pair[], (0, 10))
    @test_nowarn solve(prob, Rodas5())

    matrices, _ = get_sensitivity(sys_outer, sys_outer.plant_input)
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
    @named sys_inner = System(
        [connect(P_inner.output, :y, feedback.input2)
         connect(feedback.output, :u, P_inner.input)
         connect(ref.output, :r, feedback.input1)],
        t,
        systems = [P_inner, feedback, ref])

    P_not_broken, _ = linearize(sys_inner, :u, :y)
    @test P_not_broken.A[] == -2
    P_broken, ssys = linearize(sys_inner, :u, :y, loop_openings = [:u])
    @test isequal(defaults(ssys)[ssys.d_u], ssys.feedback.output.u)
    @test P_broken.A[] == -1
    P_broken, ssys = linearize(sys_inner, :u, :y, loop_openings = [:y])
    @test isequal(defaults(ssys)[ssys.d_y], ssys.P_inner.output.u)
    @test P_broken.A[] == -1

    Sinner = sminreal(ss(get_sensitivity(sys_inner, :u)[1]...))

    @named sys_inner = System(
        [connect(P_inner.output, :y, feedback.input2)
         connect(feedback.output, :u, P_inner.input)],
        t,
        systems = [P_inner, feedback])

    @named P_outer = FirstOrder(k = rand(), T = rand())

    @named sys_outer = System(
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
    sys = System(eqs, t, systems = [P, K], name = :hej)

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

    sys_inner = System(eqs, t, systems = [P, C, add], name = :inner)

    @named r = Constant(k = 1)
    @named F = FirstOrder(k = 1, T = 3)

    eqs = [connect(r.output, F.input)
           connect(F.output, sys_inner.add.input1)]
    sys_outer = System(eqs, t, systems = [F, sys_inner, r], name = :outer)

    matrices,
    _ = get_sensitivity(
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

    matrices,
    _ = get_comp_sensitivity(
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
    matrices,
    _ = get_looptransfer(
        sys_outer, [sys_outer.inner.plant_input, sys_outer.inner.plant_output])
    L = CS.ss(matrices...) |> sminreal
    @test tf(L[1, 1]) ≈ tf(0)
    @test tf(L[2, 2]) ≈ tf(0)
    @test sminreal(L[1, 2]) ≈ ss(-1)
    @test tf(L[2, 1]) ≈ tf(Ps)

    matrices,
    _ = linearize(
        sys_outer, [sys_outer.inner.plant_input], [nameof(sys_inner.plant_output)])
    G = CS.ss(matrices...) |> sminreal
    @test tf(G) ≈ tf(CS.feedback(Ps, Cs))
end

function normal_test_system()
    @named F1 = FirstOrder(k = 1, T = 1)
    @named F2 = FirstOrder(k = 1, T = 1)
    @named add = Blocks.Add(k1 = 1, k2 = 2)
    @named back = Feedback()

    eqs_normal = [connect(back.output, :ap, F1.input)
                  connect(back.output, F2.input)
                  connect(F1.output, add.input1)
                  connect(F2.output, add.input2)
                  connect(add.output, back.input2)]
    @named normal_inner = System(eqs_normal, t; systems = [F1, F2, add, back])

    @named step = Step()
    eqs2_normal = [
        connect(step.output, normal_inner.back.input1)
    ]
    @named sys_normal = System(eqs2_normal, t; systems = [normal_inner, step])
end

sys_normal = normal_test_system()

prob = ODEProblem(mtkcompile(sys_normal), [], (0.0, 10.0))
@test SciMLBase.successful_retcode(solve(prob, Rodas5P()))
matrices_normal, _ = get_sensitivity(sys_normal, sys_normal.normal_inner.ap)

@testset "Analysis point overriding part of connection - normal connect" begin
    @named F1 = FirstOrder(k = 1, T = 1)
    @named F2 = FirstOrder(k = 1, T = 1)
    @named add = Blocks.Add(k1 = 1, k2 = 2)
    @named back = Feedback()

    eqs = [connect(back.output, F1.input, F2.input)
           connect(F1.output, add.input1)
           connect(F2.output, add.input2)
           connect(add.output, back.input2)]
    @named inner = System(eqs, t; systems = [F1, F2, add, back])

    @named step = Step()
    eqs2 = [connect(step.output, inner.back.input1)
            connect(inner.back.output, :ap, inner.F1.input)]
    @named sys = System(eqs2, t; systems = [inner, step])

    prob = ODEProblem(mtkcompile(sys), [], (0.0, 10.0))
    @test SciMLBase.successful_retcode(solve(prob, Rodas5P()))

    matrices, _ = get_sensitivity(sys, sys.ap)
    @test matrices == matrices_normal
end

@testset "Analysis point overriding part of connection - variable connect" begin
    @named F1 = FirstOrder(k = 1, T = 1)
    @named F2 = FirstOrder(k = 1, T = 1)
    @named add = Blocks.Add(k1 = 1, k2 = 2)
    @named back = Feedback()

    eqs = [connect(back.output.u, F1.input.u, F2.input.u)
           connect(F1.output, add.input1)
           connect(F2.output, add.input2)
           connect(add.output, back.input2)]
    @named inner = System(eqs, t; systems = [F1, F2, add, back])

    @named step = Step()
    eqs2 = [connect(step.output, inner.back.input1)
            connect(inner.back.output.u, :ap, inner.F1.input.u)]
    @named sys = System(eqs2, t; systems = [inner, step])

    prob = ODEProblem(mtkcompile(sys), [], (0.0, 10.0))
    @test SciMLBase.successful_retcode(solve(prob, Rodas5P()))

    matrices, _ = get_sensitivity(sys, sys.ap)
    @test matrices == matrices_normal
end

@testset "Analysis point overriding part of connection - mixed connect" begin
    @named F1 = FirstOrder(k = 1, T = 1)
    @named F2 = FirstOrder(k = 1, T = 1)
    @named add = Blocks.Add(k1 = 1, k2 = 2)
    @named back = Feedback()

    eqs = [connect(back.output, F1.input, F2.input)
           connect(F1.output, add.input1)
           connect(F2.output, add.input2)
           connect(add.output, back.input2)]
    @named inner = System(eqs, t; systems = [F1, F2, add, back])

    @named step = Step()
    eqs2 = [connect(step.output, inner.back.input1)
            connect(inner.back.output.u, :ap, inner.F1.input.u)]
    @named sys = System(eqs2, t; systems = [inner, step])

    prob = ODEProblem(mtkcompile(sys), [], (0.0, 10.0))
    @test SciMLBase.successful_retcode(solve(prob, Rodas5P()))

    matrices, _ = get_sensitivity(sys, sys.ap)
    @test matrices == matrices_normal
end

@testset "Ignored analysis points only affect relevant connection sets" begin
    m1 = 1
    m2 = 1
    k = 1000 # Spring stiffness
    c = 10   # Damping coefficient

    @named inertia1 = Inertia(; J = m1, phi = 0, w = 0)
    @named inertia2 = Inertia(; J = m2, phi = 0, w = 0)

    @named spring = Spring(; c = k)
    @named damper = Damper(; d = c)

    @named torque = Torque(use_support = false)

    function SystemModel(u = nothing; name = :model)
        eqs = [connect(torque.flange, inertia1.flange_a)
               connect(inertia1.flange_b, spring.flange_a, damper.flange_a)
               connect(inertia2.flange_a, spring.flange_b, damper.flange_b)]
        if u !== nothing
            push!(eqs, connect(torque.tau, u.output))
            return @named model = System(
                eqs, t; systems = [torque, inertia1, inertia2, spring, damper, u])
        end
        System(eqs, t; systems = [torque, inertia1, inertia2, spring, damper], name)
    end

    @named r = Step(start_time = 1)
    @named pid = LimPID(k = 400, Ti = 0.5, Td = 1, u_max = 350)
    @named filt = SecondOrder(d = 0.9, w = 10, x = 0, xd = 0)
    @named sensor = AngleSensor()
    @named add = Add() # To add the feedback and feedforward control signals
    model = SystemModel()
    @named inverse_model = SystemModel()
    @named inverse_sensor = AngleSensor()
    connections = [connect(r.output, :r, filt.input) # Name connection r to form an analysis point
                   connect(inverse_model.inertia1.flange_b, inverse_sensor.flange) # Attach the inverse sensor to the inverse model
                   connect(filt.output, pid.reference, inverse_sensor.phi) # the filtered reference now goes to both the PID controller and the inverse model input
                   connect(inverse_model.torque.tau, add.input1)
                   connect(pid.ctr_output, add.input2)
                   connect(add.output, :u, model.torque.tau) # Name connection u to form an analysis point
                   connect(model.inertia1.flange_b, sensor.flange)
                   connect(sensor.phi, :y, pid.measurement)]
    closed_loop = System(connections, t,
        systems = [model, inverse_model, pid, filt, sensor, inverse_sensor, r, add],
        name = :closed_loop)
    # just ensure the system simplifies
    mats, _ = get_sensitivity(closed_loop, :y)
    S = CS.ss(mats...)
    fr = CS.freqrespv(S, [0.01, 1, 100])
    # https://github.com/SciML/ModelingToolkit.jl/pull/3469
    reference_fr = ComplexF64[-1.2505330104772838e-11 - 2.500062613816021e-9im,
    -0.0024688370221621625 - 0.002279011866413123im,
    1.8100018764334602 + 0.3623845793211718im]
    @test isapprox(fr, reference_fr)
end
