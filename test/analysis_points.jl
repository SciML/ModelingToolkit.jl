using ModelingToolkit, ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq, LinearAlgebra
using Test
using ModelingToolkit: t_nounits as t, D_nounits as D, AnalysisPoint, get_sensitivity,
                       get_comp_sensitivity, get_looptransfer
using Symbolics: NAMESPACE_SEPARATOR

@testset "AnalysisPoint is lowered to `connect`" begin
    @named P = FirstOrder(k = 1, T = 1)
    @named C = Gain(; k = -1)

    ap = AnalysisPoint(:plant_input)
    eqs = [connect(P.output, C.input)
           connect(C.output, ap, P.input)]
    sys_ap = ODESystem(eqs, t, systems = [P, C], name = :hej)
    sys_ap2 = @test_nowarn expand_connections(sys)

    @test all(eq -> !(eq.lhs isa AnalysisPoint), equations(sys_ap2))

    eqs = [connect(P.output, C.input)
           connect(C.output, P.input)]
    sys_normal = ODESystem(eqs, t, systems = [P, C], name = :hej)
    sys_normal2 = @test_nowarn expand_connections(sys)

    @test isequal(sys_ap2, sys_normal2)
end

# also tests `connect(input, name::Symbol, outputs...)` syntax
@testset "AnalysisPoint is accessible via `getproperty`" begin
    @named P = FirstOrder(k = 1, T = 1)
    @named C = Gain(; k = -1)

    eqs = [connect(P.output, C.input), connect(C.output, :plant_input, P.input)]
    sys_ap = ODESystem(eqs, t, systems = [P, C], name = :hej)
    ap2 = @test_nowarn sys_ap.plant_input
    @test nameof(ap2) == Symbol(join(["hej", "plant_input"], NAMESPACE_SEPARATOR))
    @named sys = ODESystem(Equation[], t; systems = [sys_ap])
    ap3 = @test_nowarn sys.hej.plant_input
    @test nameof(ap3) == Symbol(join(["sys", "hej", "plant_input"], NAMESPACE_SEPARATOR))
    sys = complete(sys)
    ap4 = sys.hej.plant_input
    @test nameof(ap4) == Symbol(join(["hej", "plant_input"], NAMESPACE_SEPARATOR))
end

### Ported from MTKStdlib

@named P = FirstOrder(k = 1, T = 1)
@named C = Gain(; k = -1)

ap = AnalysisPoint(:plant_input)
eqs = [connect(P.output, C.input), connect(C.output, ap, P.input)]
sys = ODESystem(eqs, t, systems = [P, C], name = :hej)
@named nested_sys = ODESystem(Equation[], t; systems = [sys])

@testset "simplifies and solves" begin
    ssys = structural_simplify(sys)
    prob = ODEProblem(ssys, [P.x => 1], (0, 10))
    sol = solve(prob, Rodas5())
    @test norm(sol.u[1]) >= 1
    @test norm(sol.u[end]) < 1e-6 # This fails without the feedback through C
end

@testset "get_sensitivity - $name" for (name, sys, ap) in [
    ("inner", sys, sys.plant_input),
    ("nested", nested_sys, nested_sys.hej.plant_input),
    ("inner - nonamespace", sys, :plant_input),
    ("inner - Symbol", sys, nameof(sys.plant_input)),
    ("nested - Symbol", nested_sys, nameof(nested_sys.hej.plant_input))
]
    matrices, _ = get_sensitivity(sys, ap)
    @test matrices.A[] == -2
    @test matrices.B[] * matrices.C[] == -1 # either one negative
    @test matrices.D[] == 1
end

@testset "get_comp_sensitivity - $name" for (name, sys, ap) in [
    ("inner", sys, sys.plant_input),
    ("nested", nested_sys, nested_sys.hej.plant_input),
    ("inner - nonamespace", sys, :plant_input),
    ("inner - Symbol", sys, nameof(sys.plant_input)),
    ("nested - Symbol", nested_sys, nameof(nested_sys.hej.plant_input))
]
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

@testset "get_looptransfer - $name" for (name, sys, ap) in [
    ("inner", sys, sys.plant_input),
    ("nested", nested_sys, nested_sys.hej.plant_input),
    ("inner - nonamespace", sys, :plant_input),
    ("inner - Symbol", sys, nameof(sys.plant_input)),
    ("nested - Symbol", nested_sys, nameof(nested_sys.hej.plant_input))
]
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

@testset "open_loop - $name" for (name, sys, ap) in [
    ("inner", sys, sys.plant_input),
    ("nested", nested_sys, nested_sys.hej.plant_input),
    ("inner - nonamespace", sys, :plant_input),
    ("inner - Symbol", sys, nameof(sys.plant_input)),
    ("nested - Symbol", nested_sys, nameof(nested_sys.hej.plant_input))
]
    open_sys, (du, u) = open_loop(sys, ap)
    matrices, _ = linearize(open_sys, [du], [u])
    @test matrices.A[] == -1
    @test matrices.B[] * matrices.C[] == -1 # either one negative
    @test matrices.D[] == 0
end

# Multiple analysis points

eqs = [connect(P.output, :plant_output, C.input)
       connect(C.output, :plant_input, P.input)]
sys = ODESystem(eqs, t, systems = [P, C], name = :hej)
@named nested_sys = ODESystem(Equation[], t; systems = [sys])

@testset "get_sensitivity - $name" for (name, sys, ap) in [
    ("inner", sys, sys.plant_input),
    ("nested", nested_sys, nested_sys.hej.plant_input),
    ("inner - nonamespace", sys, :plant_input),
    ("inner - Symbol", sys, nameof(sys.plant_input)),
    ("nested - Symbol", nested_sys, nameof(nested_sys.hej.plant_input))
]
    matrices, _ = get_sensitivity(sys, ap)
    @test matrices.A[] == -2
    @test matrices.B[] * matrices.C[] == -1 # either one negative
    @test matrices.D[] == 1
end

@testset "linearize - $name" for (name, sys, inputap, outputap) in [
    ("inner", sys, sys.plant_input, sys.plant_output),
    ("nested", nested_sys, nested_sys.hej.plant_input, nested_sys.hej.plant_output)
]
    @testset "input - $(typeof(input)), output - $(typeof(output))" for (input, output) in [
        (inputap, outputap),
        (nameof(inputap), outputap),
        (inputap, nameof(outputap)),
        (nameof(inputap), nameof(outputap)),
        (inputap, [outputap]),
        (nameof(inputap), [outputap]),
        (inputap, [nameof(outputap)]),
        (nameof(inputap), [nameof(outputap)])
    ]
        if input isa Symbol
            # broken because MTKStdlib defines the method for
            # `input::Union{Symbol, Vector{Symbol}}` which we can't directly call
            @test_broken linearize(sys, input, output)
            linfun, ssys = @invoke linearization_function(sys::AbstractSystem,
                input::Union{Symbol, Vector{Symbol}, AnalysisPoint, Vector{AnalysisPoint}},
                output::Any)
            matrices = linearize(ssys, linfun)
        else
            matrices, _ = linearize(sys, input, output)
        end
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

using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks: Sine, PID, SecondOrder, Step, RealOutput
using ModelingToolkit: connect

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
end
