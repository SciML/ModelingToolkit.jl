using ModelingToolkit, ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq, LinearAlgebra
using Test
using ModelingToolkit: t_nounits as t, D_nounits as D, AnalysisPoint, AbstractSystem
import ModelingToolkit as MTK
using Symbolics: NAMESPACE_SEPARATOR

@testset "AnalysisPoint is lowered to `connect`" begin
    @named P = FirstOrder(k = 1, T = 1)
    @named C = Gain(; k = -1)

    ap = AnalysisPoint(:plant_input)
    eqs = [connect(P.output, C.input)
           connect(C.output, ap, P.input)]
    sys_ap = ODESystem(eqs, t, systems = [P, C], name = :hej)
    sys_ap2 = @test_nowarn expand_connections(sys_ap)

    @test all(eq -> !(eq.lhs isa AnalysisPoint), equations(sys_ap2))

    eqs = [connect(P.output, C.input)
           connect(C.output, P.input)]
    sys_normal = ODESystem(eqs, t, systems = [P, C], name = :hej)
    sys_normal2 = @test_nowarn expand_connections(sys_normal)

    @test isequal(sys_ap2, sys_normal2)
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
sys = ODESystem(eqs, t, systems = [P, C], name = :hej)
@named nested_sys = ODESystem(Equation[], t; systems = [sys])

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
