using ModelingToolkit, CommonSolve, SymbolicIndexingInterface, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

@testset "Input derivative diagnostics" begin
    @variables x(t) = 1.0 y(t) = 1.0 u(t) = 2.0 v(t) = 0.0
    model = System(
        [D(x) ~ -x + y + v, 0 ~ x^2 + y^2 - u], t; name = :plant
    )
    lf, _ = linearization_function(
        model, [u, v], [x]; op = Dict(x => 1.0, y => 1.0, u => 2.0, v => 0.0),
        initialize = false
    )
    prob = LinearizationProblem(lf, 0.0)
    err = try
        solve(prob)
        nothing
    catch exception
        exception
    end
    @test err isa ModelingToolkit.InputDerivativeError
    if err isa ModelingToolkit.InputDerivativeError
        @test isequal(err.inputs, lf.inputs[1:1])
        @test err.num_inputs == 2
        message = sprint(showerror, err)
        @test occursin("[u(t)]", message)
        @test occursin("allow_input_derivatives = true", message)
        @test occursin("double width (4)", message)
        @test occursin("last 2 inputs", message)
        @test occursin("first 2 inputs", message)
    end
    matrices, _ = solve(prob; allow_input_derivatives = true)
    @test size(matrices.B, 2) == 4
    @test size(matrices.D, 2) == 4
    @test any(!iszero, matrices.B[:, 3:4])

    both = ModelingToolkit.InputDerivativeError(lf.inputs, 2)
    @test occursin("[u(t), v(t)]", sprint(showerror, both))
end
