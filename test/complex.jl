using ModelingToolkit
using ModelingToolkit: t_nounits as t
using Test

@mtkmodel ComplexModel begin
    @variables begin
        x(t)
        y(t)
        z(t)::Complex
    end
    @equations begin
        z ~ x + im * y
    end
end
@named mixed = ComplexModel()
@test length(equations(mixed)) == 2
