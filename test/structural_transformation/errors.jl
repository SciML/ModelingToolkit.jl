using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Test

@mtkmodel FOL begin
    @parameters begin
        τ = 3.0 # parameters
    end
    @variables begin
        x(t) = 0.0 # dependent variables
    end
    @equations begin
        D(x) ~ (1 - x) / τ
    end
end

@named rc_model = ODESystem(rc_eqs, t; systems)
sys = structural_simplify(rc_model)
@test_throws ModelingToolkit.RepeatedStructuralSimplificationError structural_simplify(sys)