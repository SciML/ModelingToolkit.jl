using ModelingToolkit
using ModelingToolkit: SciMLBase
using Test

# Note: Full JET.jl analysis is too slow for regular CI testing.
# This file contains basic type stability tests that verify critical
# functions return expected types.

@testset "Type stability checks" begin
    # Basic system setup for testing
    @independent_variables t
    @variables x(t) y(t)
    @parameters p q
    D = Differential(t)

    # Create a simple linear ODE system
    eqs = [D(x) ~ p * x + q * y, D(y) ~ -q * x + p * y]

    @testset "System construction types" begin
        sys = ODESystem(eqs, t; name = :test_sys)
        @test sys isa System
    end

    @testset "structural_simplify return type" begin
        sys = ODESystem(eqs, t; name = :test_sys)
        simp_sys = structural_simplify(sys)
        @test simp_sys isa System
    end

    @testset "Accessor function types" begin
        sys = ODESystem(eqs, t; name = :test_sys)
        simp_sys = structural_simplify(sys)

        # Core accessor functions should return proper types
        @test equations(simp_sys) isa Vector{Equation}
        @test unknowns(simp_sys) isa Vector
        @test parameters(simp_sys) isa Vector
        @test observed(simp_sys) isa Vector{Equation}
    end

    @testset "ODEProblem construction" begin
        sys = ODESystem(eqs, t; name = :test_sys)
        simp_sys = structural_simplify(sys)

        # Test that problem construction returns correct type
        # Use the non-deprecated API: merge u0 and p into a single Dict
        op = merge(Dict(x => 1.0, y => 0.0), Dict(p => -0.1, q => 2.0))
        prob = ODEProblem(simp_sys, op, (0.0, 1.0))
        @test prob isa SciMLBase.ODEProblem
    end
end
