using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Test

@testset "tspan field storage" begin
    @parameters σ ρ β
    @variables x(t) y(t) z(t)
    eqs = [D(x) ~ σ*(y-x), D(y) ~ x*(ρ-z)-y, D(z) ~ x*y - β*z]

    # Test 1: Check if tspan is stored correctly
    @named sys = ODESystem(eqs, t, [x,y,z], [σ,ρ,β]; tspan=(0.0, 10.0))
    @test getfield(sys, :tspan) == (0.0, 10.0)

    # Test 2: Check if it works without tspan (backward compatibility)
    @named sys_no_tspan = ODESystem(eqs, t, [x,y,z], [σ,ρ,β])
    @test getfield(sys_no_tspan, :tspan) === nothing
end