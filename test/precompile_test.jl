using Test
using ModelingToolkit

# Test that the precompiled ODE system works
push!(LOAD_PATH, joinpath(@__DIR__, "precompile_test"))
using ODEPrecompileTest

du = zeros(3)
u  = collect(1:3)
p  = collect(4:6)

# This case does not work, because the function gets defined in ModelingToolkit
# instead of in the compiled module!
@test parentmodule(typeof(ODEPrecompileTest.f_bad.f.f_iip).parameters[2]) == ModelingToolkit
@test parentmodule(typeof(ODEPrecompileTest.f_bad.f.f_oop).parameters[2]) == ModelingToolkit
@test_throws KeyError ODEPrecompileTest.f_bad(du, u, p, 0.1)

# This case works, because the function gets defined in the compiled module.
# @test parentmodule(typeof(ODEPrecompileTest.f_good.f.f_iip).parameters[2]) == ODEPrecompileTest
# @test parentmodule(typeof(ODEPrecompileTest.f_good.f.f_oop).parameters[2]) == ODEPrecompileTest
# @test ODEPrecompileTest.f_good(du, u, p, 0.1)