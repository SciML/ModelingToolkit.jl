using Test
using ModelingToolkit

# Test that the precompiled ODE system works
push!(LOAD_PATH, joinpath(@__DIR__, "precompile_test"))
using ODEPrecompileTest

u  = collect(1:3)
p  = collect(4:6)

# This case does not work, because "f_bad" gets defined in ModelingToolkit
# instead of in the compiled module!
@test parentmodule(typeof(ODEPrecompileTest.f_bad.f.f_iip).parameters[2]) == ModelingToolkit
@test parentmodule(typeof(ODEPrecompileTest.f_bad.f.f_oop).parameters[2]) == ModelingToolkit
@test_throws KeyError ODEPrecompileTest.f_bad(u, p, 0.1)

# This case works, because "f_good" gets defined in the precompiled module.
@test parentmodule(typeof(ODEPrecompileTest.f_good.f.f_iip).parameters[2]) == ODEPrecompileTest
@test parentmodule(typeof(ODEPrecompileTest.f_good.f.f_oop).parameters[2]) == ODEPrecompileTest
@test ODEPrecompileTest.f_good(u, p, 0.1) == [4, 0, -16]