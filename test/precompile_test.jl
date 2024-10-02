using Test
using ModelingToolkit

using Distributed

# Test that the precompiled ODE system works
@everywhere push!(LOAD_PATH, joinpath(@__DIR__, "precompile_test"))

using ODEPrecompileTest

u = collect(1:3)
p = ModelingToolkit.MTKParameters(ODEPrecompileTest.f_noeval_good.sys,
    [:σ, :ρ, :β] .=> collect(4:6))

# These cases do not work, because they get defined in the ModelingToolkit's RGF cache.
@test parentmodule(typeof(ODEPrecompileTest.f_bad.f.f_iip).parameters[2]) == ModelingToolkit
@test parentmodule(typeof(ODEPrecompileTest.f_bad.f.f_oop).parameters[2]) == ModelingToolkit
@test parentmodule(typeof(ODEPrecompileTest.f_noeval_bad.f.f_iip).parameters[2]) ==
      ModelingToolkit
@test parentmodule(typeof(ODEPrecompileTest.f_noeval_bad.f.f_oop).parameters[2]) ==
      ModelingToolkit
@test_skip begin
    @test_throws KeyError ODEPrecompileTest.f_bad(u, p, 0.1)
    @test_throws KeyError ODEPrecompileTest.f_noeval_bad(u, p, 0.1)
end

# This case works, because it gets defined with the appropriate cache and context tags.
@test parentmodule(typeof(ODEPrecompileTest.f_noeval_good.f.f_iip).parameters[2]) ==
      ODEPrecompileTest
@test parentmodule(typeof(ODEPrecompileTest.f_noeval_good.f.f_oop).parameters[2]) ==
      ODEPrecompileTest
@test ODEPrecompileTest.f_noeval_good(u, p, 0.1) == [4, 0, -16]

ODEPrecompileTest.f_eval_bad(u, p, 0.1)

@test parentmodule(typeof(ODEPrecompileTest.f_eval_good.f.f_iip)) ==
      ODEPrecompileTest
@test parentmodule(typeof(ODEPrecompileTest.f_eval_good.f.f_oop)) ==
      ODEPrecompileTest
@test ODEPrecompileTest.f_eval_good(u, p, 0.1) == [4, 0, -16]
