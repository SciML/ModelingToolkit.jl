using ModelingToolkit
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(ModelingToolkit)
end

@testset "JET" begin
    JET.test_package(ModelingToolkit; target_defined_modules = true)
end
