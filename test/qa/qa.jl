using ModelingToolkit
using Aqua
using JET
using SciMLTesting
using TOML
using Test

@testset "Downstream LinearSolve compatibility" begin
    root_project = TOML.parsefile(joinpath(@__DIR__, "..", "..", "Project.toml"))
    downstream_project = TOML.parsefile(joinpath(@__DIR__, "..", "downstream", "Project.toml"))

    @test downstream_project["deps"]["LinearSolve"] == root_project["extras"]["LinearSolve"]
    @test downstream_project["compat"]["LinearSolve"] == root_project["compat"]["LinearSolve"]
end

run_qa(
    ModelingToolkit;
    Aqua = Aqua,
    JET = JET,
    jet = true,
    jet_kwargs = (; target_defined_modules = true),
)
