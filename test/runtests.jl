using SciMLTesting
# https://github.com/JuliaLang/julia/issues/54664
import REPL

const MTKBasePath = joinpath(dirname(@__DIR__), "lib", "ModelingToolkitBase")
const RepoRoot = dirname(@__DIR__)
const EnvParents = [MTKBasePath, RepoRoot]

run_tests(;
    groups = Dict(
        "InterfaceI" => joinpath(@__DIR__, "group_interfacei.jl"),
        "Initialization" => joinpath(@__DIR__, "group_initialization.jl"),
        "InterfaceII" => joinpath(@__DIR__, "group_interfaceii.jl"),
        "SymbolicIndexingInterface" => joinpath(@__DIR__, "group_symbolicindexinginterface.jl"),
        "Downstream" => (;
            env = joinpath(@__DIR__, "downstream"),
            parent = EnvParents,
            body = joinpath(@__DIR__, "group_downstream.jl"),
        ),
        "FMI" => (;
            env = joinpath(@__DIR__, "fmi"),
            parent = EnvParents,
            body = joinpath(@__DIR__, "group_fmi.jl"),
        ),
        "Extensions" => (;
            env = joinpath(MTKBasePath, "test", "extensions"),
            parent = EnvParents,
            body = joinpath(@__DIR__, "group_extensions.jl"),
        ),
        "Optimization" => (;
            env = joinpath(MTKBasePath, "test", "optimization"),
            parent = EnvParents,
            body = joinpath(@__DIR__, "group_optimization.jl"),
        ),
    ),
    qa = (;
        env = joinpath(@__DIR__, "qa"),
        parent = EnvParents,
        body = joinpath(@__DIR__, "qa", "qa.jl"),
    ),
    all = [
        "InterfaceI",
        "Initialization",
        "InterfaceII",
        "SymbolicIndexingInterface",
        "Downstream",
        "FMI",
        "Extensions",
        "Optimization",
        "QA",
    ],
)
