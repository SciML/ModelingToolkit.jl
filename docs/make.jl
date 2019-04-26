using Documenter, ModelingToolkit, DiffEqBase


makedocs(
    sitename="ModelingToolkit.jl",
    modules=[ModelingToolkit],
    pages=[
        "Home" => "index.md",
        "api.md",
    ]
)
