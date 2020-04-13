using Documenter, ModelingToolkit

makedocs(
    sitename="ModelingToolkit.jl",
    modules=[ModelingToolkit],
    pages=[
        "Home" => "index.md",
        "api.md",
    ]
)
