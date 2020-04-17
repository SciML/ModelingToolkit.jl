using Documenter, ModelingToolkit

makedocs(
    sitename="ModelingToolkit.jl",
    authors="Chris Rackauckas",
    modules=[ModelingToolkit],
    clean=true,doctest=false,
    format = Documenter.HTML(#analytics = "UA-90474609-3",
                             assets = ["assets/favicon.ico"],
                             canonical="https://docs.sciml.ai/stable/"),
    pages=[
        "Home" => "index.md",
        "highlevel.md",
        "Systems" => Any[
            "systems/AbstractSystem.md",
            "systems/ODESystem.md",
            "systems/SDESystem.md",
            "systems/NonlinearSystem.md",
            "systems/OptimizationSystem.md",
            "systems/ReactionSystem.md",
            "systems/PDESystem.md"
        ],
        "IR.md"
    ]
)

deploydocs(
   repo = "github.com/SciML/ModelingToolkit.jl.git";
   push_preview = true
)
