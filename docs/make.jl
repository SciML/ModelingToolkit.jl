using Documenter, ModelingToolkit
using ModelingToolkit: SciMLBase
# To load docstring from extension
import FMI, CommonSolve, JumpProcesses

MTKFMIExt = Base.get_extension(ModelingToolkit, :MTKFMIExt)

# Make sure that plots don't throw a bunch of warnings / errors!
ENV["GKSwstype"] = "100"
using Plots

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

mathengine = MathJax3(Dict(:loader => Dict("load" => ["[tex]/require", "[tex]/mathtools"]),
    :tex => Dict("inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
        "packages" => [
            "base",
            "ams",
            "autoload",
            "mathtools",
            "require"
        ])))

makedocs(sitename = "ModelingToolkit.jl",
    authors = "Chris Rackauckas",
    modules = [ModelingToolkit, MTKFMIExt],
    clean = true, doctest = false, linkcheck = true,
    warnonly = [:docs_block, :missing_docs, :cross_references],
    linkcheck_ignore = [
        "https://epubs.siam.org/doi/10.1137/0903023",
        # this link tends to fail linkcheck stochastically and often takes much longer to succeed
        # even in the browser it takes ages
        "http://www.scholarpedia.org/article/Differential-algebraic_equations"
    ],
    format = Documenter.HTML(;
        assets = ["assets/favicon.ico"],
        mathengine,
        canonical = "https://docs.sciml.ai/ModelingToolkit/stable/",
        prettyurls = (get(ENV, "CI", nothing) == "true"),
        # This page gets especially big with all the problem docstrings
        size_threshold_ignore = ["API/problems.md"]),
    pages = pages)

deploydocs(repo = "github.com/SciML/ModelingToolkit.jl.git";
    push_preview = true)
