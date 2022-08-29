using Documenter, ModelingToolkit

include("pages.jl")

mathengine = MathJax3(Dict(:loader => Dict("load" => ["[tex]/require", "[tex]/mathtools"]),
                           :tex => Dict("inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                                        "packages" => [
                                            "base",
                                            "ams",
                                            "autoload",
                                            "mathtools",
                                            "require",
                                        ])))

makedocs(sitename = "ModelingToolkit.jl",
         authors = "Chris Rackauckas",
         modules = [ModelingToolkit],
         clean = true, doctest = false,
         strict = [
             :doctest,
             :linkcheck,
             :parse_error,
             :example_block,
             # Other available options are
             # :autodocs_block, :cross_references, :docs_block, :eval_block, :example_block, :footnote, :meta_block, :missing_docs, :setup_block
         ],
         format = Documenter.HTML(; analytics = "UA-90474609-3",
                                  assets = ["assets/favicon.ico"],
                                  mathengine,
                                  canonical = "https://mtk.sciml.ai/stable/",
                                  prettyurls = (get(ENV, "CI", nothing) == "true")),
         pages = pages)

deploydocs(repo = "github.com/SciML/ModelingToolkit.jl.git";
           push_preview = true)
