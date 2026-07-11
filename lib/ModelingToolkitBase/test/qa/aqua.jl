using ModelingToolkitBase
using SciMLTesting

run_qa(
    ModelingToolkitBase;
    api_docs_kwargs = (;
        rendered = true,
        docs_src = joinpath(dirname(dirname(dirname(dirname(@__DIR__)))), "docs", "src"),
    ),
)
