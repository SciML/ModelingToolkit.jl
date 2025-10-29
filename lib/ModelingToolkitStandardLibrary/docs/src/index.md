# ModelingToolkitStandardLibrary.jl

ModelingToolkitStandardLibrary.jl is a standard library for the
[ModelingToolkit](https://docs.sciml.ai/ModelingToolkit/stable/) acausal modeling system.

## Installation

To install ModelingToolkitStandardLibrary.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("ModelingToolkitStandardLibrary")
```

## Tutorials

  - [RC Circuit](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/rc_circuit/)
  - [Custom Component](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/custom_component/)
  - [Thermal Model](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/thermal_model/)
  - [DC Motor with PI-controller](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/dc_motor_pi/)

## Libraries

The following are the constituent libraries of the ModelingToolkit Standard Library.

  - [Basic Blocks](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/blocks/)
  - [Mechanical Components](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/mechanical/)
  - [Electrical Components](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/electrical/)
  - [Magnetic Components](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/magnetic/)
  - [Thermal Components](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/thermal/)
  - [Hydraulic Components](@ref hydraulic)

## Contributing

  - Please refer to the
    [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
    for guidance on PRs, issues, and other matters relating to contributing to SciML.

  - See the [SciML Style Guide](https://github.com/SciML/SciMLStyle) for common coding practices and other style decisions.
  - There are a few community forums:
    
      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Slack](https://julialang.org/slack/)
      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Zulip](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
      + On the [Julia Discourse forums](https://discourse.julialang.org)
      + See also [SciML Community page](https://sciml.ai/community/)

## Reproducibility

```@raw html
<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>
```

```@example
using Pkg # hide
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>and using this machine and Julia version.</summary>
```

```@example
using InteractiveUtils # hide
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```

```@example
using Pkg # hide
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```

```@eval
using TOML
using Markdown
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link_manifest = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
                "/assets/Manifest.toml"
link_project = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
               "/assets/Project.toml"
Markdown.parse("""You can also download the
[manifest]($link_manifest)
file and the
[project]($link_project)
file.
""")
```
