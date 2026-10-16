# Dependency Graphs

# Dependency Graph API

The dependency graph constructors return a
[`BipartiteGraph`](https://docs.sciml.ai/BipartiteGraphs/dev/api/) from
[`BipartiteGraphs.jl`](https://github.com/SciML/BipartiteGraphs.jl). The graph type and
its primitive operations are documented and versioned by BipartiteGraphs; use that API
when working with the returned graph.

## Constructing Dependency Graphs

```@docs
ModelingToolkitBase.equation_dependencies
ModelingToolkitBase.asgraph
ModelingToolkitBase.variable_dependencies
ModelingToolkitBase.asdigraph
ModelingToolkitBase.eqeq_dependencies
ModelingToolkitBase.varvar_dependencies
```

## Variable-to-Equation Mapping

```@docs
ModelingToolkit.map_variables_to_equations
```
