# Bipartite Graphs

!!! warning "Internal API"
    The functions documented on this page are internal implementation details of ModelingToolkit. They are not part of the public API and may change or be removed without notice in non-breaking releases. This documentation is provided to help contributors understand the codebase.

ModelingToolkit uses bipartite graphs to represent relationships between equations and variables in systems. These functions provide tools for working with and analyzing these graphs.

## Graph Types

```@docs
BipartiteGraph
BipartiteEdge
DiCMOBiGraph
DiffGraph
```

## Vertex Operations

```@docs
𝑠vertices
𝑑vertices
has_𝑠vertex
has_𝑑vertex
nsrcs
ndsts
```

## Neighbor Operations

```@docs
𝑠neighbors
𝑑neighbors
set_neighbors!
```

## Edge Operations

```@docs
𝑠edges
𝑑edges
```

## Graph Views and Modifications

```@docs
invview
delete_srcs!
delete_dsts!
```

## Matching Algorithms

```@docs
maximal_matching
construct_augmenting_path!
```

## Dependency Analysis

```@docs
equation_dependencies
variable_dependencies
eqeq_dependencies
varvar_dependencies
map_variables_to_equations
```

## Graph Conversion

```@docs
asgraph
asdigraph
```

## Constants

```@docs
SRC
DST
```