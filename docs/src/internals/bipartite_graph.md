# Bipartite Graphs

!!! warning "Internal API"
    The functions documented on this page are internal implementation details of ModelingToolkit. They are not part of the public API and may change or be removed without notice in non-breaking releases. This documentation is provided to help contributors understand the codebase.

ModelingToolkit uses bipartite graphs to represent relationships between equations and
variables in systems. This page documents ModelingToolkit's developer-facing dependency
graph helpers. Extensions that need the underlying graph representation and primitive
operations should depend on `BipartiteGraphs.jl` directly rather than on a ModelingToolkit
reexport.

## Underlying Graph API

The following primitive graph API is defined and versioned by the
[`BipartiteGraphs.jl` API reference](https://docs.sciml.ai/BipartiteGraphs/dev/api/):

- Types: `BipartiteGraphs.BipartiteGraph`, `BipartiteGraphs.BipartiteEdge`, and
  `BipartiteGraphs.DiCMOBiGraph`.
- Vertex operations: `BipartiteGraphs.𝑠vertices`, `BipartiteGraphs.𝑑vertices`,
  `BipartiteGraphs.has_𝑠vertex`, `BipartiteGraphs.has_𝑑vertex`,
  `BipartiteGraphs.nsrcs`, and `BipartiteGraphs.ndsts`.
- Neighbor and edge operations: `BipartiteGraphs.𝑠neighbors`,
  `BipartiteGraphs.𝑑neighbors`, `BipartiteGraphs.set_neighbors!`,
  `BipartiteGraphs.𝑠edges`, and `BipartiteGraphs.𝑑edges`.
- Views and modifications: `BipartiteGraphs.invview`,
  `BipartiteGraphs.delete_srcs!`, and `BipartiteGraphs.delete_dsts!`.
- Matching: `BipartiteGraphs.maximal_matching` and
  `BipartiteGraphs.construct_augmenting_path!`.
- Vertex kinds: `BipartiteGraphs.SRC` and `BipartiteGraphs.DST`.
