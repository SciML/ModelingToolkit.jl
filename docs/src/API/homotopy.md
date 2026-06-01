```@meta
CollapsedDocStrings = true
```

# Homotopy continuation for initialization (Modelica `homotopy`)

ModelingToolkit implements Modelica's `homotopy(actual, simplified)` operator
(Modelica specification §3.7.4). It is an initialization aid: where the `actual`
equations are hard to solve from the available guesses, the operator lets
initialization start from an easy `simplified` problem and continuously deform it
into the `actual` one.

During `complete`/`mtkcompile`, every `homotopy(actual, simplified)` node is lowered
to

```math
(1 - \lambda) \cdot \mathrm{simplified} + \lambda \cdot \mathrm{actual}
```

with a single shared parameter `__homotopy_λ` (default `1.0`). At `λ = 1` the system
reduces to `actual`; at `λ = 0` to `simplified`. For systems that contain `homotopy`
nodes, the default initialization algorithm becomes
`OverrideInit(nlsolve = TrivialThenSweep(...))`: it first attempts the trivial
single solve at `λ = 1` and, if that fails, runs a parameter-sweep continuation that
walks `λ` from 0 to 1 (mirroring OpenModelica's default user experience). Pass an
explicit `initializealg` to the problem constructor or to `solve` to override this.

At runtime, outside initialization, `homotopy(actual, simplified)` evaluates to
`actual`, as required by the specification.

!!! note
    This operator is unrelated to the polynomial `HomotopyContinuationProblem`, which
    solves polynomial systems via homotopy continuation. They only share the word
    "homotopy".

## Operator

```@docs
ModelingToolkit.homotopy
```

## Initialization algorithms

```@docs
ModelingToolkit.TrivialThenSweep
ModelingToolkit.HomotopySweep
ModelingToolkit.TrivialHomotopy
```
