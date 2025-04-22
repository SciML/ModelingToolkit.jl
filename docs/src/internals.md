# Internal Details

This is a page for detailing some of the inner workings to help future
contributors to the library.

## Observables and Variable Elimination

In the variable “elimination” algorithms, what is actually done is that variables
are removed from being unknowns and equations are moved into the `observed` category
of the system. The `observed` equations are explicit algebraic equations which
are then substituted out to completely eliminate these variables from the other
equations, allowing the system to act as though these variables no longer exist.

However, a user may want to interact with such variables, for example,
plotting their output. For this reason, these relationships are stored,
and are then used to generate the `observed` equation found in the
`SciMLFunction` interface, so that `sol[x]` lazily reconstructs the observed
variable when necessary. In this sense, there is an equivalence between
observables and the variable elimination system.

The procedure for variable elimination inside [`structural_simplify`](@ref) is

 1. [`ModelingToolkit.initialize_system_structure`](@ref).
 2. [`ModelingToolkit.alias_elimination`](@ref). This step moves equations into `observed(sys)`.
 3. [`ModelingToolkit.dae_index_lowering`](@ref) by means of [`pantelides!`](@ref) (if the system is an [`ODESystem`](@ref)).
 4. [`ModelingToolkit.tearing`](@ref).

## Preparing a system for simulation

Before a simulation or optimization can be performed, the symbolic equations stored in an [`AbstractSystem`](@ref) must be converted into executable code. This step typically occurs after the simplification explained above, and is performed when an instance of a [`SciMLBase.AbstractSciMLProblem`](@ref), such as a [`ODEProblem`](@ref), is constructed.
The call chain typically looks like this, with the function names in the case of an `ODESystem` indicated in parentheses

 1. Problem constructor ([`ODEProblem`](@ref))
 2. Build an `DEFunction` ([`process_DEProblem`](@ref) -> [`ODEFunction`](@ref)
 3. Write actual executable code ([`generate_function`](@ref) or [`generate_custom_function`](@ref))

Apart from [`generate_function`](@ref), which generates the dynamics function, `ODEFunction` also builds functions for observed equations (`build_explicit_observed_function`) and Jacobians (`generate_jacobian`) etc. These are all stored in the `ODEFunction`.

## Creating an `MTKParameters` object

It may be useful to create a parameter object without creating the problem. For this
purpose, the `MTKParameters` constructor is exposed as public API.

```@docs
MTKParameters
```
