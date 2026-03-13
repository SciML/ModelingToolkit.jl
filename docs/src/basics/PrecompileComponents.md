# Building Precompilation-Friendly Components

How components are _defined_ has a significant impact on how well ModelingToolkit precompiles.
The core principle is **type stability**: the `System` constructor and the functions it calls
can only precompile specialized, fast code paths when the types of all arguments are fully
known at compile time. Type instability causes the compiler to fall back to dynamic dispatch,
generating generic `Any`-typed code that is slow and hard to compile.

## Build `vars` and `params` as `SymbolicT[]` vectors using `push!`

The macro forms `@variables` and `@parameters` return heterogeneous types. For example, a scalar
variable `x(t)` is a `Num`, while an array variable `y(t)[1:2]` is an `Arr{Num,1}`. Collecting
them with a literal array expression like `[x, y]` or `SymbolicT[x, y]` produces a `Vector{Any}`
or triggers a specialized `getindex` method for each unique combination of argument types.

Instead, declare the collection as `SymbolicT[]` and add elements one at a time with `push!`.
`SymbolicT` is a public API type alias from Symbolics.jl that wraps both `Num` and `Arr` uniformly:

```julia
using Symbolics: SymbolicT

@variables begin
    x(t_nounits)
    y(t_nounits)
end
vars = SymbolicT[]
push!(vars, x)
push!(vars, y)

@parameters begin
    α = α
    β = β
end
params = SymbolicT[]
push!(params, α)
push!(params, β)
```

This guarantees that `vars` and `params` are `Vector{SymbolicT}` — a concrete, fully-typed
collection — so the `System` constructor receives an argument of a known type and can compile a
specialized, cacheable method for it.

## Build `initial_conditions` and `guesses` as `Dict{SymbolicT, SymbolicT}`

Passing initial conditions or guesses inline as keyword arguments (e.g.
`initial_conditions = [x => 1.0, y => ones(2)]`) creates a
`Vector{Any}`, and causes type-instabilities in the `System` constructor.

Instead, build an explicit `Dict{SymbolicT, SymbolicT}` and `push!` pairs into it:

```julia
initial_conditions = Dict{SymbolicT, SymbolicT}()
push!(initial_conditions, x => 3.1)
push!(initial_conditions, y => 1)

guesses = Dict{SymbolicT, SymbolicT}()
# Similar
```

An equivalent formulation is:

```julia
initial_conditions = Dict{SymbolicT, SymbolicT}()
initial_conditions[x] = 3.1
initial_conditions[y] = 1
```

This avoids the problematic `Dict{SymbolicT, SymbolicT}(pairs...)` constructor, which loops
over the `pairs` tuple and can cause issues when each element in `pairs` is a different type.

## Build the equations vector as `Equation[]`

Equations should be collected into an `Equation[]` vector and populated with `push!`:

```julia
eqs = Equation[]
push!(eqs, D_nounits(x) ~ α * x - β * x * y)
push!(eqs, D_nounits(y) ~ -δ * y + γ * x * y)
```

The main intention here is to ensure that `eqs` is a `Vector{Equation}`. The same may be
achieved using the array literal syntax:

```julia
eqs = Equation[
    D_nounits(x) ~ α * x - β * x * y,
    D_nounits(y) ~ -δ * y + γ * x * y,
]
```

Though this does compile a unique method for `Base.vect` (or similar methods) depending on
the number of values in the array literal.

## Ensure the list of subsystems is a `Vector{System}`

Pass the subsystems list as an explicit `System[]` rather than omitting it or passing an empty
array (`[]`). This removes ambiguity about the element type:

```julia
return System(eqs, t_nounits, vars, params;
    systems = System[], initial_conditions, guesses, name)
```

## Complete example

Putting it all together, here is a component definition that follows all of the guidelines above
and is designed to hit ModelingToolkit's well-compiled code paths:

```julia
using ModelingToolkit
using ModelingToolkit: SymbolicT

@component function LotkaVolterra(;
        name, α = 1.3, β = 0.9, γ = 0.8, δ = 1.8)
    @parameters begin
        α = α
        β = β
        γ = γ
        δ = δ
    end
    params = SymbolicT[]
    push!(params, α)
    push!(params, β)
    push!(params, γ)
    push!(params, δ)

    @variables begin
        x(t_nounits)
        y(t_nounits)
    end
    vars = SymbolicT[]
    push!(vars, x)
    push!(vars, y)

    initial_conditions = Dict{SymbolicT, SymbolicT}()
    push!(initial_conditions, x => 3.1)
    push!(initial_conditions, y => 1.5)

    guesses = Dict{SymbolicT, SymbolicT}()

    eqs = Equation[]
    push!(eqs, D_nounits(x) ~ α * x - β * x * y)
    push!(eqs, D_nounits(y) ~ -δ * y + γ * x * y)

    return System(eqs, t_nounits, vars, params;
        systems = System[], initial_conditions, guesses, name)
end
```

## Use PrecompileTools.jl

Exercising all/many of the component constructors will help ensure they are precompiled
appropriately.

## Diagnosing type instabilities

If you want to verify that your component definition is type-stable, [Cthulhu.jl](https://github.com/JuliaDebug/Cthulhu.jl)
is the most reliable tool. Use `@descend` on the component constructor and look for red
(runtime-dispatch) calls. In Cthulhu, press `T` to switch to Typed IR, `o` to enable
optimizations, and `w` to highlight unstable code — this view is harder to read but gives the
most accurate picture of what the compiler actually infers. The Julia flags `--trace-compile`
and `--trace-compile-timing` are also useful for identifying methods that were not 
precompiled or were invalidated at load time.

!!! note "Julia 1.12 and precompilation"
    The `System` constructor, `complete`, and `mtkcompile` are specifically optimized to precompile
    as much as possible on Julia 1.12 and later. While the improvements also benefit earlier Julia
    versions, 1.12's compiler infrastructure (in particular, changes to how inference handles mutual
    recursion) allows significantly more of the compiled code to be cached. If precompilation
    performance is a priority, testing on Julia 1.12 will give the most informative results.

!!! note
    The discussion in [ModelingToolkit#4270](https://github.com/SciML/ModelingToolkit.jl/issues/4270)
    can also be of use here.
