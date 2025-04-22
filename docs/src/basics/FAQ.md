# Frequently Asked Questions

## Why are my parameters some obscure object?

In ModelingToolkit.jl version 9, the parameter vector was replaced with a custom
`MTKParameters` object, whose internals are intentionally undocumented and subject
to change without a breaking release. This enables us to efficiently store and generate
code for parameters of multiple types. To obtain parameter values use
[SymbolicIndexingInterface.jl](https://github.com/SciML/SymbolicIndexingInterface.jl/) or
[SciMLStructures.jl](https://github.com/SciML/SciMLStructures.jl/). For example:

```julia
prob.ps[lorenz.β] # obtains the value of parameter `β`. Note the `.ps` instead of `.p`
getβ = getp(prob, lorenz.β) # returns a function that can fetch the value of `β`
getβ(sol) # can be used on any object that is based off of the same system
getβ(prob)
```

Indexes into the `MTKParameters` object take the form of `ParameterIndex` objects, which
are similarly undocumented. Following is the list of behaviors that should be relied on for
`MTKParameters`:

  - It implements the SciMLStructures interface.
  - It can be queried for parameters using functions returned from
    `SymbolicIndexingInterface.getp`.
  - `getindex(::MTKParameters, ::ParameterIndex)` can be used to obtain the value of a
    parameter with the given index.
  - `setindex!(::MTKParameters, value, ::ParameterIndex)` can be used to set the value of a
    parameter with the given index.
  - `parameter_values(sys, sym)` will return a `ParameterIndex` object if `sys` has been
    `complete`d (through `structural_simplify`, `complete` or `@mtkbuild`).
  - `copy(::MTKParameters)` is defined and duplicates the parameter object, including the
    memory used by the underlying buffers.

Any other behavior of `MTKParameters` (other `getindex`/`setindex!` methods, etc.) is an
undocumented internal and should not be relied upon.

## How do I use non-numeric/array-valued parameters?

In ModelingToolkit.jl version 9, parameters are required to have a `symtype` matching
the type of their values. For example, this will error during problem construction:

```julia
@parameters p = [1, 2, 3]
```

Since by default parameters have a `symtype` of `Real` (which is interpreted as `Float64`)
but the default value given to it is a `Vector{Int}`. For array-valued parameters, use the
following syntax:

```julia
@parameters p[1:n, 1:m]::T # `T` is the `eltype` of the parameter array
@parameters p::T # `T` is the type of the array
```

The former approach is preferred, since the size of the array is known. If the array is not
a `Base.Array` or the size is not known during model construction, the second syntax is
required.

The same principle applies to any parameter type that is not `Float64`.

```julia
@parameters p1::Int # integer-valued
@parameters p2::Bool # boolean-valued
@parameters p3::MyCustomStructType # non-numeric
@parameters p4::ComponentArray{...} # non-standard array
```

## Getting the index for a symbol

Ordering of symbols is not guaranteed after symbolic transformations, and parameters
are now stored in a custom `MTKParameters` object instead of a vector. Thus, values
should be referred to by their name. For example `sol[lorenz.x]`. To obtain the index,
use the following functions from
[SymbolicIndexingInterface.jl](https://github.com/SciML/SymbolicIndexingInterface.jl/):

```julia
variable_index(sys, sym)
parameter_index(sys, sym)
```

Note that while the variable index will be an integer, the parameter index is a struct of
type `ParameterIndex` whose internals should not be relied upon.

## Can I index with strings?

Strings are not considered symbolic variables, and thus cannot directly be used for symbolic
indexing. However, ModelingToolkit does provide a method to parse the string representation of
a variable, given the system in which that variable exists.

```@docs
ModelingToolkit.parse_variable
```

## Transforming value maps to arrays

ModelingToolkit.jl allows (and recommends) input maps like `[x => 2.0, y => 3.0]`
because symbol ordering is not guaranteed. However, what if you want to get the
lowered array? You can use the internal function `varmap_to_vars` for unknowns.
and the `MTKParameters` constructor for parameters. For example:

```julia
unew = varmap_to_vars([x => 1.0, y => 2.0, z => 3.0], unknowns(sys))
pnew = ModelingToolkit.MTKParameters(sys, [β => 3.0, c => 10.0, γ => 2.0], unew)
```

## How do I handle `if` statements in my symbolic forms?

For statements that are in the `if then else` form, use `Base.ifelse` from the
to represent the code in a functional form. For handling direct `if` statements,
you can use equivalent boolean mathematical expressions. For example, `if x > 0 ...`
can be implemented as just `(x > 0) * `, where if `x <= 0` then the boolean will
evaluate to `0` and thus the term will be excluded from the model.

## ERROR: TypeError: non-boolean (Num) used in boolean context?

If you see the error:

```
ERROR: TypeError: non-boolean (Num) used in boolean context
```

then it's likely you are trying to trace through a function which cannot be
directly represented in Julia symbols. The techniques to handle this problem,
such as `@register_symbolic`, are described in detail
[in the Symbolics.jl documentation](https://symbolics.juliasymbolics.org/dev/manual/faq/#Transforming-my-function-to-a-symbolic-equation-has-failed.-What-do-I-do?-1).

## Using ModelingToolkit with Optimization / Automatic Differentiation

If you are using ModelingToolkit inside a loss function and are having issues with
mixing MTK with automatic differentiation, getting performance, etc… don't! Instead, use
MTK outside the loss function to generate the code, and then use the generated code
inside the loss function.

For example, let's say you were building ODEProblems in the loss function like:

```julia
function loss(p)
    prob = ODEProblem(sys, [], [p1 => p[1], p2 => p[2]])
    sol = solve(prob, Tsit5())
    sum(abs2, sol)
end
```

Since `ODEProblem` on a MTK `sys` will have to generate code, this will be slower than
caching the generated code, and will require automatic differentiation to go through the
code generation process itself. All of this is unnecessary. Instead, generate the problem
once outside the loss function, and update the parameter values inside the loss function:

```julia
prob = ODEProblem(sys, [], [p1 => p[1], p2 => p[2]])
function loss(p)
    # update parameters
    sol = solve(prob, Tsit5())
    sum(abs2, sol)
end
```

If a subset of the parameters are optimized, `setp` from SymbolicIndexingInterface.jl
should be used to generate an efficient function for setting parameter values. For example:

```julia
using SymbolicIndexingInterface

prob = ODEProblem(sys, [], [p1 => p[1], p2 => p[2]])
setter! = setp(sys, [p1, p2])
function loss(p)
    setter!(prob, p)
    sol = solve(prob, Tsit5())
    sum(abs2, sol)
end
```

[SciMLStructures.jl](https://github.com/SciML/SciMLStructures.jl/) can be leveraged to
obtain all the parameters for optimization using the `Tunable` portion. By default, all
numeric or numeric array parameters are marked as tunable, unless explicitly marked as
`tunable = false` in the variable metadata.

```julia
using SciMLStructures: replace!, Tunable

prob = ODEProblem(sys, [], [p1 => p[1], p2 => p[2]])
function loss(p)
    replace!(Tunable(), prob.p, p)
    sol = solve(prob, Tsit5())
    sum(abs2, sol)
end

p, replace, alias = SciMLStructures.canonicalize(Tunable(), prob.p)
# p is an `AbstractVector` which can be optimized
# if `alias == true`, then `p` aliases the memory used by `prob.p`, so
# changes to the array will be reflected in parameter values
```

# ERROR: ArgumentError: SymbolicUtils.BasicSymbolic{Real}[xˍt(t)] are missing from the variable map.

This error can come up after running `structural_simplify` on a system that generates dummy derivatives (i.e. variables with `ˍt`).  For example, here even though all the variables are defined with initial values, the `ODEProblem` generation will throw an error that defaults are missing from the variable map.

```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

sts = @variables x1(t)=0.0 x2(t)=0.0 x3(t)=0.0 x4(t)=0.0
eqs = [x1 + x2 + 1 ~ 0
       x1 + x2 + x3 + 2 ~ 0
       x1 + D(x3) + x4 + 3 ~ 0
       2 * D(D(x1)) + D(D(x2)) + D(D(x3)) + D(x4) + 4 ~ 0]
@named sys = ODESystem(eqs, t)
sys = structural_simplify(sys)
prob = ODEProblem(sys, [], (0, 1))
```

We can solve this problem by using the `missing_variable_defaults()` function

```julia
prob = ODEProblem(sys, ModelingToolkit.missing_variable_defaults(sys), (0, 1))
```

This function provides 0 for the default values, which is a safe assumption for dummy derivatives of most models.  However, the 2nd argument allows for a different default value or values to be used if needed.

```
julia> ModelingToolkit.missing_variable_defaults(sys, [1,2,3])
3-element Vector{Pair}:
  x1ˍt(t) => 1
 x2ˍtt(t) => 2
 x3ˍtt(t) => 3
```

## Change the unknown variable vector type

Use the `u0_constructor` keyword argument to map an array to the desired
container type. For example:

```julia
using ModelingToolkit, StaticArrays
using ModelingToolkit: t_nounits as t, D_nounits as D

sts = @variables x1(t) = 0.0
eqs = [D(x1) ~ 1.1 * x1]
@mtkbuild sys = ODESystem(eqs, t)
prob = ODEProblem{false}(sys, [], (0, 1); u0_constructor = x -> SVector(x...))
```

## Using a custom independent variable

When possible, we recommend `using ModelingToolkit: t_nounits as t, D_nounits as D` as the independent variable and its derivative.
However, if you want to use your own, you can do so:

```julia
using ModelingToolkit

@independent_variables x
D = Differential(x)
@variables y(x)
@named sys = ODESystem([D(y) ~ x], x)
```

## Ordering of tunable parameters

Tunable parameters are floating point parameters, not used in callbacks and not marked with `tunable = false` in their metadata. These are expected to be used with AD
and optimization libraries. As such, they are stored together in one `Vector{T}`. To obtain the ordering of tunable parameters in this buffer, use:

```@docs
tunable_parameters
```

If you have an array in which a particular dimension is in the order of tunable parameters (e.g. the jacobian with respect to tunables) then that dimension of the
array can be reordered into the required permutation using the symbolic variables:

```@docs
reorder_dimension_by_tunables!
reorder_dimension_by_tunables
```

For example:

```@example reorder
using ModelingToolkit

@parameters p q[1:3] r[1:2, 1:2]

@named sys = ODESystem(Equation[], ModelingToolkit.t_nounits, [], [p, q, r])
sys = complete(sys)
```

The canonicalized tunables portion of `MTKParameters` will be in the order of tunables:

```@example reorder
using SciMLStructures: canonicalize, Tunable

ps = MTKParameters(sys, [p => 1.0, q => [2.0, 3.0, 4.0], r => [5.0 6.0; 7.0 8.0]])
arr = canonicalize(Tunable(), ps)[1]
```

We can reorder this to contain the value for `p`, then all values for `q`, then for `r` using:

```@example reorder
reorder_dimension_by_tunables(sys, arr, [p, q, r])
```

This also works with interleaved subarrays of symbolics:

```@example reorder
reorder_dimension_by_tunables(sys, arr, [q[1], r[1, :], q[2], r[2, :], q[3], p])
```

And arbitrary dimensions of higher dimensional arrays:

```@example reorder
highdimarr = stack([i * arr for i in 1:5]; dims = 1)
```

```@example reorder
reorder_dimension_by_tunables(sys, highdimarr, [q[1:2], r[1, :], q[3], r[2, :], p]; dim = 2)
```
