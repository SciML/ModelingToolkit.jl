# Frequently Asked Questions

## Getting the index for a symbol

Since **ordering of symbols is not guaranteed after symbolic transformations**,
one should normally refer to values by their name. For example, `sol[lorenz.x]`
from the solution. But what if you need to get the index? The following helper
function will do the trick:

```julia
indexof(sym,syms) = findfirst(isequal(sym),syms)
indexof(σ,parameters(sys))
```

## Transforming value maps to arrays

ModelingToolkit.jl allows (and recommends) input maps like `[x => 2.0, y => 3.0]`
because symbol ordering is not guaranteed. However, what if you want to get the
lowered array? You can use the internal function `varmap_to_vars`. For example:

```julia
pnew = varmap_to_vars([β=>3.0, c=>10.0, γ=>2.0],parameters(sys))
```

## How do I handle `if` statements in my symbolic forms?

For statements that are in the `if then else` form, use `IfElse.ifelse` from the
[IfElse.jl](https://github.com/SciML/IfElse.jl) package to represent the code in a
functional form. For handling direct `if` statements, you can use equivalent boolean
mathematical expressions. For example `if x > 0 ...` can be implemented as just
`(x > 0) * `, where if `x <= 0` then the boolean will evaluate to `0` and thus the
term will be excluded from the model.

## ERROR: TypeError: non-boolean (Num) used in boolean context?

If you see the error:

```julia
ERROR: TypeError: non-boolean (Num) used in boolean context
```

then it's likely you are trying to trace through a function which cannot be
directly represented in Julia symbols. The techniques to handle this problem,
such as `@register_symbolic`, are described in detail
[in the Symbolics.jl documentation](https://symbolics.juliasymbolics.org/dev/manual/faq/#Transforming-my-function-to-a-symbolic-equation-has-failed.-What-do-I-do?-1).

## Using ModelingToolkit with Optimization / Automatic Differentiation

If you are using ModelingToolkit inside a loss function and are having issues with
mixing MTK with automatic differentiation, getting performance, etc... don't! Instead, use
MTK outside the loss function to generate the code, and then use the generated code
inside the loss function.

For example, let's say you were building ODEProblems in the loss function like:

```julia
function loss(p)
    prob = ODEProblem(sys, [], [p1 => p[1], p2 => p[2]])
    sol = solve(prob, Tsit5())
    sum(abs2,sol)
end
```

Since `ODEProblem` on a MTK `sys` will have to generate code, this will be slower than
caching the generated code, and will require automatic differentiation to go through the
code generation process itself. All of this is unnecessary. Instead, generate the problem
once outside the loss function, and remake the prob inside of the loss function:

```julia
prob = ODEProblem(sys, [], [p1 => p[1], p2 => p[2]])
function loss(p)
    remake(prob,p = ...)
    sol = solve(prob, Tsit5())
    sum(abs2,sol)
end
```

Now, one has to be careful with `remake` to ensure that the parameters are in the right
order. One can use the previously mentioned indexing functionality to generate index
maps for reordering `p` like:

```julia
p = @parameters x y z
idxs = ModelingToolkit.varmap_to_vars([p[1] => 1, p[2] => 2, p[3] => 3], p)
p[idxs]
```

Using this, the fixed index map can be used in the loss function. This would look like:

```julia
prob = ODEProblem(sys, [], [p1 => p[1], p2 => p[2]])
idxs = Int.(ModelingToolkit.varmap_to_vars([p1 => 1, p2 => 2], p))
function loss(p)
    remake(prob,p = p[idxs])
    sol = solve(prob, Tsit5())
    sum(abs2,sol)
end
```
