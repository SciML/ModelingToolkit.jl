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
mathematical expressions. For example `if x > 0 ...` can be implementated as just
`(x > 0) * `, where if `x <= 0` then the boolean will evaluate to `0` and thus the
term will be excluded from the model.

## ERROR: TypeError: non-boolean (Num) used in boolean context?

If you see the error:

```julia
ERROR: TypeError: non-boolean (Num) used in boolean context
```

then it's likely you are trying to trace through a function which cannot be
directly represented in Julia symbols. The techniques to handle this problem,
such as `@register`, are described in detail 
[in the Symbolics.jl documentation](https://symbolics.juliasymbolics.org/dev/manual/faq/#Transforming-my-function-to-a-symbolic-equation-has-failed.-What-do-I-do?-1).
