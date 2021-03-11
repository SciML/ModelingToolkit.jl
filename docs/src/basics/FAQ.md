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
