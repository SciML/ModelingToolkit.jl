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

## Embedding data into a symbolic model

Let's say for example you want to embed data for the timeseries of some
forcing equations into the right-hand side of and ODE, or data into a PDE. What
you would do in these cases is use the `@register` function over an interpolator.
For example, [DataInterpolations.jl](https://github.com/PumasAI/DataInterpolations.jl)
is a good library for interpolating the data. Then you can do:

```julia
spline = CubicSpline(data,datat)
f(t) = spline(t)
@register f(t)
```

This will make `f(t)` be a function that Symbolics.jl will not attempt to trace.
One should also consider defining the derivative to the function, if available.
