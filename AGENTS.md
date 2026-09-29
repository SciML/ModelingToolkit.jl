# AGENTS.md

Repository-specific conventions for ModelingToolkit.jl. See `CONTRIBUTING.md` for the
style and formatting rules that apply to every SciML repository.

## Code generation

### Never `gensym` a name that ends up in generated code

Names embedded in an `Expr` that is handed to `eval_or_rgf` / `RuntimeGeneratedFunction`
must be fixed symbols, not `gensym`s. A `gensym` embeds a process-global counter, so the
same system lowers to a different `Expr` in the precompile process than in the user
session. That defeats the `RuntimeGeneratedFunctions` Expr-hash cache: precompiled bodies
are never hit, every `Expr` is a distinct type, and constructing the same problem twice
recompiles the generated function.

Use a fixed sentinel instead, prefixed `__mtk_` (see `generated_argument_name`) or
suffixed `ₘₜₖ` (see `HOMOTOPY_LAMBDA`, `__log_assertions_ₘₜₖ`). Both make a collision with
a user-chosen symbol implausible; where a collision would silently produce wrong code
rather than an error, guard against it explicitly as `lower_homotopy` does.

`gensym` remains fine for names that never reach generated code, such as the default
`name` of a system built by `modelingtoolkitize`.

### Build `Expr`s by pushing, not splatting

Prefer `push!`/`append!` onto `expr.args` over `Expr(head, xs...)`. The splat is a
dynamic call whose argument count is unknown to the compiler, so it inflates codegen time
and infers poorly:

```julia
# no
Expr(:tuple, buffers...)

# yes
tup = Expr(:tuple)
append!(tup.args, buffers)
```

### Keep codegen-side containers concretely typed vectors

Collections that codegen iterates over — parameter groups, buffer lists, argument lists —
should be a concretely typed `Vector`, not a tuple of tuples. Tuples force the whole
surrounding loop to specialize on the shape of each individual system, which is
catastrophic for inference and compile time. Reserve tuples for values that genuinely
need to be heterogeneous or statically sized in the *generated* code.
