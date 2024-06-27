# Working with Precompilation and Binary Building

## tl;dr, I just want precompilation to work

The tl;dr is, if you want to make precompilation work then instead of

```julia
ODEProblem(sys, u0, tspan, p)
```

use:

```julia
ODEProblem(sys, u0, tspan, p, eval_module = @__MODULE__, eval_expression = true)
```

As a full example, here's an example of a module that would precompile effectively:

```julia
module PrecompilationMWE
using ModelingToolkit

@variables x(ModelingToolkit.t_nounits)
@named sys = ODESystem([ModelingToolkit.D_nounits(x) ~ -x + 1], ModelingToolkit.t_nounits)
prob = ODEProblem(structural_simplify(sys), [x => 30.0], (0, 100), [],
    eval_expression = true, eval_module = @__MODULE__)

end
```

If you use that in your package's code then 99% of the time that's the right answer to get
precompilation working.

## I'm doing something fancier and need a bit more of an explanation

Oh you dapper soul, time for the bigger explanation. Julia's `eval` function evaluates a
function into a module at a specified world-age. If you evaluate a function within a function
and try to call it from within that same function, you will hit a world-age error. This looks like:

```julia
function worldageerror()
    f = eval(:((x) -> 2x))
    f(2)
end
```

```
julia> worldageerror()
ERROR: MethodError: no method matching (::var"#5#6")(::Int64)

Closest candidates are:
  (::var"#5#6")(::Any) (method too new to be called from this world context.)
   @ Main REPL[12]:2
```

This is done for many reasons, in particular if the code that is called within a function could change
at any time, then Julia functions could not ever properly optimize because the meaning of any function
or dispatch could always change and you would lose performance by guarding against that. For a full
discussion of world-age, see [this paper](https://arxiv.org/abs/2010.07516).

However, this would be greatly inhibiting to standard ModelingToolkit usage because then something as
simple as building an ODEProblem in a function and then using it would get a world age error:

```julia
function wouldworldage()
    prob = ODEProblem(sys, [], (0.0, 1.0))
    sol = solve(prob)
end
```

The reason is because `prob.f` would be constructed via `eval`, and thus `prob.f` could not be called
in the function, which means that no solve could ever work in the same function that generated the
problem. That does mean that:

```julia
function wouldworldage()
    prob = ODEProblem(sys, [], (0.0, 1.0))
end
sol = solve(prob)
```

is fine, or putting

```julia
prob = ODEProblem(sys, [], (0.0, 1.0))
sol = solve(prob)
```

at the top level of a module is perfectly fine too. They just cannot happen in the same function.

This would be a major limitation to ModelingToolkit, and thus we developed
[RuntimeGeneratedFunctions](https://github.com/SciML/RuntimeGeneratedFunctions.jl) to get around
this limitation. It will not be described beyond that, it is dark art and should not be investigated.
But it does the job. But that does mean that it plays... oddly with Julia's compilation.

There are ways to force RuntimeGeneratedFunctions to perform their evaluation and caching within
a given module, but that is not recommended because it does not play nicely with Julia v1.9's
introduction of package images for binary caching.

Thus when trying to make things work with precompilation, we recommend using `eval`. This is
done by simply adding `eval_expression=true` to the problem constructor. However, this is not
a silver bullet because the moment you start using eval, all potential world-age restrictions
apply, and thus it is recommended this is simply used for evaluating at the top level of modules
for the purpose of precompilation and ensuring binaries of your MTK functions are built correctly.

However, there is one caveat that `eval` in Julia works depending on the module that it is given.
If you have `MyPackage` that you are precompiling into, or say you are using `juliac` or PackageCompiler
or some other static ahead-of-time (AOT) Julia compiler, then you don't want to accidentally `eval`
that function to live in ModelingToolkit and instead want to make sure it is `eval`'d to live in `MyPackage`
(since otherwise it will not cache into the binary). ModelingToolkit cannot know that in advance, and thus
you have to pass in the module you wish for the functions to "live" in. This is done via the `eval_module`
argument.

Hence `ODEProblem(sys, u0, tspan, p, eval_module=@__MODULE__, eval_expression=true)` will work if you
are running this expression in the scope of the module you wish to be precompiling. However, if you are
attempting to AOT compile a different module, this means that `eval_module` needs to be appropriately
chosen. And, because `eval_expression=true`, all caveats of world-age apply.
