# Parsing Expressions into Solvable Systems

Many times when creating DSLs or creating ModelingToolkit extensions to read new file formats,
it can become imperative to parse expressions. In many cases, it can be easy to use `Base.parse`
to take things to standard Julia expressions, but how can you take a `Base.Expr` and generate
symbolic forms from that? For example, say we had the following system we wanted to solve:

```@example parsing
ex = [:(y ~ x)
      :(y ~ -2x + 3 / z)
      :(z ~ 2)]
```

We can use the function `parse_expr_to_symbolic` from Symbolics.jl to generate the symbolic
form of the expression:

```@example parsing
using Symbolics
eqs = parse_expr_to_symbolic.(ex, (Main,))
```

From there, we can use ModelingToolkit to transform the symbolic equations into a numerical
nonlinear solve:

```@example parsing
using ModelingToolkit, SymbolicIndexingInterface, NonlinearSolve
vars = union(ModelingToolkit.vars.(eqs)...)
@mtkbuild ns = NonlinearSystem(eqs, vars, [])

varmap = Dict(SymbolicIndexingInterface.getname.(vars) .=> vars)
prob = NonlinearProblem(ns, [varmap[:x] => 1.0, varmap[:y] => 1.0, varmap[:z] => 1.0])
sol = solve(prob, NewtonRaphson())
```
