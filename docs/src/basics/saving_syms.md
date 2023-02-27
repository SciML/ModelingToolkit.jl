# [Saving Only Certain Symbols](@id sym_save_idxs)

It is possible to specify symbolically which states of an ODESystem to save, with the `save_idxs` keyword argument
to the `solve` call. This may be important to you to save memory for large problems.

Take care to disable the `dense_output` flag when constructing the problem to ensure that observed
variables can be properly reconstructed. Failing to do this may cause incorrect construction of observed variables.

```julia
using ModelingToolkit, OrdinaryDiffEq, Test
@parameters t
@variables a(t) b(t) c(t) d(t) e(t)

D = Differential(t)

eqs = [D(a) ~ a,
    D(b) ~ b,
    D(c) ~ c,
    D(d) ~ d,
    e ~ d]

@named sys = ODESystem(eqs, t, [a, b, c, d, e], [];
                        defaults = Dict([a => 1.0,
                                            b => 1.0,
                                            c => 1.0,
                                            d => 1.0,
                                            e => 1.0]))
sys = structural_simplify(sys)
prob = ODEProblem(sys, [], (0, 1.0))
prob_sym = ODEProblem(sys, [], (0, 1.0), dense_output = false)

sol = solve(prob, Tsit5())
sol_sym = solve(prob_sym, Tsit5(), save_idxs = [a, c, e])

@test sol_sym[a] ≈ sol[a]
@test sol_sym[c] ≈ sol[c]
@test sol_sym[d] ≈ sol[d] # Dependency `d` of the observed variable `e` is automatically saved too.
@test sol_sym[e] ≈ sol[e]

@test sol.u != sol_sym.u

@test_throws Exception sol_sym[b]
```