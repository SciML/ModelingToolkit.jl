using ModelingToolkit, Test
@parameters α β δ
expr = (((1 / β - 1) + δ) / α) ^ (1 / (α - 1))
ref = [β, δ, α]
sol = Num.(ModelingToolkit.get_variables(expr))
@test all(simplify, sol[i] == ref[i] for i in 1:3)

@parameters γ
s = α => γ
expr = (((1 / β - 1) + δ) / α) ^ (1 / (α - 1))
sol = ModelingToolkit.substitute(expr, s)
new = (((1 / β - 1) + δ) / γ) ^ (1 / (γ - 1))
@test iszero(sol - new)

# test namespace_expr
@parameters t a p(t)
pterm = p.val
pnsp = ModelingToolkit.namespace_expr(pterm, :namespace, :t)
@test typeof(pterm) == typeof(pnsp)
@test ModelingToolkit.getname(pnsp) == Symbol("namespace₊p")
asym = a.val
ansp = ModelingToolkit.namespace_expr(asym, :namespace, :t)
@test typeof(asym) == typeof(ansp)
@test ModelingToolkit.getname(ansp) == Symbol("namespace₊a")
