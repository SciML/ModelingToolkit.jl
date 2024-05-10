using ModelingToolkit, Test
using ModelingToolkit: value
using SymbolicUtils: <ₑ
@parameters α β δ
expr = (((1 / β - 1) + δ) / α)^(1 / (α - 1))
ref = sort([β, δ, α], lt = <ₑ)
sol = sort(Num.(ModelingToolkit.get_variables(expr)), lt = <ₑ)
@test all(x -> x isa Num, sol[i] == ref[i] for i in 1:3)
@test all(simplify ∘ value, sol[i] == ref[i] for i in 1:3)

@parameters γ
s = α => γ
expr = (((1 / β - 1) + δ) / α)^(1 / (α - 1))
sol = ModelingToolkit.substitute(expr, s)
new = (((1 / β - 1) + δ) / γ)^(1 / (γ - 1))
@test iszero(sol - new)

# Continuous
using ModelingToolkit: isdifferential, vars, collect_differential_variables,
                       collect_ivs
@variables t u(t) y(t)
using ModelingToolkit: D
eq = D(y) ~ u
v = vars(eq)
@test v == Set([D(y), u])

ov = collect_differential_variables(eq)
@test ov == Set(Any[y])

aov = ModelingToolkit.collect_applied_operators(eq, Differential)
@test aov == Set(Any[D(y)])

ts = collect_ivs([eq])
@test ts == Set([t])

# Test units of diff vars of `Term` type.
using ModelingToolkit: t, get_unit, default_toterm
using DynamicQuantities
@variables k(t) [unit = u"kg"]
k2 = value(D(D(k)))
@test get_unit(default_toterm(k2)) == get_unit(k2)

@variables l(t) [unit = u"mg"]
l2 = value(D(D(l)))
@test_logs (:warn,
    """Ignoring the unit while converting `Differential(t)(Differential(t)(l(t)))` to a term.
1.0e-6 kg uses non SI unit. Please use SI unit only.""") default_toterm(value(l2))
