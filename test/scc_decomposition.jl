using ModelingToolkit, Symbolics, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables x(t) a(t) b(t) c(t)
@parameters p q
@named sys = System(
    [
        D(x) ~ -x + c,
        0 ~ a + 2b - p,
        0 ~ a - b + x,
        0 ~ c^3 + c - q * a,
    ], t
)
@test_throws ArgumentError ModelingToolkit.scc_decomposition(complete(sys))

sys = mtkcompile(sys)
dec = ModelingToolkit.scc_decomposition(sys)
dvs = unknowns(sys)
eqs = equations(sys)
eqvars(i) = Set(Symbolics.get_variables(eqs[only(dec.eq_sccs[i])]))

# The differential equation is excluded; the linear block in `a` precedes the
# nonlinear block in `c`, which depends on it.
@test length(dec.var_sccs) == length(dec.eq_sccs) == 2
@test isequal(dvs[only(dec.var_sccs[1])], a)
@test isequal(dvs[only(dec.var_sccs[2])], c)
@test a in eqvars(1) && !(c in eqvars(1))
@test c in eqvars(2) && a in eqvars(2)
@test dec.islinear == [true, false]
