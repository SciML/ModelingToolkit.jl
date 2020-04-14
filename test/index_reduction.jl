using ModelingToolkit
using ModelingToolkit: sys2bigraph
using DiffEqBase
using Test

# Define some variables
@parameters t L g
@variables x(t) y(t) w(t) z(t) T(t) x_t(t) y_t(t)
@derivatives D'~t

eqs2 = [D(D(x)) ~ T*x,
        D(D(y)) ~ T*y - g,
        0 ~ x^2 + y^2 - L^2]
pendulum2 = ODESystem(eqs2, t, [x, y, T], [L, g], name=:pendulum)
lowered_sys = ModelingToolkit.ode_order_lowering(pendulum2)

lowered_eqs = [D(x_t) ~ T*x,
               D(y_t) ~ T*y - g,
               0 ~ x^2 + y^2 - L^2,
               D(x) ~ x_t,
               D(y) ~ x_t]
@test_skip ODESystem(lowered_eqs) == lowered_sys # not gonna work
@test_broken isequal(lowered_sys.eqs, lowered_eqs)

# Simple pendulum in cartesian coordinates
eqs = [D(x) ~ w,
       D(y) ~ z,
       D(w) ~ T*x,
       D(z) ~ T*y - g,
       0 ~ x^2 + y^2 - L^2]
pendulum = ODESystem(eqs, t, [x, y, w, z, T], [L, g], name=:pendulum)

# V-nodes D(x), D(y), D(w), D(z), T
# E-nodes
vars, edges = sys2bigraph(pendulum)

@test ModelingToolkit.matching(edges, length(vars)) == [1, 2, 3, 4, 0]
