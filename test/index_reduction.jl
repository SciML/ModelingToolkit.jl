using ModelingToolkit
using ModelingToolkit: sys2bigraph
using DiffEqBase
using Test

# Define some variables
@parameters t L g
@variables x(t) y(t) w(t) z(t) T(t)
@derivatives D'~t

eqs2 = [D(D(x)) ~ T*x,
        D(D(y)) ~ T*y - g,
        0 ~ x^2 + y^2 - L^2]
pendulum2 = ODESystem(eqs, t, [x, y, T], [L, g], name=:pendulum)
ModelingToolkit.ode_order_lowering(pendulum2)

# Simple pendulum in cartesian coordinates
eqs = [D(x) ~ w,
       D(y) ~ z,
       D(w) ~ T*x,
       D(z) ~ T*y - g,
       0 ~ x^2 + y^2 - L^2]
pendulum = ODESystem(eqs, t, [x, y, w, z, T], [L, g], name=:pendulum)

sys2bigraph(pendulum)
