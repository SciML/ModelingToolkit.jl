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

edges, vars, vars_asso = sys2bigraph(pendulum)
@test ModelingToolkit.matching(edges, length(vars), vars_asso .== 0) == [0, 0, 0, 0, 1, 2, 3, 4, 0]

edges, assign, vars_asso, eqs_asso = ModelingToolkit.pantelides(pendulum)

@test edges == [
 [5, 3],               # 1
 [6, 4],               # 2
 [7, 9, 1],            # 3
 [8, 9, 2],            # 4
 [2, 1],               # 5
 [2, 1, 6, 5],         # 6
 [5, 3, 10, 7],        # 7
 [6, 4, 11, 8],        # 8
 [2, 1, 6, 5, 11, 10], # 9
]
#                  [1, 2, 3, 4, 5,  6,  7,  8,  9, 10,   11]
#                  [x, y, w, z, x', y', w', z', T, x'', y'']
@test vars_asso == [5, 6, 7, 8, 10, 11, 0,  0,  0,  0,    0]
#1: D(x) ~ w
#2: D(y) ~ z
#3: D(w) ~ T*x
#4: D(z) ~ T*y - g
#5: 0 ~ x^2 + y^2 - L^2
# ----
#6: D(5) -> 0 ~ 2xx'+ 2yy'
#7: D(1) -> D(D(x)) ~ D(w)
#8: D(2) -> D(D(y)) ~ D(z)
#9: D(6) -> 0 ~ 2xx'' + 2x'x' + 2yy'' + 2y'y'
#                 [1, 2, 3, 4, 5, 6, 7, 8, 9]
@test eqs_asso == [7, 8, 0, 0, 6, 9, 0, 0, 0]
