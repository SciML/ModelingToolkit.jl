using ModelingToolkit
using ModelingToolkit: sys2bigraph
using DiffEqBase
using Test

# Define some variables
@parameters t L g
@variables x(t) y(t) w(t) z(t) T(t) xˍt(t) yˍt(t)
@derivatives D'~t

eqs2 = [D(D(x)) ~ T*x,
        D(D(y)) ~ T*y - g,
        0 ~ x^2 + y^2 - L^2]
pendulum2 = ODESystem(eqs2, t, [x, y, T], [L, g], name=:pendulum)
lowered_sys = ModelingToolkit.ode_order_lowering(pendulum2)

lowered_eqs = [D(xˍt) ~ T*x,
               D(yˍt) ~ T*y - g,
               0 ~ x^2 + y^2 - L^2,
               D(x) ~ xˍt,
               D(y) ~ yˍt]
@test ODESystem(lowered_eqs, t, [xˍt, yˍt, x, y, T], [L, g]) == lowered_sys
@test isequal(lowered_sys.eqs, lowered_eqs)

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

function new_system(sys, assign, vars_asso, eqs_asso)
end

@test sort.(edges) == sort.([
 [5, 3],               # 1
 [6, 4],               # 2
 [7, 9, 1],            # 3
 [8, 9, 2],            # 4
 [2, 1],               # 5
 [2, 1, 6, 5],         # 6
 [5, 3, 10, 7],        # 7
 [6, 4, 11, 8],        # 8
 [2, 1, 6, 5, 11, 10], # 9
])
#                  [1, 2, 3, 4, 5,  6,  7,  8,  9, 10,   11]
#                  [x, y, w, z, x', y', w', z', T, x'', y''] -- how can I get this vector of symbols?
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

using ModelingToolkit
@parameters t L g
@variables x(t) y(t) w(t) z(t) T(t) xˍt(t) yˍt(t)
@derivatives D'~t
idx1_pendulum = [D(x) ~ w,
                 D(y) ~ z,
                 D(w) ~ T*x,
                 D(z) ~ T*y - g,
                 #0 ~ x^2 + y^2 - L^2,
                 #0 ~ 2x*w + 2y*z,
                 #D(xˍt) ~ D(w),
                 D(xˍt) ~ T*x,
                 # D(D(x)) ~ D(w) and substitute the rhs
                 #D(D(x)) ~ T*x,
                 # D(D(y)) ~ D(z) and substitute the rhs
                 #D(D(y)) ~ T*y - g,
                 D(yˍt) ~ T*y - g,
                 # 2x*D(D(x)) + 2*D(x)*D(x) + 2y*D(D(y)) + 2*D(y)*D(y) and
                 # substitute the rhs
                 0 ~ 2x*(T*x) + 2*xˍt*xˍt + 2y*(T*y - g) + 2*yˍt*yˍt]
idx1_pendulum = ODESystem(idx1_pendulum, t, [x, y, w, z, xˍt, yˍt, T], [L, g])
first_order_idx1_pendulum = ode_order_lowering(idx1_pendulum)

using OrdinaryDiffEq
using LinearAlgebra
#             [x, y, w, z, xˍt, yˍt, T]
#M = Diagonal([1, 1, 1, 1,   1,   1, 0])
prob = ODEProblem(ODEFunction(first_order_idx1_pendulum),
        #  [x, y, w, z, xˍt, yˍt, T]
           [1, 0, 0, 0,   0,   0, 0.0],# 0, 0, 0, 0],
           (0, 100.0),
           [1, 9.8],
           mass_matrix=calculate_massmatrix(first_order_idx1_pendulum))
sol = solve(prob, Rodas5())
