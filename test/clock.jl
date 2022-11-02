using ModelingToolkit, Test
dt = 0.1
@variables t x(t) y(t) u(t) yd(t) ud(t) r(t)
@parameters kp
D = Differential(t)

eqs = [yd ~ Sample(t, dt)(y)
       ud ~ kp * (r - yd)

       # plant (time continuous part)
       u ~ Hold(ud)
       D(x) ~ -x + u
       y ~ x]
@named sys = ODESystem(eqs)
# compute equation and variables' time domains

ts = TearingState(sys)
ci = ModelingToolkit.ClockInference(ts)
ModelingToolkit.infer_clocks!(ci)
varmap = Dict(ci.ts.fullvars .=> ci.var_domain)
eqmap = ci.eq_domain

d = Clock(t, dt)
# Note that TearingState reorders the equations
@test eqmap[1] == Continuous()
@test eqmap[2] == d
@test eqmap[3] == d
@test eqmap[4] == Continuous()
@test eqmap[5] == Continuous()

@test varmap[yd] == d
@test varmap[ud] == d
@test varmap[r] == d
@test varmap[x] == Continuous()
@test varmap[y] == Continuous()
@test varmap[u] == Continuous()
