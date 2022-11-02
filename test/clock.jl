using ModelingToolkit, Test

function infer_clocks(sys)
    ts = TearingState(sys)
    ci = ModelingToolkit.ClockInference(ts)
    ModelingToolkit.infer_clocks!(ci), Dict(ci.ts.fullvars .=> ci.var_domain)
end

@info "Testing hybrid system"
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

ci, varmap = infer_clocks(sys)
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

@info "Testing multi-rate hybrid system"
dt = 0.1
dt2 = 0.2
@variables t x(t) y(t) u(t) r(t) yd1(t) ud1(t) yd2(t) ud2(t)
@parameters kp
D = Differential(t)

eqs = [
       # controller (time discrete part `dt=0.1`)
       yd1 ~ Sample(t, dt)(y)
       ud1 ~ kp * (Sample(t, dt)(r) - yd1)
       yd2 ~ Sample(t, dt2)(y)
       ud2 ~ kp * (Sample(t, dt2)(r) - yd2)

       # plant (time continuous part)
       u ~ Hold(ud1) + Hold(ud2)
       D(x) ~ -x + u
       y ~ x]
@named sys = ODESystem(eqs)
ci, varmap = infer_clocks(sys)

d = Clock(t, dt)
d2 = Clock(t, dt2)
#@test get_eq_domain(eqs[1]) == d
#@test get_eq_domain(eqs[3]) == d2

@test varmap[yd1] == d
@test varmap[ud1] == d
@test varmap[yd2] == d2
@test varmap[ud2] == d2
@test varmap[r] == Continuous()
@test varmap[x] == Continuous()
@test varmap[y] == Continuous()
@test varmap[u] == Continuous()
