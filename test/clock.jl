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
#TODO: test linearize

#=
 Differential(t)(x(t)) ~ u(t) - x(t)
 0 ~ Sample(Clock(t, 0.1))(y(t)) - yd(t)
 0 ~ kp*(r(t) - yd(t)) - ud(t)
 0 ~ Hold()(ud(t)) - u(t)
 0 ~ x(t) - y(t)

====
By inference:

 Differential(t)(x(t)) ~ u(t) - x(t)
 0 ~ Hold()(ud(t)) - u(t) # Hold()(ud(t)) is constant except in an event
 0 ~ x(t) - y(t)

 0 ~ Sample(Clock(t, 0.1))(y(t)) - yd(t)
 0 ~ kp*(r(t) - yd(t)) - ud(t)

====

 Differential(t)(x(t)) ~ u(t) - x(t)
 0 ~ Hold()(ud(t)) - u(t)
 0 ~ x(t) - y(t)

 yd(t) := Sample(Clock(t, 0.1))(y(t))
 ud(t) := kp*(r(t) - yd(t))
=#

#=
     D(x) ~ Shift(x, 0, dt) + 1 # this should never meet with continous variables
=>   (Shift(x, 0, dt) - Shift(x, -1, dt))/dt ~ Shift(x, 0, dt) + 1
=>   Shift(x, 0, dt) - Shift(x, -1, dt) ~ Shift(x, 0, dt) * dt + dt
=>   Shift(x, 0, dt) - Shift(x, 0, dt) * dt ~ Shift(x, -1, dt) + dt
=>   (1 - dt) * Shift(x, 0, dt) ~ Shift(x, -1, dt) + dt
=>   Shift(x, 0, dt) := (Shift(x, -1, dt) + dt) / (1 - dt) # Discrete system
=#

ci, varmap = infer_clocks(sys)
eqmap = ci.eq_domain
tss, io = ModelingToolkit.split_system(deepcopy(ci))
ts_c = deepcopy(tss[1])
@set! ts_c.structure.solvable_graph = nothing
sss, = ModelingToolkit.structural_simplify!(ts_c, io)
@test equations(sss) == [D(x) ~ u - x]

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

@info "test composed systems"

dt = 0.5
@variables t
d = Clock(t, dt)
k = ShiftIndex(d)
timevec = 0:0.1:4

function plant(; name)
    @variables x(t)=1 u(t)=0 y(t)=0
    D = Differential(t)
    eqs = [D(x) ~ -x + u
           y ~ x]
    ODESystem(eqs, t; name = name)
end

function filt(; name)
    @variables x(t)=0 u(t)=0 y(t)=0
    a = 1 / exp(dt)
    eqs = [x(k + 1) ~ a * x + (1 - a) * u(k)
           y ~ x]
    ODESystem(eqs, t, name = name)
end

function controller(kp; name)
    @variables y(t)=0 r(t)=0 ud(t)=0 yd(t)=0
    @parameters kp = kp
    eqs = [yd ~ Sample(y)
           ud ~ kp * (r - yd)]
    ODESystem(eqs, t; name = name)
end

@named f = filt()
@named c = controller(1)
@named p = plant()

connections = [f.u ~ -1#(t >= 1)  # step input
               f.y ~ c.r # filtered reference to controller reference
               Hold(c.ud) ~ p.u # controller output to plant input
               p.y ~ c.y]

@named cl = ODESystem(connections, t, systems = [f, c, p])

ci, varmap = infer_clocks(cl)

@test varmap[f.x] == Clock(t, 0.5)
@test varmap[p.x] == Continuous()
@test varmap[p.y] == Continuous()
@test varmap[c.ud] == Clock(t, 0.5)
@test varmap[c.yd] == Clock(t, 0.5)
@test varmap[c.y] == Continuous()
@test varmap[f.y] == Clock(t, 0.5)
@test varmap[f.u] == Clock(t, 0.5)
@test varmap[p.u] == Continuous()
@test varmap[c.r] == Clock(t, 0.5)
