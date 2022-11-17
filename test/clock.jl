using ModelingToolkit, Test, Setfield, OrdinaryDiffEq, DiffEqCallbacks

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
# u(n + 1) := f(u(n))

eqs = [yd ~ Sample(t, dt)(y)
       ud ~ kp * (r - yd)
       r ~ 1.0

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

using ModelingToolkit.SystemStructures
ci, varmap = infer_clocks(sys)
eqmap = ci.eq_domain
tss, inputs = ModelingToolkit.split_system(deepcopy(ci))
sss, = SystemStructures._structural_simplify!(deepcopy(tss[1]), (inputs[1], ()))
@test equations(sss) == [D(x) ~ u - x]
sss, = SystemStructures._structural_simplify!(deepcopy(tss[2]), (inputs[2], ()))
@test isempty(equations(sss))
@test observed(sss) == [r ~ 1.0; yd ~ Sample(t, dt)(y); ud ~ kp * (r - yd)]

d = Clock(t, dt)
# Note that TearingState reorders the equations
@test eqmap[1] == Continuous()
@test eqmap[2] == d
@test eqmap[3] == d
@test eqmap[4] == d
@test eqmap[5] == Continuous()
@test eqmap[6] == Continuous()

@test varmap[yd] == d
@test varmap[ud] == d
@test varmap[r] == d
@test varmap[x] == Continuous()
@test varmap[y] == Continuous()
@test varmap[u] == Continuous()

@info "Testing shift normalization"
dt = 0.1
@variables t x(t) y(t) u(t) yd(t) ud(t) r(t) z(t)
@parameters kp
D = Differential(t)
d = Clock(t, dt)
k = ShiftIndex(d)

eqs = [yd ~ Sample(t, dt)(y)
       ud ~ kp * (r - yd) + z(k)
       r ~ 1.0

       # plant (time continuous part)
       u ~ Hold(ud)
       D(x) ~ -x + u
       y ~ x
       z(k + 2) ~ z(k) + yd
       #=
       z(k + 2) ~ z(k) + yd
       =>
       z′(k + 1) ~ z(k) + yd
       z(k + 1)  ~ z′(k)
       =#
       ]
@named sys = ODESystem(eqs)
ss = structural_simplify(sys)
prob = ODEProblem(ss, [x => 0.0, y => 0.0], (0.0, 1.0),
                  [kp => 1.0; z => 0.0; z(k + 1) => 0.0])
sol = solve(prob, Tsit5(), kwargshandle = KeywordArgSilent)
# For all inputs in parameters, just initialize them to 0.0, and then set them
# in the callback.

# kp is the only real parameter
function foo!(du, u, p, t)
    x = u[1]
    ud = p[2]
    du[1] = -x + ud
end
function affect!(integrator, saved_values)
    kp = integrator.p[1]
    yd = integrator.u[1]
    z_t = integrator.p[3]
    z = integrator.p[4]
    r = 1.0
    ud = kp * (r - yd) + z
    push!(saved_values.t, integrator.t)
    push!(saved_values.saveval, [integrator.p[4], integrator.p[3]])
    integrator.p[2] = ud
    integrator.p[3] = z + yd
    integrator.p[4] = z_t
    nothing
end
saved_values = SavedValues(Float64, Vector{Float64});
cb = PeriodicCallback(Base.Fix2(affect!, saved_values), 0.1)
prob = ODEProblem(foo!, [0.0], (0.0, 1.0), [1.0, 0.0, 0.0, 0.0], callback = cb)
sol2 = solve(prob, Tsit5())
@test sol.u == sol2.u
@test saved_values.t == sol.prob.kwargs[:disc_saved_values][1].t
@test saved_values.saveval == sol.prob.kwargs[:disc_saved_values][1].saveval

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
