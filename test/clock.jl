using ModelingToolkit, Test, Setfield, OrdinaryDiffEq, DiffEqCallbacks
using ModelingToolkit: Continuous
using ModelingToolkit: t_nounits as t, D_nounits as D

function infer_clocks(sys)
    ts = TearingState(sys)
    ci = ModelingToolkit.ClockInference(ts)
    ModelingToolkit.infer_clocks!(ci), Dict(ci.ts.fullvars .=> ci.var_domain)
end

@info "Testing hybrid system"
dt = 0.1
@variables x(t) y(t) u(t) yd(t) ud(t) r(t)
@parameters kp
# u(n + 1) := f(u(n))

eqs = [yd ~ Sample(dt)(y)
       ud ~ kp * (r - yd)
       r ~ 1.0

       # plant (time continuous part)
       u ~ Hold(ud)
       D(x) ~ -x + u
       y ~ x]
@named sys = ODESystem(eqs, t)
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
     D(x) ~ Shift(x, 0, dt) + 1 # this should never meet with continuous variables
=>   (Shift(x, 0, dt) - Shift(x, -1, dt))/dt ~ Shift(x, 0, dt) + 1
=>   Shift(x, 0, dt) - Shift(x, -1, dt) ~ Shift(x, 0, dt) * dt + dt
=>   Shift(x, 0, dt) - Shift(x, 0, dt) * dt ~ Shift(x, -1, dt) + dt
=>   (1 - dt) * Shift(x, 0, dt) ~ Shift(x, -1, dt) + dt
=>   Shift(x, 0, dt) := (Shift(x, -1, dt) + dt) / (1 - dt) # Discrete system
=#

ci, varmap = infer_clocks(sys)
eqmap = ci.eq_domain
tss, inputs, continuous_id = ModelingToolkit.split_system(deepcopy(ci))
sss, = ModelingToolkit._structural_simplify!(
    deepcopy(tss[continuous_id]), (inputs[continuous_id], ()))
@test equations(sss) == [D(x) ~ u - x]
sss, = ModelingToolkit._structural_simplify!(deepcopy(tss[1]), (inputs[1], ()))
@test isempty(equations(sss))
d = Clock(dt)
k = ShiftIndex(d)
@test observed(sss) == [yd(k + 1) ~ Sample(dt)(y); r(k + 1) ~ 1.0;
       ud(k + 1) ~ kp * (r(k + 1) - yd(k + 1))]

d = Clock(dt)
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
@variables x(t) y(t) u(t) yd(t) ud(t)
@parameters kp
d = Clock(dt)
k = ShiftIndex(d)

eqs = [yd ~ Sample(dt)(y)
       ud ~ kp * yd + ud(k - 2)

       # plant (time continuous part)
       u ~ Hold(ud)
       D(x) ~ -x + u
       y ~ x]
@named sys = ODESystem(eqs, t)
@test_throws ModelingToolkit.HybridSystemNotSupportedException ss=structural_simplify(sys);

@test_skip begin
    Tf = 1.0
    prob = ODEProblem(ss, [x => 0.1], (0.0, Tf),
        [kp => 1.0; ud(k - 1) => 2.1; ud(k - 2) => 2.0])
    # create integrator so callback is evaluated at t=0 and we can test correct param values
    int = init(prob, Tsit5(); kwargshandle = KeywordArgSilent)
    @test sort(vcat(int.p...)) == [0.1, 1.0, 2.1, 2.1, 2.1] # yd, kp, ud(k-1), ud, Hold(ud)
    prob = ODEProblem(ss, [x => 0.1], (0.0, Tf),
        [kp => 1.0; ud(k - 1) => 2.1; ud(k - 2) => 2.0]) # recreate problem to empty saved values
    sol = solve(prob, Tsit5(), kwargshandle = KeywordArgSilent)

    ss_nosplit = structural_simplify(sys; split = false)
    prob_nosplit = ODEProblem(ss_nosplit, [x => 0.1], (0.0, Tf),
        [kp => 1.0; ud(k - 1) => 2.1; ud(k - 2) => 2.0])
    int = init(prob_nosplit, Tsit5(); kwargshandle = KeywordArgSilent)
    @test sort(int.p) == [0.1, 1.0, 2.1, 2.1, 2.1] # yd, kp, ud(k-1), ud, Hold(ud)
    prob_nosplit = ODEProblem(ss_nosplit, [x => 0.1], (0.0, Tf),
        [kp => 1.0; ud(k - 1) => 2.1; ud(k - 2) => 2.0]) # recreate problem to empty saved values
    sol_nosplit = solve(prob_nosplit, Tsit5(), kwargshandle = KeywordArgSilent)
    # For all inputs in parameters, just initialize them to 0.0, and then set them
    # in the callback.

    # kp is the only real parameter
    function foo!(du, u, p, t)
        x = u[1]
        ud = p[2]
        du[1] = -x + ud
    end
    function affect!(integrator, saved_values)
        yd = integrator.u[1]
        kp = integrator.p[1]
        ud = integrator.p[2]
        udd = integrator.p[3]

        integrator.p[2] = kp * yd + udd
        integrator.p[3] = ud

        push!(saved_values.t, integrator.t)
        push!(saved_values.saveval, [integrator.p[2], integrator.p[3]])

        nothing
    end
    saved_values = SavedValues(Float64, Vector{Float64})
    cb = PeriodicCallback(
        Base.Fix2(affect!, saved_values), 0.1; final_affect = true, initial_affect = true)
    #                                           kp   ud
    prob = ODEProblem(foo!, [0.1], (0.0, Tf), [1.0, 2.1, 2.0], callback = cb)
    sol2 = solve(prob, Tsit5())
    @test sol.u == sol2.u
    @test sol_nosplit.u == sol2.u
    @test saved_values.t == sol.prob.kwargs[:disc_saved_values][1].t
    @test saved_values.t == sol_nosplit.prob.kwargs[:disc_saved_values][1].t
    @test saved_values.saveval == sol.prob.kwargs[:disc_saved_values][1].saveval
    @test saved_values.saveval == sol_nosplit.prob.kwargs[:disc_saved_values][1].saveval

    @info "Testing multi-rate hybrid system"
    dt = 0.1
    dt2 = 0.2
    @variables x(t) y(t) u(t) r(t) yd1(t) ud1(t) yd2(t) ud2(t)
    @parameters kp

    eqs = [
           # controller (time discrete part `dt=0.1`)
           yd1 ~ Sample(dt)(y)
           ud1 ~ kp * (Sample(dt)(r) - yd1)
           yd2 ~ Sample(dt2)(y)
           ud2 ~ kp * (Sample(dt2)(r) - yd2)

           # plant (time continuous part)
           u ~ Hold(ud1) + Hold(ud2)
           D(x) ~ -x + u
           y ~ x]
    @named sys = ODESystem(eqs, t)
    ci, varmap = infer_clocks(sys)

    d = Clock(dt)
    d2 = Clock(dt2)
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
    d = Clock(dt)
    k = ShiftIndex(d)
    timevec = 0:0.1:4

    function plant(; name)
        @variables x(t)=1 u(t)=0 y(t)=0
        eqs = [D(x) ~ -x + u
               y ~ x]
        ODESystem(eqs, t; name = name)
    end

    function filt(; name)
        @variables x(t)=0 u(t)=0 y(t)=0
        a = 1 / exp(dt)
        eqs = [x ~ a * x(k - 1) + (1 - a) * u(k - 1)
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

    @test varmap[f.x] == Clock(0.5)
    @test varmap[p.x] == Continuous()
    @test varmap[p.y] == Continuous()
    @test varmap[c.ud] == Clock(0.5)
    @test varmap[c.yd] == Clock(0.5)
    @test varmap[c.y] == Continuous()
    @test varmap[f.y] == Clock(0.5)
    @test varmap[f.u] == Clock(0.5)
    @test varmap[p.u] == Continuous()
    @test varmap[c.r] == Clock(0.5)

    ## Multiple clock rates
    @info "Testing multi-rate hybrid system"
    dt = 0.1
    dt2 = 0.2
    @variables x(t)=0 y(t)=0 u(t)=0 yd1(t)=0 ud1(t)=0 yd2(t)=0 ud2(t)=0
    @parameters kp=1 r=1

    eqs = [
           # controller (time discrete part `dt=0.1`)
           yd1 ~ Sample(dt)(y)
           ud1 ~ kp * (r - yd1)
           # controller (time discrete part `dt=0.2`)
           yd2 ~ Sample(dt2)(y)
           ud2 ~ kp * (r - yd2)

           # plant (time continuous part)
           u ~ Hold(ud1) + Hold(ud2)
           D(x) ~ -x + u
           y ~ x]

    @named cl = ODESystem(eqs, t)

    d = Clock(dt)
    d2 = Clock(dt2)

    ci, varmap = infer_clocks(cl)
    @test varmap[yd1] == d
    @test varmap[ud1] == d
    @test varmap[yd2] == d2
    @test varmap[ud2] == d2
    @test varmap[x] == Continuous()
    @test varmap[y] == Continuous()
    @test varmap[u] == Continuous()

    ss = structural_simplify(cl)
    ss_nosplit = structural_simplify(cl; split = false)

    if VERSION >= v"1.7"
        prob = ODEProblem(ss, [x => 0.0], (0.0, 1.0), [kp => 1.0])
        prob_nosplit = ODEProblem(ss_nosplit, [x => 0.0], (0.0, 1.0), [kp => 1.0])
        sol = solve(prob, Tsit5(), kwargshandle = KeywordArgSilent)
        sol_nosplit = solve(prob_nosplit, Tsit5(), kwargshandle = KeywordArgSilent)

        function foo!(dx, x, p, t)
            kp, ud1, ud2 = p
            dx[1] = -x[1] + ud1 + ud2
        end

        function affect1!(integrator)
            kp = integrator.p[1]
            y = integrator.u[1]
            r = 1.0
            ud1 = kp * (r - y)
            integrator.p[2] = ud1
            nothing
        end
        function affect2!(integrator)
            kp = integrator.p[1]
            y = integrator.u[1]
            r = 1.0
            ud2 = kp * (r - y)
            integrator.p[3] = ud2
            nothing
        end
        cb1 = PeriodicCallback(affect1!, dt; final_affect = true, initial_affect = true)
        cb2 = PeriodicCallback(affect2!, dt2; final_affect = true, initial_affect = true)
        cb = CallbackSet(cb1, cb2)
        #                                           kp   ud1  ud2
        prob = ODEProblem(foo!, [0.0], (0.0, 1.0), [1.0, 1.0, 1.0], callback = cb)
        sol2 = solve(prob, Tsit5())

        @test sol.u≈sol2.u atol=1e-6
        @test sol_nosplit.u≈sol2.u atol=1e-6
    end

    ##
    @info "Testing hybrid system with components"
    using ModelingToolkitStandardLibrary.Blocks

    dt = 0.05
    d = Clock(dt)
    k = ShiftIndex(d)

    @mtkmodel DiscretePI begin
        @components begin
            input = RealInput()
            output = RealOutput()
        end
        @parameters begin
            kp = 1, [description = "Proportional gain"]
            ki = 1, [description = "Integral gain"]
        end
        @variables begin
            x(t) = 0, [description = "Integral state"]
            u(t)
            y(t)
        end
        @equations begin
            x(k) ~ x(k - 1) + ki * u(k) * SampleTime() / dt
            output.u(k) ~ y(k)
            input.u(k) ~ u(k)
            y(k) ~ x(k - 1) + kp * u(k)
        end
    end

    @mtkmodel Sampler begin
        @components begin
            input = RealInput()
            output = RealOutput()
        end
        @equations begin
            output.u ~ Sample(dt)(input.u)
        end
    end

    @mtkmodel ZeroOrderHold begin
        @extend u, y = siso = Blocks.SISO()
        @equations begin
            y ~ Hold(u)
        end
    end

    @mtkmodel ClosedLoop begin
        @components begin
            plant = FirstOrder(k = 0.3, T = 1)
            sampler = Sampler()
            holder = ZeroOrderHold()
            controller = DiscretePI(kp = 2, ki = 2)
            feedback = Feedback()
            ref = Constant(k = 0.5)
        end
        @equations begin
            connect(ref.output, feedback.input1)
            connect(feedback.output, controller.input)
            connect(controller.output, holder.input)
            connect(holder.output, plant.input)
            connect(plant.output, sampler.input)
            connect(sampler.output, feedback.input2)
        end
    end

    ##
    @named model = ClosedLoop()
    _model = complete(model)

    ci, varmap = infer_clocks(expand_connections(_model))

    @test varmap[_model.plant.input.u] == Continuous()
    @test varmap[_model.plant.u] == Continuous()
    @test varmap[_model.plant.x] == Continuous()
    @test varmap[_model.plant.y] == Continuous()
    @test varmap[_model.plant.output.u] == Continuous()
    @test varmap[_model.holder.output.u] == Continuous()
    @test varmap[_model.sampler.input.u] == Continuous()
    @test varmap[_model.controller.u] == d
    @test varmap[_model.holder.input.u] == d
    @test varmap[_model.controller.output.u] == d
    @test varmap[_model.controller.y] == d
    @test varmap[_model.feedback.input1.u] == d
    @test varmap[_model.ref.output.u] == d
    @test varmap[_model.controller.input.u] == d
    @test varmap[_model.controller.x] == d
    @test varmap[_model.sampler.output.u] == d
    @test varmap[_model.feedback.output.u] == d
    @test varmap[_model.feedback.input2.u] == d

    ssys = structural_simplify(model)

    Tf = 0.2
    timevec = 0:(d.dt):Tf

    import ControlSystemsBase as CS
    import ControlSystemsBase: c2d, tf, feedback, lsim
    # z = tf('z', d.dt)
    # P = c2d(tf(0.3, [1, 1]), d.dt)
    P = c2d(CS.ss([-1], [0.3], [1], 0), d.dt)
    C = CS.ss([1], [2], [1], [2], d.dt)

    # Test the output of the continuous partition
    G = feedback(P * C)
    res = lsim(G, (x, t) -> [0.5], timevec)
    y = res.y[:]

    # plant = FirstOrder(k = 0.3, T = 1)
    # controller = DiscretePI(kp = 2, ki = 2)
    # ref = Constant(k = 0.5)

    # ; model.controller.x(k-1) => 0.0
    prob = ODEProblem(ssys,
        [model.plant.x => 0.0; model.controller.kp => 2.0; model.controller.ki => 2.0],
        (0.0, Tf))
    int = init(prob, Tsit5(); kwargshandle = KeywordArgSilent)
    @test_broken int.ps[Hold(ssys.holder.input.u)] == 2 # constant output * kp issue https://github.com/SciML/ModelingToolkit.jl/issues/2356
    @test int.ps[ssys.controller.x] == 1 # c2d
    @test int.ps[Sample(d)(ssys.sampler.input.u)] == 0 # disc state
    sol = solve(prob,
        Tsit5(),
        kwargshandle = KeywordArgSilent,
        abstol = 1e-8,
        reltol = 1e-8)
    @test_skip begin
        # plot([y sol(timevec, idxs = model.plant.output.u)], m = :o, lab = ["CS" "MTK"])

        ##

        @test sol(timevec, idxs = model.plant.output.u)≈y rtol=1e-8 # The output of the continuous partition is delayed exactly one sample

        # Test the output of the discrete partition
        G = feedback(C, P)
        res = lsim(G, (x, t) -> [0.5], timevec)
        y = res.y[:]

        @test_broken sol(timevec .+ 1e-10, idxs = model.controller.output.u)≈y rtol=1e-8 # Broken due to discrete observed
        # plot([y sol(timevec .+ 1e-12, idxs=model.controller.output.u)], lab=["CS" "MTK"])

        # TODO: test the same system, but with the PI controller implemented as
        # x(k) ~ x(k-1) + ki * u
        # y ~ x(k-1) + kp * u
        # Instead. This should be equivalent to the above, but gve me an error when I tried
    end

    ## Test continuous clock

    c = ModelingToolkit.SolverStepClock()
    k = ShiftIndex(c)

    @mtkmodel CounterSys begin
        @variables begin
            count(t) = 0
            u(t) = 0
            ud(t) = 0
        end
        @equations begin
            ud ~ Sample(c)(u)
            count ~ ud(k - 1)
        end
    end

    @mtkmodel FirstOrderSys begin
        @variables begin
            x(t) = 0
        end
        @equations begin
            D(x) ~ -x + sin(t)
        end
    end

    @mtkmodel FirstOrderWithStepCounter begin
        @components begin
            counter = CounterSys()
            firstorder = FirstOrderSys()
        end
        @equations begin
            counter.u ~ firstorder.x
        end
    end

    @mtkbuild model = FirstOrderWithStepCounter()
    prob = ODEProblem(model, [], (0.0, 10.0))
    sol = solve(prob, Tsit5(), kwargshandle = KeywordArgSilent)

    @test sol.prob.kwargs[:disc_saved_values][1].t == sol.t[1:2:end] # Test that the discrete-time system executed at every step of the continuous solver. The solver saves each time step twice, one state value before discrete affect and one after.
    @test_nowarn ModelingToolkit.build_explicit_observed_function(
        model, model.counter.ud)(sol.u[1], prob.p, sol.t[1])

    @variables x(t)=1.0 y(t)=1.0
    eqs = [D(y) ~ Hold(x)
           x ~ x(k - 1) + x(k - 2)]
    @mtkbuild sys = ODESystem(eqs, t)
    prob = ODEProblem(sys, [], (0.0, 10.0))
    int = init(prob, Tsit5(); kwargshandle = KeywordArgSilent)
    @test int.ps[x] == 2.0
    @test int.ps[x(k - 1)] == 1.0

    @test_throws ErrorException ODEProblem(sys, [], (0.0, 10.0), [x => 2.0])
    prob = ODEProblem(sys, [], (0.0, 10.0), [x(k - 1) => 2.0])
    int = init(prob, Tsit5(); kwargshandle = KeywordArgSilent)
    @test int.ps[x] == 3.0
    @test int.ps[x(k - 1)] == 2.0
end
