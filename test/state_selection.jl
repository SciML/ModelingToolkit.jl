using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

sts = @variables x1(t) x2(t) x3(t) x4(t)
params = @parameters u1(t) u2(t) u3(t) u4(t)
eqs = [x1 + x2 + u1 ~ 0
       x1 + x2 + x3 + u2 ~ 0
       x1 + D(x3) + x4 + u3 ~ 0
       2 * D(D(x1)) + D(D(x2)) + D(D(x3)) + D(x4) + u4 ~ 0]
@named sys = ODESystem(eqs, t)

let dd = dummy_derivative(sys)
    has_dx1 = has_dx2 = false
    for eq in equations(dd)
        vars = ModelingToolkit.vars(eq)
        has_dx1 |= D(x1) in vars || D(D(x1)) in vars
        has_dx2 |= D(x2) in vars || D(D(x2)) in vars
    end
    @test has_dx1 ⊻ has_dx2 # only one of x1 and x2 can be a dummy derivative
    @test length(unknowns(dd)) == length(equations(dd)) < 9
end

@test_skip let pss = partial_state_selection(sys)
    @test length(equations(pss)) == 1
    @test length(unknowns(pss)) == 2
    @test length(equations(ode_order_lowering(pss))) == 2
end

@parameters σ ρ β
@variables x(t) y(t) z(t) a(t) u(t) F(t)

eqs = [D(x) ~ σ * (y - x)
       D(y) ~ x * (ρ - z) - y + β
       0 ~ z - x + y
       0 ~ a + z
       u ~ z + a]

lorenz1 = ODESystem(eqs, t, name = :lorenz1)
let al1 = alias_elimination(lorenz1)
    let lss = partial_state_selection(al1)
        @test length(equations(lss)) == 2
    end
end

# 1516
let
    @connector function Fluid_port(; name, p = 101325.0, m = 0.0, T = 293.15)
        sts = @variables p(t) [guess = p] m(t) [guess = m, connect = Flow] T(t) [
            guess = T, connect = Stream]
        ODESystem(Equation[], t, sts, []; name = name)
    end

    #this one is for latter
    @connector function Heat_port(; name, Q = 0.0, T = 293.15)
        sts = @variables T(t) [guess = T] Q(t) [guess = Q, connect = Flow]
        ODESystem(Equation[], t, sts, []; name = name)
    end

    # like ground but for fluid systems (fluid_port.m is expected to be zero in closed loop)
    function Compensator(; name, p = 101325.0, T_back = 273.15)
        @named fluid_port = Fluid_port()
        ps = @parameters p=p T_back=T_back
        eqs = [fluid_port.p ~ p
               fluid_port.T ~ T_back]
        compose(ODESystem(eqs, t, [], ps; name = name), fluid_port)
    end

    function Source(; name, delta_p = 100, T_feed = 293.15)
        @named supply_port = Fluid_port() # expected to feed connected pipe -> m<0
        @named return_port = Fluid_port() # expected to receive from connected pipe -> m>0
        ps = @parameters delta_p=delta_p T_feed=T_feed
        eqs = [supply_port.m ~ -return_port.m
               supply_port.p ~ return_port.p + delta_p
               supply_port.T ~ instream(supply_port.T)
               return_port.T ~ T_feed]
        compose(ODESystem(eqs, t, [], ps; name = name), [supply_port, return_port])
    end

    function Substation(; name, T_return = 343.15)
        @named supply_port = Fluid_port() # expected to receive from connected pipe -> m>0
        @named return_port = Fluid_port() # expected to feed connected pipe -> m<0
        ps = @parameters T_return = T_return
        eqs = [supply_port.m ~ -return_port.m
               supply_port.p ~ return_port.p # zero pressure loss for now
               supply_port.T ~ instream(supply_port.T)
               return_port.T ~ T_return]
        compose(ODESystem(eqs, t, [], ps; name = name), [supply_port, return_port])
    end

    function Pipe(; name, L = 1000, d = 0.1, N = 100, rho = 1000, f = 1)
        @named fluid_port_a = Fluid_port()
        @named fluid_port_b = Fluid_port()
        ps = @parameters L=L d=d rho=rho f=f N=N
        sts = @variables v(t) [guess = 0.0] dp_z(t) [guess = 0.0]
        eqs = [fluid_port_a.m ~ -fluid_port_b.m
               fluid_port_a.T ~ instream(fluid_port_a.T)
               fluid_port_b.T ~ fluid_port_a.T
               v * pi * d^2 / 4 * rho ~ fluid_port_a.m
               dp_z ~ abs(v) * v * 0.5 * rho * L / d * f  # pressure loss
               D(v) * rho * L ~ (fluid_port_a.p - fluid_port_b.p - dp_z)]
        compose(ODESystem(eqs, t, sts, ps; name = name), [fluid_port_a, fluid_port_b])
    end
    function System(; name, L = 10.0)
        @named compensator = Compensator()
        @named source = Source()
        @named substation = Substation()
        @named supply_pipe = Pipe(L = L)
        @named return_pipe = Pipe(L = L)
        subs = [compensator, source, substation, supply_pipe, return_pipe]
        ps = @parameters L = L
        eqs = [connect(compensator.fluid_port, source.supply_port)
               connect(source.supply_port, supply_pipe.fluid_port_a)
               connect(supply_pipe.fluid_port_b, substation.supply_port)
               connect(substation.return_port, return_pipe.fluid_port_b)
               connect(return_pipe.fluid_port_a, source.return_port)]
        compose(ODESystem(eqs, t, [], ps; name = name), subs)
    end

    @named system = System(L = 10)
    @unpack supply_pipe, return_pipe = system
    sys = structural_simplify(system)
    u0 = [
        sys.supply_pipe.v => 0.1, sys.return_pipe.v => 0.1, D(supply_pipe.v) => 0.0,
        D(return_pipe.fluid_port_a.m) => 0.0,
        D(supply_pipe.fluid_port_a.m) => 0.0]
    prob1 = ODEProblem(sys, [], (0.0, 10.0), [], guesses = u0)
    prob2 = ODEProblem(sys, [], (0.0, 10.0), [], guesses = u0)
    prob3 = DAEProblem(sys, D.(unknowns(sys)) .=> 0.0, [], (0.0, 10.0), guesses = u0)
    @test solve(prob1, FBDF()).retcode == ReturnCode.Success
    #@test solve(prob2, FBDF()).retcode == ReturnCode.Success
    @test solve(prob3, DFBDF()).retcode == ReturnCode.Success
end

# 1537
let
    @variables begin
        p_1(t)
        p_2(t)
        rho_1(t)
        rho_2(t)
        rho_3(t)
        u_1(t)
        u_2(t)
        u_3(t)
        mo_1(t)
        mo_2(t)
        mo_3(t)
        Ek_1(t)
        Ek_2(t)
        Ek_3(t)
    end

    @parameters dx=100 f=0.3 pipe_D=0.4

    eqs = [p_1 ~ 1.2e5
           p_2 ~ 1e5
           u_1 ~ 10
           mo_1 ~ u_1 * rho_1
           mo_2 ~ u_2 * rho_2
           mo_3 ~ u_3 * rho_3
           Ek_1 ~ rho_1 * u_1 * u_1
           Ek_2 ~ rho_2 * u_2 * u_2
           Ek_3 ~ rho_3 * u_3 * u_3
           rho_1 ~ p_1 / 273.11 / 300
           rho_2 ~ (p_1 + p_2) * 0.5 / 273.11 / 300
           rho_3 ~ p_2 / 273.11 / 300
           D(rho_2) ~ (mo_1 - mo_3) / dx
           D(mo_2) ~ (Ek_1 - Ek_3 + p_1 - p_2) / dx - f / 2 / pipe_D * u_2 * u_2]

    @named trans = ODESystem(eqs, t)

    sys = structural_simplify(trans)

    n = 3
    u = 0 * ones(n)
    rho = 1.2 * ones(n)

    u0 = [p_1 => 1.2e5
          p_2 => 1e5
          u_1 => 0
          u_2 => 0.1
          u_3 => 0.2
          rho_1 => 1.1
          rho_2 => 1.2
          rho_3 => 1.3
          mo_1 => 0
          mo_2 => 1
          mo_3 => 2
          Ek_2 => 2
          Ek_3 => 3]
    prob1 = ODEProblem(sys, [], (0.0, 0.1), guesses = u0)
    prob2 = ODEProblem(sys, [], (0.0, 0.1), guesses = u0)
    @test solve(prob1, FBDF()).retcode == ReturnCode.Success
    @test solve(prob2, FBDF()).retcode == ReturnCode.Success
end

let
    # constant parameters ----------------------------------------------------
    A_1f = 0.0908
    A_2f = 0.036
    p_1f_0 = 1.8e6
    p_2f_0 = p_1f_0 * A_1f / A_2f
    m_total = 3245
    K1 = 4.60425e-5
    K2 = 0.346725
    K3 = 0
    density = 876
    bulk = 1.2e9
    l_1f = 0.7
    x_f_fullscale = 0.025
    p_s = 200e5
    # --------------------------------------------------------------------------

    # modelingtoolkit setup ----------------------------------------------------
    params = @parameters l_2f=0.7 damp=1e3
    vars = @variables begin
        p1(t)
        p2(t)
        dp1(t) = 0
        dp2(t) = 0
        xf(t) = 0
        rho1(t)
        rho2(t)
        drho1(t) = 0
        drho2(t) = 0
        V1(t)
        V2(t)
        dV1(t) = 0
        dV2(t) = 0
        w(t) = 0
        dw(t) = 0
        ddw(t) = 0
    end

    defs = [p1 => p_1f_0
            p2 => p_2f_0
            rho1 => density * (1 + p_1f_0 / bulk)
            rho2 => density * (1 + p_2f_0 / bulk)
            V1 => l_1f * A_1f
            V2 => l_2f * A_2f
            D(p1) => dp1
            D(p2) => dp2
            D(w) => dw
            D(dw) => ddw]

    # equations ------------------------------------------------------------------
    # sqrt -> log as a hack
    flow(x, dp) = K1 * abs(dp) * abs(x) + K2 * log(abs(dp)) * abs(x) + K3 * abs(dp) * x^2
    xm = xf / x_f_fullscale
    Δp1 = p_s - p1
    Δp2 = p2

    eqs = [+flow(xm, Δp1) ~ rho1 * dV1 + drho1 * V1
           0 ~ ifelse(w > 0.5,
               (0) - (rho2 * dV2 + drho2 * V2),
               (-flow(xm, Δp2)) - (rho2 * dV2 + drho2 * V2))
           V1 ~ (l_1f + w) * A_1f
           V2 ~ (l_2f - w) * A_2f
           dV1 ~ +dw * A_1f
           dV2 ~ -dw * A_2f
           rho1 ~ density * (1.0 + p1 / bulk)
           rho2 ~ density * (1.0 + p2 / bulk)
           drho1 ~ density * (dp1 / bulk)
           drho2 ~ density * (dp2 / bulk)
           D(p1) ~ dp1
           D(p2) ~ dp2
           D(w) ~ dw
           D(dw) ~ ddw
           xf ~ 20e-3 * (1 - cos(2 * π * 5 * t))
           0 ~ ifelse(w > 0.5,
               (m_total * ddw) - (p1 * A_1f - p2 * A_2f - damp * dw),
               (m_total * ddw) - (p1 * A_1f - p2 * A_2f))]
    # ----------------------------------------------------------------------------

    # solution -------------------------------------------------------------------
    @named catapult = ODESystem(eqs, t, vars, params, defaults = defs)
    sys = structural_simplify(catapult)
    prob = ODEProblem(sys, [], (0.0, 0.1), [l_2f => 0.55, damp => 1e7]; jac = true)
    @test solve(prob, Rodas4()).retcode == ReturnCode.Success
end
