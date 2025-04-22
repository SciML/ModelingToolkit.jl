using Test
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@connector function TwoPhaseFluidPort(; name, P = 0.0, m_flow = 0.0, h_outflow = 0.0)
    pars = @parameters begin
        rho
        bulk
        viscosity
    end

    vars = @variables begin
        (h_outflow(t) = h_outflow), [connect = Stream]
        (m_flow(t) = m_flow), [connect = Flow]
        P(t) = P
    end

    ODESystem(Equation[], t, vars, pars; name = name)
end

@connector function TwoPhaseFluid(; name, R, B, V)
    pars = @parameters begin
        rho = R
        bulk = B
        viscosity = V
    end

    vars = @variables begin
        m_flow(t), [connect = Flow]
    end

    # equations ---------------------------
    eqs = Equation[m_flow ~ 0]

    ODESystem(eqs, t, vars, pars; name)
end

function MassFlowSource_h(; name,
        h_in = 420e3,
        m_flow_in = -0.01)
    pars = @parameters begin
        h_in = h_in
        m_flow_in = m_flow_in
    end

    vars = @variables begin
        P(t)
    end

    @named port = TwoPhaseFluidPort()

    subs = [port]

    eqns = Equation[]

    push!(eqns, port.P ~ P)
    push!(eqns, port.m_flow ~ -m_flow_in)
    push!(eqns, port.h_outflow ~ h_in)

    compose(ODESystem(eqns, t, vars, pars; name = name), subs)
end

# Simplified components.
function AdiabaticStraightPipe(; name,
        kwargs...)
    vars = []
    pars = []

    @named port_a = TwoPhaseFluidPort()
    @named port_b = TwoPhaseFluidPort()

    subs = [port_a; port_b]

    eqns = Equation[]

    push!(eqns, connect(port_a, port_b))
    sys = ODESystem(eqns, t, vars, pars; name = name)
    sys = compose(sys, subs)
end

function SmallBoundary_Ph(; name,
        P_in = 1e6,
        h_in = 400e3)
    vars = []

    pars = @parameters begin
        P = P_in
        h = h_in
    end

    @named port1 = TwoPhaseFluidPort()

    subs = [port1]

    eqns = Equation[]

    push!(eqns, port1.P ~ P)
    push!(eqns, port1.h_outflow ~ h)

    compose(ODESystem(eqns, t, vars, pars; name = name), subs)
end

# N1M1 model and test code.
function N1M1(; name,
        P_in = 1e6,
        h_in = 400e3,
        kwargs...)
    @named port_a = TwoPhaseFluidPort()
    @named source = SmallBoundary_Ph(P_in = P_in, h_in = h_in)

    subs = [port_a; source]

    eqns = Equation[]

    push!(eqns, connect(source.port1, port_a))

    sys = ODESystem(eqns, t, [], [], name = name)
    sys = compose(sys, subs)
end

@named fluid = TwoPhaseFluid(; R = 876, B = 1.2e9, V = 0.034)
@named n1m1 = N1M1()
@named pipe = AdiabaticStraightPipe()
@named sink = MassFlowSource_h(m_flow_in = -0.01, h_in = 400e3)

eqns = [connect(n1m1.port_a, pipe.port_a)
        connect(pipe.port_b, sink.port)]

@named sys = ODESystem(eqns, t)

eqns = [domain_connect(fluid, n1m1.port_a)
        connect(n1m1.port_a, pipe.port_a)
        connect(pipe.port_b, sink.port)]

@named n1m1Test = ODESystem(eqns, t, [], []; systems = [fluid, n1m1, pipe, sink])

@test_nowarn structural_simplify(n1m1Test)
@unpack source, port_a = n1m1
ssort(eqs) = sort(eqs, by = string)
@test ssort(equations(expand_connections(n1m1))) == ssort([0 ~ port_a.m_flow
             0 ~ source.port1.m_flow - port_a.m_flow
             source.port1.P ~ port_a.P
             source.port1.P ~ source.P
             source.port1.h_outflow ~ port_a.h_outflow
             source.port1.h_outflow ~ source.h])
@unpack port_a, port_b = pipe
@test ssort(equations(expand_connections(pipe))) ==
      ssort([0 ~ -port_a.m_flow - port_b.m_flow
             0 ~ port_a.m_flow
             0 ~ port_b.m_flow
             port_a.P ~ port_b.P
             port_a.h_outflow ~ instream(port_b.h_outflow)
             port_b.h_outflow ~ instream(port_a.h_outflow)])
@test ssort(equations(expand_connections(sys))) ==
      ssort([0 ~ n1m1.port_a.m_flow + pipe.port_a.m_flow
             0 ~ pipe.port_b.m_flow + sink.port.m_flow
             n1m1.port_a.P ~ pipe.port_a.P
             pipe.port_b.P ~ sink.port.P])
@test ssort(equations(expand_connections(n1m1Test))) ==
      ssort([0 ~ -pipe.port_a.m_flow - pipe.port_b.m_flow
             0 ~ n1m1.source.port1.m_flow - n1m1.port_a.m_flow
             0 ~ n1m1.port_a.m_flow + pipe.port_a.m_flow
             0 ~ pipe.port_b.m_flow + sink.port.m_flow
             fluid.m_flow ~ 0
             n1m1.port_a.P ~ pipe.port_a.P
             n1m1.source.port1.P ~ n1m1.port_a.P
             n1m1.source.port1.P ~ n1m1.source.P
             n1m1.source.port1.h_outflow ~ n1m1.port_a.h_outflow
             n1m1.source.port1.h_outflow ~ n1m1.source.h
             pipe.port_a.P ~ pipe.port_b.P
             pipe.port_a.h_outflow ~ sink.port.h_outflow
             pipe.port_b.P ~ sink.port.P
             pipe.port_b.h_outflow ~ n1m1.port_a.h_outflow
             sink.port.P ~ sink.P
             sink.port.h_outflow ~ sink.h_in
             sink.port.m_flow ~ -sink.m_flow_in])

# N1M2 model and test code.
function N1M2(; name,
        P_in = 1e6,
        h_in = 400e3,
        kwargs...)
    @named port_a = TwoPhaseFluidPort()
    @named port_b = TwoPhaseFluidPort()

    @named source = SmallBoundary_Ph(P_in = P_in, h_in = h_in)

    subs = [port_a; port_b; source]

    eqns = Equation[]

    push!(eqns, connect(source.port1, port_a))
    push!(eqns, connect(source.port1, port_b))

    sys = ODESystem(eqns, t, [], [], name = name)
    sys = compose(sys, subs)
end

@named n1m2 = N1M2()
@named sink1 = MassFlowSource_h(m_flow_in = -0.01, h_in = 400e3)
@named sink2 = MassFlowSource_h(m_flow_in = -0.01, h_in = 400e3)

eqns = [connect(n1m2.port_a, sink1.port)
        connect(n1m2.port_b, sink2.port)]

@named sys = ODESystem(eqns, t)
@named n1m2Test = compose(sys, n1m2, sink1, sink2)
@test_nowarn structural_simplify(n1m2Test)

@named n1m2 = N1M2()
@named pipe1 = AdiabaticStraightPipe()
@named pipe2 = AdiabaticStraightPipe()
@named sink1 = MassFlowSource_h(m_flow_in = -0.01, h_in = 400e3)
@named sink2 = MassFlowSource_h(m_flow_in = -0.01, h_in = 400e3)

eqns = [connect(n1m2.port_a, pipe1.port_a)
        connect(pipe1.port_b, sink1.port)
        connect(n1m2.port_b, pipe2.port_a)
        connect(pipe2.port_b, sink2.port)]

@named sys = ODESystem(eqns, t)
@named n1m2AltTest = compose(sys, n1m2, pipe1, pipe2, sink1, sink2)
@test_nowarn structural_simplify(n1m2AltTest)

# N2M2 model and test code.
function N2M2(; name,
        kwargs...)
    @named port_a = TwoPhaseFluidPort()
    @named port_b = TwoPhaseFluidPort()
    @named pipe = AdiabaticStraightPipe()

    subs = [port_a; port_b; pipe]

    eqns = Equation[]

    push!(eqns, connect(port_a, pipe.port_a))
    push!(eqns, connect(pipe.port_b, port_b))

    sys = ODESystem(eqns, t, [], [], name = name)
    sys = compose(sys, subs)
end

@named n2m2 = N2M2()
@named source = MassFlowSource_h(m_flow_in = -0.01, h_in = 400e3)
@named sink = SmallBoundary_Ph(P_in = 1e6, h_in = 400e3)

eqns = [connect(source.port, n2m2.port_a)
        connect(n2m2.port_b, sink.port1)]

@named sys = ODESystem(eqns, t)
@named n2m2Test = compose(sys, n2m2, source, sink)
@test_nowarn structural_simplify(n2m2Test)

# stream var
@named sp1 = TwoPhaseFluidPort()
@named sp2 = TwoPhaseFluidPort()
@named sys = ODESystem([connect(sp1, sp2)], t)
sys_exp = expand_connections(compose(sys, [sp1, sp2]))
@test ssort(equations(sys_exp)) == ssort([0 ~ -sp1.m_flow - sp2.m_flow
             0 ~ sp1.m_flow
             0 ~ sp2.m_flow
             sp1.P ~ sp2.P
             sp1.h_outflow ~ ModelingToolkit.instream(sp2.h_outflow)
             sp2.h_outflow ~ ModelingToolkit.instream(sp1.h_outflow)])

# array var
@connector function VecPin(; name)
    sts = @variables v(t)[1:2]=[1.0, 0.0] i(t)[1:2]=1.0 [connect = Flow]
    ODESystem(Equation[], t, [sts...;], []; name = name)
end

@named vp1 = VecPin()
@named vp2 = VecPin()
@named vp3 = VecPin()

@named simple = ODESystem([connect(vp1, vp2, vp3)], t)
sys = expand_connections(compose(simple, [vp1, vp2, vp3]))
@test ssort(equations(sys)) == ssort([0 .~ collect(vp1.i)
             0 .~ collect(vp2.i)
             0 .~ collect(vp3.i)
             vp1.v[1] ~ vp2.v[1]
             vp1.v[2] ~ vp2.v[2]
             vp1.v[1] ~ vp3.v[1]
             vp1.v[2] ~ vp3.v[2]
             0 ~ -vp1.i[1] - vp2.i[1] - vp3.i[1]
             0 ~ -vp1.i[2] - vp2.i[2] - vp3.i[2]])

@connector function VectorHeatPort(; name, N = 100, T0 = 0.0, Q0 = 0.0)
    @variables (T(t))[1:N]=T0 (Q(t))[1:N]=Q0 [connect = Flow]
    ODESystem(Equation[], t, [T; Q], []; name = name)
end

@test_nowarn @named a = VectorHeatPort()

# --------------------------------------------------
# Test the new Domain feature

sys_ = expand_connections(n1m1Test)
sys_defs = ModelingToolkit.defaults(sys_)
csys = complete(n1m1Test)
@test Symbol(sys_defs[csys.pipe.port_a.rho]) == Symbol(csys.fluid.rho)
@test Symbol(sys_defs[csys.pipe.port_b.rho]) == Symbol(csys.fluid.rho)

# Testing the domain feature with non-stream system...

@connector function HydraulicPort(; P, name)
    pars = @parameters begin
        p_int = P
        rho
        bulk
        viscosity
    end

    vars = @variables begin
        p(t) = p_int
        dm(t), [connect = Flow]
    end

    # equations ---------------------------
    eqs = Equation[]

    ODESystem(eqs, t, vars, pars; name, defaults = [dm => 0])
end

@connector function Fluid(; name, R, B, V)
    pars = @parameters begin
        rho = R
        bulk = B
        viscosity = V
    end

    vars = @variables begin
        dm(t), [connect = Flow]
    end

    # equations ---------------------------
    eqs = [
        dm ~ 0
    ]

    ODESystem(eqs, t, vars, pars; name)
end

function StepSource(; P, name)
    pars = @parameters begin
        p_int = P
    end

    vars = []

    # nodes -------------------------------
    systems = @named begin
        H = HydraulicPort(; P = p_int)
    end

    # equations ---------------------------
    eqs = [
        H.p ~ p_int * (t > 0.01)
    ]

    ODESystem(eqs, t, vars, pars; name, systems)
end

function StaticVolume(; P, V, name)
    pars = @parameters begin
        p_int = P
        vol = V
    end

    vars = @variables begin
        p(t) = p_int
        vrho(t)
        drho(t) = 0
    end

    # nodes -------------------------------
    systems = @named begin
        H = HydraulicPort(; P = p_int)
    end

    # fluid props ------------------------
    rho_0 = H.rho

    # equations ---------------------------
    eqs = [D(vrho) ~ drho
           vrho ~ rho_0 * (1 + p / H.bulk)
           H.p ~ p
           H.dm ~ drho * V]

    ODESystem(eqs, t, vars, pars; name, systems,
        defaults = [vrho => rho_0 * (1 + p_int / H.bulk)])
end

function PipeBase(; P, R, name)
    pars = @parameters begin
        p_int = P
        resistance = R
    end

    vars = []

    # nodes -------------------------------
    systems = @named begin
        HA = HydraulicPort(; P = p_int)
        HB = HydraulicPort(; P = p_int)
    end

    # equations ---------------------------
    eqs = [HA.p - HB.p ~ HA.dm * resistance / HA.viscosity
           0 ~ HA.dm + HB.dm
           domain_connect(HA, HB)]

    ODESystem(eqs, t, vars, pars; name, systems)
end

function Pipe(; P, R, name)
    pars = @parameters begin
        p_int = P
        resistance = R
    end

    vars = []

    systems = @named begin
        HA = HydraulicPort(; P = p_int)
        HB = HydraulicPort(; P = p_int)
        p12 = PipeBase(; P = p_int, R = resistance)
        v1 = StaticVolume(; P = p_int, V = 0.01)
        v2 = StaticVolume(; P = p_int, V = 0.01)
    end

    eqs = [connect(v1.H, p12.HA, HA)
           connect(v2.H, p12.HB, HB)]

    ODESystem(eqs, t, vars, pars; name, systems)
end

function TwoFluidSystem(; name)
    pars = []
    vars = []

    # nodes -------------------------------
    systems = @named begin
        fluid_a = Fluid(; R = 876, B = 1.2e9, V = 0.034)
        source_a = StepSource(; P = 10e5)
        pipe_a = Pipe(; P = 0, R = 1e6)
        volume_a = StaticVolume(; P = 0, V = 0.1)

        fluid_b = Fluid(; R = 1000, B = 2.5e9, V = 0.00034)
        source_b = StepSource(; P = 10e5)
        pipe_b = Pipe(; P = 0, R = 1e6)
        volume_b = StaticVolume(; P = 0, V = 0.1)
    end

    # equations ---------------------------
    eqs = [connect(fluid_a, source_a.H)
           connect(source_a.H, pipe_a.HA)
           connect(pipe_a.HB, volume_a.H)
           connect(fluid_b, source_b.H)
           connect(source_b.H, pipe_b.HA)
           connect(pipe_b.HB, volume_b.H)]

    ODESystem(eqs, t, vars, pars; name, systems)
end

@named two_fluid_system = TwoFluidSystem()
sys = expand_connections(two_fluid_system)

sys_defs = ModelingToolkit.defaults(sys)
csys = complete(two_fluid_system)

@test Symbol(sys_defs[csys.volume_a.H.rho]) == Symbol(csys.fluid_a.rho)
@test Symbol(sys_defs[csys.volume_b.H.rho]) == Symbol(csys.fluid_b.rho)

@test_nowarn structural_simplify(two_fluid_system)

function OneFluidSystem(; name)
    pars = []
    vars = []

    # nodes -------------------------------
    systems = @named begin
        fluid = Fluid(; R = 876, B = 1.2e9, V = 0.034)

        source_a = StepSource(; P = 10e5)
        pipe_a = Pipe(; P = 0, R = 1e6)
        volume_a = StaticVolume(; P = 0, V = 0.1)

        source_b = StepSource(; P = 20e5)
        pipe_b = Pipe(; P = 0, R = 1e6)
        volume_b = StaticVolume(; P = 0, V = 0.1)
    end

    # equations ---------------------------
    eqs = [connect(fluid, source_a.H, source_b.H)
           connect(source_a.H, pipe_a.HA)
           connect(pipe_a.HB, volume_a.H)
           connect(source_b.H, pipe_b.HA)
           connect(pipe_b.HB, volume_b.H)]

    ODESystem(eqs, t, vars, pars; name, systems)
end

@named one_fluid_system = OneFluidSystem()
sys = expand_connections(one_fluid_system)

sys_defs = ModelingToolkit.defaults(sys)
csys = complete(one_fluid_system)

@test Symbol(sys_defs[csys.volume_a.H.rho]) == Symbol(csys.fluid.rho)
@test Symbol(sys_defs[csys.volume_b.H.rho]) == Symbol(csys.fluid.rho)

@test_nowarn structural_simplify(one_fluid_system)
