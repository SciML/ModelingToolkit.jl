using Test
using ModelingToolkit
@variables t

@connector function TwoPhaseFluidPort(; name, P = 0.0, m_flow = 0.0, h_outflow = 0.0)
    vars = @variables h_outflow(t)=h_outflow [connect = Stream] m_flow(t)=m_flow [
        connect = Flow,
    ] P(t)=P
    ODESystem(Equation[], t, vars, []; name = name)
end

function MassFlowSource_h(; name,
                          h_in = 420e3,
                          m_flow_in = -0.01)
    pars = @parameters begin
        h_in = h_in
        m_flow_in = m_flow_in
    end

    vars = @variables begin P(t) end

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

@named n1m1 = N1M1()
@named pipe = AdiabaticStraightPipe()
@named sink = MassFlowSource_h(m_flow_in = -0.01, h_in = 400e3)

eqns = [connect(n1m1.port_a, pipe.port_a)
        connect(pipe.port_b, sink.port)]

@named sys = ODESystem(eqns, t)
@named n1m1Test = compose(sys, n1m1, pipe, sink)
@test_nowarn structural_simplify(n1m1Test)
@unpack source, port_a = n1m1
@test sort(equations(expand_connections(n1m1)), by = string) == [0 ~ port_a.m_flow
       0 ~ source.port1.m_flow - port_a.m_flow
       source.port1.P ~ port_a.P
       source.port1.P ~ source.P
       source.port1.h_outflow ~ port_a.h_outflow
       source.port1.h_outflow ~ source.h]
@unpack port_a, port_b = pipe
@test sort(equations(expand_connections(pipe)), by = string) ==
      [0 ~ -port_a.m_flow - port_b.m_flow
       0 ~ port_a.m_flow
       0 ~ port_b.m_flow
       port_a.P ~ port_b.P
       port_a.h_outflow ~ instream(port_b.h_outflow)
       port_b.h_outflow ~ instream(port_a.h_outflow)]
@test sort(equations(expand_connections(sys)), by = string) ==
      [0 ~ n1m1.port_a.m_flow + pipe.port_a.m_flow
       0 ~ pipe.port_b.m_flow + sink.port.m_flow
       n1m1.port_a.P ~ pipe.port_a.P
       pipe.port_b.P ~ sink.port.P]
@test sort(equations(expand_connections(n1m1Test)), by = string) ==
      [0 ~ -pipe.port_a.m_flow - pipe.port_b.m_flow
       0 ~ n1m1.port_a.m_flow + pipe.port_a.m_flow
       0 ~ n1m1.source.port1.m_flow - n1m1.port_a.m_flow
       0 ~ pipe.port_b.m_flow + sink.port.m_flow
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
       sink.port.m_flow ~ -sink.m_flow_in]

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
@test sort(equations(sys_exp), by = string) == [0 ~ -sp1.m_flow - sp2.m_flow
       0 ~ sp1.m_flow
       0 ~ sp2.m_flow
       sp1.P ~ sp2.P
       sp1.h_outflow ~ ModelingToolkit.instream(sp2.h_outflow)
       sp2.h_outflow ~ ModelingToolkit.instream(sp1.h_outflow)]

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
@test sort(equations(sys), by = string) == sort([0 .~ collect(vp1.i)
            0 .~ collect(vp2.i)
            0 .~ collect(vp3.i)
            vp1.v[1] ~ vp2.v[1]
            vp1.v[2] ~ vp2.v[2]
            vp1.v[1] ~ vp3.v[1]
            vp1.v[2] ~ vp3.v[2]
            0 ~ -vp1.i[1] - vp2.i[1] - vp3.i[1]
            0 ~ -vp1.i[2] - vp2.i[2] - vp3.i[2]], by = string)

@connector function VectorHeatPort(; name, N = 100, T0 = 0.0, Q0 = 0.0)
    @variables (T(t))[1:N]=T0 (Q(t))[1:N]=Q0 [connect = Flow]
    ODESystem(Equation[], t, [T; Q], []; name = name)
end

@test_nowarn @named a = VectorHeatPort()
