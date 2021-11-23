using Test
using ModelingToolkit
@variables t

@connector function TwoPhaseFluidPort(;name, P=0.0, m_flow=0.0, h_outflow=0.0)
    vars = @variables h_outflow(t)=h_outflow [connect=Stream] m_flow(t)=m_flow [connect=Flow] P(t)=P
    ODESystem(Equation[], t, vars, []; name=name)
end

function MassFlowSource_h(;name,
        h_in=420e3,
        m_flow_in=-0.01,
    )

    pars = @parameters begin
        h_in=h_in
        m_flow_in=m_flow_in
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

    compose(ODESystem(eqns, t, vars, pars; name=name), subs)
end

# Simplified components.
function AdiabaticStraightPipe(;name,
        kwargs...,
    )

    vars = []
    pars = []

    @named port_a = TwoPhaseFluidPort()
    @named port_b = TwoPhaseFluidPort()

    subs = [port_a; port_b]

    eqns = Equation[]

    #=
    push!(eqns, port_a.P ~ port_b.P)
    push!(eqns, 0 ~ port_a.m_flow + port_b.m_flow)
    push!(eqns, port_b.h_outflow ~ instream(port_a.h_outflow))
    push!(eqns, port_a.h_outflow ~ instream(port_b.h_outflow))
    =#

    push!(eqns, connect(port_a, port_b))
    sys = ODESystem(eqns, t, vars, pars; name=name)
    sys = compose(sys, subs)
end


function SmallBoundary_Ph(;name,
        P_in=1e6,
        h_in=400e3,
    )

    vars = []

    pars = @parameters begin
        P=P_in
        h=h_in
    end

    @named port1 = TwoPhaseFluidPort()

    subs = [port1]

    eqns = Equation[]

    push!(eqns, port1.P ~ P)
    push!(eqns, port1.h_outflow ~ h)

    compose(ODESystem(eqns, t, vars, pars; name=name), subs)
end


# N1M1 model and test code.
function N1M1(;name,
        P_in=1e6,
        h_in=400e3,
        kwargs...,
    )

    @named port_a = TwoPhaseFluidPort()
    @named source = SmallBoundary_Ph(P_in=P_in, h_in=h_in)

    subs = [port_a; source]

    eqns = Equation[]

    push!(eqns, connect(source.port1, port_a))

    sys = ODESystem(eqns, t, [], [], name=name)
    sys = compose(sys, subs)
end

@named n1m1 = N1M1()
@named pipe = AdiabaticStraightPipe()
@named sink = MassFlowSource_h(m_flow_in=-0.01, h_in=400e3)

streams_a = [n1m1.port_a, pipe.port_a]
streams_b = [pipe.port_b, sink.port]

eqns = [
        connect(n1m1.port_a, pipe.port_a)
        connect(pipe.port_b, sink.port)
       ]

@named sys = ODESystem(eqns, t)
@named n1m1Test = compose(sys, n1m1, pipe, sink)
@test_nowarn structural_simplify(n1m1Test)


# N1M2 model and test code.
function N1M2(;name,
        P_in=1e6,
        h_in=400e3,
        kwargs...,
    )

    @named port_a = TwoPhaseFluidPort()
    @named port_b = TwoPhaseFluidPort()

    @named source = SmallBoundary_Ph(P_in=P_in, h_in=h_in)

    subs = [port_a; port_b; source]

    eqns = Equation[]

    push!(eqns, connect(source.port1, port_a))
    push!(eqns, connect(source.port1, port_b))

    sys = ODESystem(eqns, t, [], [], name=name)
    sys = compose(sys, subs)
end

@named n1m2 = N1M2()
@named sink1 = MassFlowSource_h(m_flow_in=-0.01, h_in=400e3)
@named sink2 = MassFlowSource_h(m_flow_in=-0.01, h_in=400e3)

eqns = [
        connect(n1m2.port_a, sink1.port)
        connect(n1m2.port_b, sink2.port)
       ]

@named sys = ODESystem(eqns, t)
@named n1m2Test = compose(sys, n1m2, sink1, sink2)
@test_nowarn structural_simplify(n1m2Test)


@named n1m2 = N1M2()
@named pipe1 = AdiabaticStraightPipe()
@named pipe2 = AdiabaticStraightPipe()
@named sink1 = MassFlowSource_h(m_flow_in=-0.01, h_in=400e3)
@named sink2 = MassFlowSource_h(m_flow_in=-0.01, h_in=400e3)

eqns = [
        connect(n1m2.port_a, pipe1.port_a)
        connect(pipe1.port_b, sink1.port)

        connect(n1m2.port_b, pipe2.port_a)
        connect(pipe2.port_b, sink2.port)
       ]

@named sys = ODESystem(eqns, t)
@named n1m2AltTest = compose(sys, n1m2, pipe1, pipe2, sink1, sink2)
@test_nowarn structural_simplify(n1m2AltTest)


# N2M2 model and test code.
function N2M2(;name,
        kwargs...,
    )

    @named port_a = TwoPhaseFluidPort()
    @named port_b = TwoPhaseFluidPort()
    @named pipe = AdiabaticStraightPipe()

    streams_a = [port_a, pipe.port_a]
    streams_b = [pipe.port_b, port_b]

    subs = [port_a; port_b; pipe]

    eqns = Equation[]

    push!(eqns, connect(port_a, pipe.port_a))
    push!(eqns, connect(pipe.port_b, port_b))

    sys = ODESystem(eqns, t, [], [], name=name)
    sys = compose(sys, subs)
end

@named n2m2 = N2M2()
@named source = MassFlowSource_h(m_flow_in=-0.01, h_in=400e3)
@named sink = SmallBoundary_Ph(P_in=1e6, h_in=400e3)

eqns = [
        connect(source.port, n2m2.port_a)
        connect(n2m2.port_b, sink.port1)
       ]

@named sys = ODESystem(eqns, t)
@named n2m2Test = compose(sys, n2m2, source, sink)
@test_nowarn structural_simplify(n2m2Test)

# array var
@connector function VecPin(;name)
    sts = @variables v[1:2](t)=[1.0,0.0] i[1:2](t)=1.0 [connect = Flow]
    ODESystem(Equation[], t, [sts...;], []; name=name)
end

@named vp1 = VecPin()
@named vp2 = VecPin()

@named simple = ODESystem([connect(vp1, vp2)], t)
sys = expand_connections(compose(simple, [vp1, vp2]))
@test equations(sys) == [
                         vp1.v[1] ~ vp2.v[1]
                         vp1.v[2] ~ vp2.v[2]
                         0 ~ -vp1.i[1] - vp2.i[1]
                         0 ~ -vp1.i[2] - vp2.i[2]
                        ]
