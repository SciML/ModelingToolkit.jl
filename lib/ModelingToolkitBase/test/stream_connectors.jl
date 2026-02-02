using Test
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using OrdinaryDiffEq

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

    System(Equation[], t, vars, pars; name = name)
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

    System(eqs, t, vars, pars; name)
end

function MassFlowSource_h(;
        name,
        h_in = 420.0e3,
        m_flow_in = -0.01
    )
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

    return compose(System(eqns, t, vars, pars; name = name), subs)
end

# Simplified components.
function AdiabaticStraightPipe(;
        name,
        kwargs...
    )
    vars = []
    pars = []

    @named port_a = TwoPhaseFluidPort()
    @named port_b = TwoPhaseFluidPort()

    subs = [port_a; port_b]

    eqns = Equation[]

    push!(eqns, connect(port_a, port_b))
    sys = System(eqns, t, vars, pars; name = name)
    return sys = compose(sys, subs)
end

function SmallBoundary_Ph(;
        name,
        P_in = 1.0e6,
        h_in = 400.0e3
    )
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

    return compose(System(eqns, t, vars, pars; name = name), subs)
end

# N1M1 model and test code.
function N1M1(;
        name,
        P_in = 1.0e6,
        h_in = 400.0e3,
        kwargs...
    )
    @named port_a = TwoPhaseFluidPort()
    @named source = SmallBoundary_Ph(P_in = P_in, h_in = h_in)

    subs = [port_a; source]

    eqns = Equation[]

    push!(eqns, connect(source.port1, port_a))

    sys = System(eqns, t, [], [], name = name)
    return sys = compose(sys, subs)
end

@named fluid = TwoPhaseFluid(; R = 876, B = 1.2e9, V = 0.034)
@named n1m1 = N1M1()
@named pipe = AdiabaticStraightPipe()
@named sink = MassFlowSource_h(m_flow_in = -0.01, h_in = 400.0e3)

eqns = [
    connect(n1m1.port_a, pipe.port_a)
    connect(pipe.port_b, sink.port)
]

@named sys = System(eqns, t; systems = [n1m1, pipe, sink])

eqns = [
    domain_connect(fluid, n1m1.port_a)
    connect(n1m1.port_a, pipe.port_a)
    connect(pipe.port_b, sink.port)
]

@named n1m1Test = System(eqns, t, [], []; systems = [fluid, n1m1, pipe, sink])

@test_nowarn mtkcompile(n1m1Test)
@unpack source, port_a = n1m1
ssort(eqs) = sort(eqs, by = string)
@test ssort(equations(expand_connections(n1m1))) == ssort(
    [
        0 ~ port_a.m_flow
        0 ~ source.port1.m_flow - port_a.m_flow
        source.port1.P ~ port_a.P
        source.port1.P ~ source.P
        port_a.h_outflow ~ source.port1.h_outflow
        source.port1.h_outflow ~ source.h
    ]
)
@unpack port_a, port_b = pipe
@test ssort(equations(expand_connections(pipe))) ==
    ssort(
    [
        0 ~ -port_a.m_flow - port_b.m_flow
        0 ~ port_a.m_flow
        0 ~ port_b.m_flow
        port_a.P ~ port_b.P
        port_a.h_outflow ~ instream(port_b.h_outflow)
        port_b.h_outflow ~ instream(port_a.h_outflow)
    ]
)
@test equations(expand_connections(sys)) ⊇
    [
    0 ~ n1m1.port_a.m_flow + pipe.port_a.m_flow
    0 ~ pipe.port_b.m_flow + sink.port.m_flow
    n1m1.port_a.P ~ pipe.port_a.P
    pipe.port_b.P ~ sink.port.P
]
@test ssort(equations(expand_connections(n1m1Test))) ==
    ssort(
    [
        0 ~ -pipe.port_a.m_flow - pipe.port_b.m_flow
        0 ~ n1m1.source.port1.m_flow - n1m1.port_a.m_flow
        0 ~ n1m1.port_a.m_flow + pipe.port_a.m_flow
        0 ~ pipe.port_b.m_flow + sink.port.m_flow
        fluid.m_flow ~ 0
        n1m1.port_a.P ~ pipe.port_a.P
        n1m1.source.port1.P ~ n1m1.port_a.P
        n1m1.source.port1.P ~ n1m1.source.P
        n1m1.port_a.h_outflow ~ n1m1.source.port1.h_outflow
        n1m1.source.port1.h_outflow ~ n1m1.source.h
        pipe.port_a.P ~ pipe.port_b.P
        pipe.port_a.h_outflow ~ sink.port.h_outflow
        pipe.port_b.P ~ sink.port.P
        pipe.port_b.h_outflow ~ n1m1.port_a.h_outflow
        sink.port.P ~ sink.P
        sink.port.h_outflow ~ sink.h_in
        sink.port.m_flow ~ -sink.m_flow_in
    ]
)

# N1M2 model and test code.
function N1M2(;
        name,
        P_in = 1.0e6,
        h_in = 400.0e3,
        kwargs...
    )
    @named port_a = TwoPhaseFluidPort()
    @named port_b = TwoPhaseFluidPort()

    @named source = SmallBoundary_Ph(P_in = P_in, h_in = h_in)

    subs = [port_a; port_b; source]

    eqns = Equation[]

    push!(eqns, connect(source.port1, port_a))
    push!(eqns, connect(source.port1, port_b))

    sys = System(eqns, t, [], [], name = name)
    return sys = compose(sys, subs)
end

@named n1m2 = N1M2()
@named sink1 = MassFlowSource_h(m_flow_in = -0.01, h_in = 400.0e3)
@named sink2 = MassFlowSource_h(m_flow_in = -0.01, h_in = 400.0e3)

eqns = [
    connect(n1m2.port_a, sink1.port)
    connect(n1m2.port_b, sink2.port)
]

@named sys = System(eqns, t)
@named n1m2Test = compose(sys, n1m2, sink1, sink2)
@test_nowarn mtkcompile(n1m2Test)

@named n1m2 = N1M2()
@named pipe1 = AdiabaticStraightPipe()
@named pipe2 = AdiabaticStraightPipe()
@named sink1 = MassFlowSource_h(m_flow_in = -0.01, h_in = 400.0e3)
@named sink2 = MassFlowSource_h(m_flow_in = -0.01, h_in = 400.0e3)

eqns = [
    connect(n1m2.port_a, pipe1.port_a)
    connect(pipe1.port_b, sink1.port)
    connect(n1m2.port_b, pipe2.port_a)
    connect(pipe2.port_b, sink2.port)
]

@named sys = System(eqns, t)
@named n1m2AltTest = compose(sys, n1m2, pipe1, pipe2, sink1, sink2)
@test_nowarn mtkcompile(n1m2AltTest)

# N2M2 model and test code.
function N2M2(;
        name,
        kwargs...
    )
    @named port_a = TwoPhaseFluidPort()
    @named port_b = TwoPhaseFluidPort()
    @named pipe = AdiabaticStraightPipe()

    subs = [port_a; port_b; pipe]

    eqns = Equation[]

    push!(eqns, connect(port_a, pipe.port_a))
    push!(eqns, connect(pipe.port_b, port_b))

    sys = System(eqns, t, [], [], name = name)
    return sys = compose(sys, subs)
end

@named n2m2 = N2M2()
@named source = MassFlowSource_h(m_flow_in = -0.01, h_in = 400.0e3)
@named sink = SmallBoundary_Ph(P_in = 1.0e6, h_in = 400.0e3)

eqns = [
    connect(source.port, n2m2.port_a)
    connect(n2m2.port_b, sink.port1)
]

@named sys = System(eqns, t)
@named n2m2Test = compose(sys, n2m2, source, sink)
@test_nowarn mtkcompile(n2m2Test)

# stream var
@named sp1 = TwoPhaseFluidPort()
@named sp2 = TwoPhaseFluidPort()
@named sys = System([connect(sp1, sp2)], t)
sys_exp = expand_connections(compose(sys, [sp1, sp2]))
@test ssort(equations(sys_exp)) == ssort(
    [
        0 ~ -sp1.m_flow - sp2.m_flow
        0 ~ sp1.m_flow
        0 ~ sp2.m_flow
        sp1.P ~ sp2.P
        sp1.h_outflow ~ ModelingToolkitBase.instream(sp2.h_outflow)
        sp2.h_outflow ~ ModelingToolkitBase.instream(sp1.h_outflow)
    ]
)

# array var
@connector function VecPin(; name)
    sts = @variables v(t)[1:2] = [1.0, 0.0] i(t)[1:2] = ones(2) [connect = Flow]
    System(Equation[], t, [sts...;], []; name = name)
end

@named vp1 = VecPin()
@named vp2 = VecPin()
@named vp3 = VecPin()

@named simple = System([connect(vp1, vp2, vp3)], t)
sys = expand_connections(compose(simple, [vp1, vp2, vp3]))
@test ssort(equations(sys)) == ssort(
    [
        0 .~ collect(vp1.i)
        0 .~ collect(vp2.i)
        0 .~ collect(vp3.i)
        collect(vp1.v) .~ collect(vp2.v)
        collect(vp1.v) .~ collect(vp3.v)
        0 ~ -vp1.i[1] - vp2.i[1] - vp3.i[1]
        0 ~ -vp1.i[2] - vp2.i[2] - vp3.i[2]
    ]
)

@connector function VectorHeatPort(; name, N = 100, T0 = 0.0, Q0 = 0.0)
    @variables (T(t))[1:N] = T0 * ones(N) (Q(t))[1:N] = Q0 * ones(N) [connect = Flow]
    System(Equation[], t, [T; Q], []; name = name)
end

@test_nowarn @named a = VectorHeatPort()

# --------------------------------------------------
# Test the new Domain feature

sys_ = expand_connections(n1m1Test)
sys_defs = ModelingToolkitBase.bindings(sys_)
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

    System(eqs, t, vars, pars; name, initial_conditions = [dm => 0])
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
        dm ~ 0,
    ]

    System(eqs, t, vars, pars; name)
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
        H.p ~ p_int * (t > 0.01),
    ]

    return System(eqs, t, vars, pars; name, systems)
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
    eqs = [
        D(vrho) ~ drho
        vrho ~ rho_0 * (1 + p / H.bulk)
        H.p ~ p
        H.dm ~ drho * V
    ]

    return System(
        eqs, t, vars, pars; name, systems,
        initial_conditions = [vrho => rho_0 * (1 + p_int / H.bulk)]
    )
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
    eqs = [
        HA.p - HB.p ~ HA.dm * resistance / HA.viscosity
        0 ~ HA.dm + HB.dm
        domain_connect(HA, HB)
    ]

    return System(eqs, t, vars, pars; name, systems)
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

    eqs = [
        connect(v1.H, p12.HA, HA)
        connect(v2.H, p12.HB, HB)
    ]

    return System(eqs, t, vars, pars; name, systems)
end

function TwoFluidSystem(; name)
    pars = []
    vars = []

    # nodes -------------------------------
    systems = @named begin
        fluid_a = Fluid(; R = 876, B = 1.2e9, V = 0.034)
        source_a = StepSource(; P = 10.0e5)
        pipe_a = Pipe(; P = 0, R = 1.0e6)
        volume_a = StaticVolume(; P = 0, V = 0.1)

        fluid_b = Fluid(; R = 1000, B = 2.5e9, V = 0.00034)
        source_b = StepSource(; P = 10.0e5)
        pipe_b = Pipe(; P = 0, R = 1.0e6)
        volume_b = StaticVolume(; P = 0, V = 0.1)
    end

    # equations ---------------------------
    eqs = [
        connect(fluid_a, source_a.H)
        connect(source_a.H, pipe_a.HA)
        connect(pipe_a.HB, volume_a.H)
        connect(fluid_b, source_b.H)
        connect(source_b.H, pipe_b.HA)
        connect(pipe_b.HB, volume_b.H)
    ]

    return System(eqs, t, vars, pars; name, systems)
end

@named two_fluid_system = TwoFluidSystem()
sys = expand_connections(two_fluid_system)

sys_defs = ModelingToolkitBase.bindings(sys)
csys = complete(two_fluid_system)

@test Symbol(sys_defs[csys.volume_a.H.rho]) == Symbol(csys.fluid_a.rho)
@test Symbol(sys_defs[csys.volume_b.H.rho]) == Symbol(csys.fluid_b.rho)

@test_nowarn mtkcompile(two_fluid_system)

function OneFluidSystem(; name)
    pars = []
    vars = []

    # nodes -------------------------------
    systems = @named begin
        fluid = Fluid(; R = 876, B = 1.2e9, V = 0.034)

        source_a = StepSource(; P = 10.0e5)
        pipe_a = Pipe(; P = 0, R = 1.0e6)
        volume_a = StaticVolume(; P = 0, V = 0.1)

        source_b = StepSource(; P = 20.0e5)
        pipe_b = Pipe(; P = 0, R = 1.0e6)
        volume_b = StaticVolume(; P = 0, V = 0.1)
    end

    # equations ---------------------------
    eqs = [
        connect(fluid, source_a.H, source_b.H)
        connect(source_a.H, pipe_a.HA)
        connect(pipe_a.HB, volume_a.H)
        connect(source_b.H, pipe_b.HA)
        connect(pipe_b.HB, volume_b.H)
    ]

    return System(eqs, t, vars, pars; name, systems)
end

@named one_fluid_system = OneFluidSystem()
sys = expand_connections(one_fluid_system)

sys_defs = ModelingToolkitBase.bindings(sys)
csys = complete(one_fluid_system)

@test Symbol(sys_defs[csys.volume_a.H.rho]) == Symbol(csys.fluid.rho)
@test Symbol(sys_defs[csys.volume_b.H.rho]) == Symbol(csys.fluid.rho)

if @isdefined(ModelingToolkit)
    @test_nowarn mtkcompile(one_fluid_system)
end

@testset "Issue#4219: `instream` of array variables handled with `instream_rt`" begin
    @connector function FluidPort(; name, Ns = 1)
        vars = @variables begin
            p(t), [guess = 0.0, description = "Pressure, Pa"]
            md(t), [connect = Flow, guess = 0.0, description = "Mass flow, kg/s"]
            x(t)[1:Ns], [connect = Stream, guess = 1 / Ns, description = "Mass fractions, -"]
        end
        System(Equation[], t, vars, []; name = name)
    end

    @component function PressureBoundary(; name, Ns = 1, p_set = 101325.0)
        @named a = FluidPort(Ns = Ns)
        x_p = ones(Ns) ./ Ns
        pars = @parameters begin
            p = p_set
            x[1:Ns] = x_p
        end
        eqs = [
            a.p ~ p
            a.x ~ x
        ]
        compose(System(eqs, t, [], pars; name = name), a)
    end

    @component function PressureReference(; name, Ns = 1, p_set = 101325.0)
        @named a = FluidPort(Ns = Ns)
        x_p = ones(Ns) ./ Ns
        pars = @parameters begin
            p = p_set
            x[1:Ns] = x_p
        end
        eqs = [
            a.p ~ p
            a.md ~ 0.0
            a.x ~ x
        ]
        compose(System(eqs, t, [], pars; name = name), a)
    end

    @component function PumpFixed(; name, md_set = 1.0, x_set = [1])
        Ns = length(x_set)
        @named a = FluidPort(Ns = Ns)
        pars = @parameters begin
            md = md_set, [description = "Specified mass flow rate"]
            x[1:Ns] = x_set, [description = "Specified mass fractions"]
        end
        eqs = [
            # Inflow to system means outflow from component (negative md)
            a.md ~ -md,
            # Define the stream variables for the port
            a.x ~ x,
        ]
        compose(System(eqs, t, [], pars; name = name), a)
    end

    @component function Valve(; name, Ns = 1, Kv = 1.0)
        @named a = FluidPort(Ns = Ns)
        @named b = FluidPort(Ns = Ns)
        pars = @parameters begin
            # Units must be (mass flow) / (pressure^0.5)
            K = Kv, [description = "Flow coefficient"]
            eps = 1.0e-4, [description = "Regularization factor"]
        end
        vars = @variables begin
            Dp(t), [description = "Pressure drop a-b"]
        end
        eqs = [
            # 1. Pressure drop
            Dp ~ a.p - b.p
            # 2. Valve quasi-turbulent flow
            a.md ~ K * Dp / (Dp^2 + eps)^(1 / 4)
            # 3. Continuity: No storage
            a.md + b.md ~ 0
            # 4. Stream Propagation: Same composition on both sides
            a.x ~ instream(b.x)
            b.x ~ instream(a.x)
        ]
        compose(System(eqs, t, vars, pars; name = name), a, b)
    end

    @component function Tank(; name, Ap = 1.0, rho_p = 1000.0, m0 = [1])
        Ns = length(m0)
        m0 = abs.(m0) .+ 1.0e-9 # Prevent zero initial mass
        @named a = FluidPort(Ns = Ns) # Top Port (Headspace)
        @named b = FluidPort(Ns = Ns) # Bottom Port (Drain)
        gp = 9.81
        pars = @parameters begin
            A = Ap
            rho = rho_p
            g = gp
        end
        # Symbolic variables
        vars = @variables begin
            m(t)[1:Ns] = m0
            h(t)
            x(t)[1:Ns]
        end
        xa = [ifelse(a.md > 0, instream(a.x[i]), x[i]) for i in 1:Ns]
        xb = [ifelse(b.md > 0, instream(b.x[i]), x[i]) for i in 1:Ns]
        eqs = [
            # --- Geometry ---
            h ~ sum(m) / (rho * A)
            # --- Pressure
            b.p ~ a.p + rho * g * h
            # --- Mass fractions ---
            x ~ m ./ sum(m)
            # --- Streams ---
            a.x ~ x
            b.x ~ x
            # --- Dynamics ---
            D(m) ~ a.md * xa + b.md * xb
        ]
        compose(System(eqs, t, vars, pars; name = name), a, b)
    end
    m0_1 = [2, 1]
    Ns = length(m0_1)
    md_in_1 = 2.0
    x_in_1 = [0.8, 0.2]
    p_a = 101325.0
    A_1 = 0.3
    rho = 1000.0
    K_12 = 0.1
    @named mp_1 = PumpFixed(; md_set = md_in_1, x_set = x_in_1)
    @named pr_1 = PressureReference(; Ns = Ns, p_set = p_a)
    @named pb_1 = PressureBoundary(; Ns = Ns, p_set = p_a)
    @named va_1e = Valve(; Ns = Ns, Kv = K_12)
    @named tnk_1 = Tank(; Ap = A_1, rho_p = rho, m0 = m0_1)

    eqs = [
        connect(mp_1.a, pr_1.a, tnk_1.a)
        connect(tnk_1.b, va_1e.a)
        connect(va_1e.b, pb_1.a)
    ]

    @named __sys = System(eqs, t)
    @named _sys = compose(__sys, [mp_1, pr_1, pb_1, va_1e, tnk_1])
    @test_nowarn expand_connections(_sys)

    if @isdefined(ModelingToolkit)
        sys = mtkcompile(_sys)
        prob = ODEProblem(sys, nothing, (0.0, 30.0))
        sol = solve(prob, Tsit5())
        @test SciMLBase.successful_retcode(sol)
    end
end
