using Test
using ModelingToolkitBase, OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using ModelingToolkitBase:
                           t_nounits as t, D_nounits as D, generate_connection_set,
                           scalarize, ConnectionVertex
using Symbolics, Graphs
import SymbolicUtils as SU
import ModelingToolkitBase as MTK
import SciMLBase

@connector function Pin(; name)
    vars = @variables begin
        v(t)
        i(t), [connect = Flow]
    end
    System(Equation[], t, vars, []; name)
end

@component function Resistor(; name, R = 1.0)
    @parameters R = R
    @named p = Pin()
    @named n = Pin()
    @variables v(t) i(t)
    eqs = [v ~ p.v - n.v; i ~ p.i; p.i + n.i ~ 0; v ~ i * R]
    System(eqs, t, [v, i], [R]; systems = [p, n], name)
end

@component function Capacitor(; name, C = 1.0)
    @parameters C = C
    @named p = Pin()
    @named n = Pin()
    @variables v(t) i(t)
    eqs = [v ~ p.v - n.v; i ~ p.i; p.i + n.i ~ 0; D(v) ~ i / C]
    System(eqs, t, [v, i], [C]; systems = [p, n], name)
end

@component function Voltage(; name, V = 1.0)
    @parameters V = V
    @named p = Pin()
    @named n = Pin()
    @variables v(t)
    eqs = [v ~ p.v - n.v; v ~ V; p.i + n.i ~ 0]
    System(eqs, t, [v], [V]; systems = [p, n], name)
end

@component function Ground(; name)
    @named g = Pin()
    System([g.v ~ 0], t; systems = [g], name)
end

# node with a single connector for the `portmap = nothing` case
@component function SinglePinNode(; name)
    @named port = Pin()
    @variables x(t)
    System([D(x) ~ x], t, [x], []; systems = [port], name)
end

@connector function RealInput(; name)
    vars = @variables u(t) [input = true]
    System(Equation[], t, vars, []; name)
end

@connector function RealOutput(; name)
    vars = @variables u(t) [output = true]
    System(Equation[], t, vars, []; name)
end

@component function Gain(; name, k = 1.0)
    @parameters k = k
    @named input = RealInput()
    @named output = RealOutput()
    System([output.u ~ k * input.u], t, [], [k]; systems = [input, output], name)
end

@component function Source(; name, y = 1.0)
    @parameters y = y
    @named output = RealOutput()
    System([output.u ~ y], t, [], [y]; systems = [output], name)
end

@component function Sink(; name)
    @named input = RealInput()
    @variables u(t)
    System([D(u) ~ input.u], t, [u], []; systems = [input], name)
end

function cset_key(cvert::ConnectionVertex)
    return (join(cvert.name, "."), collect(cvert.idx), cvert.isouter, cvert.type)
end

function cset_keys(sys)
    _, (csets, _) = generate_connection_set(sys)
    return Set(Set(cset_key(v) for v in cset) for cset in csets)
end

function parallel_rc_model(multiconnect_form::Bool)
    resistor1 = Resistor(name = :resistor1)
    resistor2 = Resistor(name = :resistor2, R = 2.0)
    capacitor = Capacitor(name = :capacitor)
    source = Voltage(name = :source)
    ground = Ground(name = :ground)
    nodes = [resistor1, resistor2, capacitor, source, ground]
    if multiconnect_form
        # nets: {source.p, resistor1.p, resistor2.p}, {resistor1.n, resistor2.n,
        # capacitor.p}, {capacitor.n, source.n, ground.g}
        edges = ConnectionEdge[
            ConnectionEdge(4, 1, :p, :p),
            ConnectionEdge(4, 2, :p, :p),
            ConnectionEdge(1, 3, :n, :p),
            ConnectionEdge(2, 3, :n, :p),
            ConnectionEdge(3, 4, :n, :n),
            ConnectionEdge(3, 5, :n, :g)
        ]
        eqs = [multiconnect(nodes, edges)]
    else
        eqs = Equation[
            connect(source.p, resistor1.p),
            connect(source.p, resistor2.p),
            connect(resistor1.n, capacitor.p),
            connect(resistor2.n, capacitor.p),
            connect(capacitor.n, source.n),
            connect(capacitor.n, ground.g)
        ]
    end
    @named sys = System(eqs, t)
    return compose(sys, nodes)
end

function summands(ex)
    return iscall(ex) && operation(ex) === (+) ? SU.arguments(ex) : SU.SArgsT((ex,))
end

@testset "Subsystem network construction and validation" begin
    resistor1 = Resistor(name = :r1)
    resistor2 = Resistor(name = :r2)
    nodes = [resistor1, resistor2]
    edges = [ConnectionEdge(1, 2, :n, :p)]
    eq = multiconnect(nodes, edges)
    @test eq isa Equation
    @test MTK.value(eq.lhs) isa Connection
    net = MTK.value(eq.rhs).systems
    @test net isa ConnectionNetwork
    @test net.nodes == nodes
    @test net.edges == edges
    @test_nowarn show(IOBuffer(), eq)

    # unique node names required
    dup = [Resistor(name = :r), Resistor(name = :r)]
    @test_throws ["unique"] multiconnect(dup, edges)

    # edge indices must be in bounds
    @test_throws ["only 2 nodes"] multiconnect(nodes, [ConnectionEdge(1, 3, :n, :p)])
    @test_throws ["only 2 nodes"] multiconnect(nodes, [ConnectionEdge(0, 2, :n, :p)])

    # ports must exist
    @test_throws ["no connector or variable named `q`"] multiconnect(
        nodes, [ConnectionEdge(1, 2, :n, :q)]
    )

    # no self-connections on the same port
    @test_throws ["connects a port to itself"] multiconnect(
        nodes, [ConnectionEdge(1, 1, :n, :n)]
    )

    # connector ports cannot be paired with variable ports
    @test_throws ["different kinds"] multiconnect(
        nodes, [ConnectionEdge(1, 2, :n, :i)]
    )
end

@testset "Subsystem network expands to the same connection sets" begin
    sys_mc = parallel_rc_model(true)
    sys_man = parallel_rc_model(false)
    @test cset_keys(sys_mc) == cset_keys(sys_man)
end

@testset "Subsystem network solves identically to manual connections" begin
    sys_mc = mtkcompile(parallel_rc_model(true))
    prob_mc = ODEProblem(sys_mc, [sys_mc.capacitor.v => 0.0], (0.0, 5.0))
    sol_mc = solve(prob_mc, Rodas5P(); reltol = 1.0e-8, abstol = 1.0e-8)
    @test SciMLBase.successful_retcode(sol_mc)

    sys_man = mtkcompile(parallel_rc_model(false))
    prob_man = ODEProblem(sys_man, [sys_man.capacitor.v => 0.0], (0.0, 5.0))
    sol_man = solve(prob_man, Rodas5P(); reltol = 1.0e-8, abstol = 1.0e-8)
    @test sol_mc[sys_mc.capacitor.v] ≈ sol_man[sys_man.capacitor.v] atol = 1.0e-6
end

@testset "Graphs.jl interface with port maps" begin
    resistor1 = Resistor(name = :r1)
    capacitor = Capacitor(name = :c1)
    ground = Ground(name = :gnd)
    nodes = [resistor1, capacitor, ground]
    g = Graphs.path_graph(3)  # r1 - c1 - gnd

    # callable port map
    portmap = function (i, j)
        i == 1 && return (:n, :p)
        return (:n, :g)
    end
    eq = multiconnect(nodes, g, portmap)
    net = MTK.value(eq.rhs).systems
    @test net.edges == [
        ConnectionEdge(1, 2, :n, :p), ConnectionEdge(2, 3, :n, :g)
    ]

    # vertex count must match the number of nodes
    @test_throws ["4 vertices"] multiconnect(nodes, Graphs.path_graph(4), portmap)

    # `nothing` portmap uses the single connector of each node
    single_nodes = [
        SinglePinNode(name = :n1), SinglePinNode(name = :n2),
        SinglePinNode(name = :n3)
    ]
    eq_auto = multiconnect(single_nodes, Graphs.path_graph(3))
    @test MTK.value(eq_auto.rhs).systems.edges ==
          [ConnectionEdge(1, 2, :port, :port), ConnectionEdge(2, 3, :port, :port)]

    # ambiguous `nothing` portmap errors with the list of connectors
    @test_throws ["Cannot infer ports"] multiconnect(nodes, Graphs.path_graph(3))

    # dictionary portmap
    dictmap = Dict(Graphs.Edge(1, 2) => (:n, :p), Graphs.Edge(2, 3) => (:n, :g))
    eq_dict = multiconnect(nodes, g, dictmap)
    @test MTK.value(eq_dict.rhs).systems.edges == net.edges
    @test_throws ["No port mapping"] multiconnect(
        nodes, g, Dict(Graphs.Edge(1, 2) => (:n, :p))
    )
end

@testset "Causal variable and connector ports" begin
    source = Source(name = :src)
    gain1 = Gain(name = :g1)
    gain2 = Gain(name = :g2)
    sink = Sink(name = :sink)
    nodes = [source, gain1, gain2, sink]
    g = Graphs.path_digraph(4)
    eq = multiconnect(nodes, g, :output => :input)
    @named sys_mc = System([eq], t)
    sys_mc = compose(sys_mc, nodes)

    source2 = Source(name = :src)
    gain12 = Gain(name = :g1)
    gain22 = Gain(name = :g2)
    sink2 = Sink(name = :sink)
    eqs = Equation[
        connect(source2.output, gain12.input),
        connect(gain12.output, gain22.input),
        connect(gain22.output, sink2.input)
    ]
    @named sys_man = System(eqs, t)
    sys_man = compose(sys_man, [source2, gain12, gain22, sink2])

    @test cset_keys(sys_mc) == cset_keys(sys_man)

    # solving the causal chain: sink.u integrates 1.0 * 2.0 * 3.0
    gain1 = Gain(name = :g1, k = 2.0)
    gain2 = Gain(name = :g2, k = 3.0)
    nodes = [Source(name = :src), gain1, gain2, Sink(name = :sink)]
    eq = multiconnect(nodes, g, :output => :input)
    @named sys = System([eq], t)
    sys = mtkcompile(compose(sys, nodes))
    prob = ODEProblem(sys, [sys.sink.u => 0.0], (0.0, 1.0))
    sol = solve(prob, Rodas5P(); reltol = 1.0e-8, abstol = 1.0e-8)
    @test sol[sys.sink.u][end] ≈ 6.0 atol = 1.0e-6

    # causal connections on bare variables (not wrapped in connectors)
    blk = function (; name, gain = 1.0)
        @parameters k = gain
        vars = @variables begin
            u(t), [input = true]
            y(t), [output = true]
        end
        System([y ~ k * u], t, vars, [k]; name)
    end
    b1 = blk(name = :b1, gain = 2.0)
    b2 = blk(name = :b2, gain = 3.0)
    eq = multiconnect([b1, b2], Graphs.path_digraph(2), :y => :u)
    @named sys = System([eq], t)
    sys = compose(sys, [b1, b2])

    b1m = blk(name = :b1, gain = 2.0)
    b2m = blk(name = :b2, gain = 3.0)
    @named sys_man = System([connect(b1m.y, b2m.u)], t)
    sys_man = compose(sys_man, [b1m, b2m])
    @test cset_keys(sys) == cset_keys(sys_man)
end

@testset "Connector-as-node and nested networks" begin
    # nodes may be connectors themselves via `port = nothing`
    p1 = Pin(name = :p1)
    resistor = Resistor(name = :r)
    nodes = [p1, resistor]
    eq = multiconnect(nodes, [ConnectionEdge(1, 2, nothing, :p)])
    @named sys = System([eq], t)
    sys = compose(sys, nodes)

    p1m = Pin(name = :p1)
    rm = Resistor(name = :r)
    @named sys_man = System([connect(p1m, rm.p)], t)
    sys_man = compose(sys_man, [p1m, rm])
    @test cset_keys(sys) == cset_keys(sys_man)

    # a `multiconnect` inside a subsystem gets correct namespacing
    r1 = Resistor(name = :r1)
    r2 = Resistor(name = :r2)
    mid = System(
        [multiconnect([r1, r2], [ConnectionEdge(1, 2, :n, :p)])], t; name = :mid
    )
    mid = compose(mid, [r1, r2])
    top = System(Equation[], t; name = :top)
    top = compose(top, [mid])

    r1m = Resistor(name = :r1)
    r2m = Resistor(name = :r2)
    mid_man = System([connect(r1m.n, r2m.p)], t; name = :mid)
    mid_man = compose(mid_man, [r1m, r2m])
    top_man = System(Equation[], t; name = :top)
    top_man = compose(top_man, [mid_man])
    @test cset_keys(top) == cset_keys(top_man)
end

@testset "Array-variable form: equality and flow" begin
    @variables v(t)[1:4] f(t)[1:4]
    v, f = unwrap(v), unwrap(f)
    f = MTK.setmetadata(f, MTK.VariableConnectType, MTK.Flow)

    g = Graphs.star_graph(4)  # edges 1-2, 1-3, 1-4: one net of all 4 ports
    eqs = multiconnect((v, f) => (v, f), g)
    @test length(eqs) == 2
    eq_eq, eq_fl = eqs
    # the single compact equality equation expands to the per-edge equalities
    @test issetequal(scalarize(eq_eq), Equation[v[1] ~ v[2], v[1] ~ v[3], v[1] ~ v[4]])
    # all four ports share the net since node 1's port is common to all edges
    @test isequal(eq_fl.lhs, Symbolics.COMMON_ZERO)
    @test issetequal(summands(eq_fl.rhs), [f[i] for i in 1:4])

    # disjoint edges form separate nets
    g2 = Graphs.SimpleGraph(4)
    add_edge!(g2, 1, 2)
    add_edge!(g2, 3, 4)
    eqs2 = multiconnect(f => f, g2)
    nets = Set(Set{Any}(summands(eq.rhs)) for eq in eqs2)
    @test nets == Set([Set{Any}([f[1], f[2]]), Set{Any}([f[3], f[4]])])

    # different variables on each side are distinct ports
    @variables fa(t)[1:4] fb(t)[1:4]
    fa, fb = MTK.setmetadata.(unwrap.([fa, fb]), MTK.VariableConnectType, MTK.Flow)
    g3 = Graphs.path_graph(4)
    eqs3 = multiconnect(fa => fb, g3)
    # nets: {fa[1],fb[2]}, {fa[2],fb[3]}, {fa[3],fb[4]}, {fa[4]}, {fb[1]}
    @test length(eqs3) == 5
    @test any(eq -> isequal(eq.rhs, fb[1]), eqs3)
    @test any(eq -> isequal(eq.rhs, fa[4]), eqs3)
end

@testset "Array-variable form: causal" begin
    @variables o(t)[1:4] u(t)[1:4] w(t)[1:4]
    o, u, w = unwrap.([o, u, w])
    o = MTK.setmetadata(o, MTK.VariableOutput, true)
    u = MTK.setmetadata(u, MTK.VariableInput, true)

    dg = Graphs.path_digraph(4)
    eqs = multiconnect(o => u, dg)
    @test only(eqs) isa Equation
    @test issetequal(scalarize(only(eqs)), Equation[o[1] ~ u[2], o[2] ~ u[3], o[3] ~ u[4]])

    # an input may not be driven by two outputs
    bad = Graphs.DiGraph(4)
    add_edge!(bad, 1, 3)
    add_edge!(bad, 2, 3)
    @test_throws ["driven by multiple edges"] multiconnect(o => u, bad)

    # output-input orientation is validated
    @test_throws ["pair an output variable with an input"] multiconnect(o => o, dg)
    @test_throws ["pair an output variable with an input"] multiconnect(u => u, dg)

    # a port variable without causal metadata cannot pair with an io-typed one
    @test_throws ["pair an output variable with an input"] multiconnect(o => w, dg)

    # port variables must be 1-D arrays of length nv(g)
    @variables x(t) arr2(t)[1:3]
    x, arr2 = unwrap.([x, arr2])
    @test_throws ["one-dimensional symbolic array"] multiconnect(x => x, dg)
    @test_throws ["expected `nv(graph) = 4`"] multiconnect(arr2 => arr2, dg)
end

@testset "Connection networks survive system transformations" begin
    sys_mc = parallel_rc_model(true)
    newiv = sys_mc.capacitor.v
    sys2 = MTK.change_independent_variable(sys_mc, newiv)
    conneq = only(
        filter(e -> MTK.value(e.lhs) isa Connection, MTK.get_eqs(sys2))
    )
    net = MTK.value(conneq.rhs).systems
    @test net isa ConnectionNetwork
    @test nameof.(net.nodes) ==
          [:resistor1, :resistor2, :capacitor, :source, :ground]
    # the network must reference the transformed node systems, not the stale
    # pre-transformation objects
    @test all(n -> isequal(MTK.get_iv(n), MTK.get_iv(sys2)), net.nodes)

    # `toexpr` emits a reconstructible `multiconnect` expression
    io = IOBuffer()
    write(io, sys_mc)
    str = String(take!(io))
    @test occursin("multiconnect)([resistor1, resistor2, capacitor, source, ground]", str)
    @test occursin("ConnectionEdge)(4, 1, :p, :p)", str)
end

@testset "Array-variable network through mtkcompile" begin
    # A series RC loop written entirely in array variables. Node 1 is an ideal
    # voltage source tied to ground, nodes 2 and 3 are resistors, and node 4 is a
    # capacitor. Every node has a p-port and an n-port; the graph closes the loop
    # n_i -> p_{i+1} (cyclically), so `multiconnect` generates all four junction
    # equalities and conservation equations from two compact symbolic equations.
    @parameters V=12.0 R=1.0 C=1.0
    @variables vp(t)[1:4] vn(t)[1:4] fp(t)[1:4] fn(t)[1:4] vc(t)
    vp, vn, fp, fn = unwrap.([vp, vn, fp, fn])
    fp = MTK.setmetadata(fp, MTK.VariableConnectType, MTK.Flow)
    fn = MTK.setmetadata(fn, MTK.VariableConnectType, MTK.Flow)
    eqs = Equation[
        vp[1] ~ 0.0, vn[1] ~ V,                       # source + ground reference
        vp[2] - vn[2] ~ R * fp[2], fp[2] + fn[2] ~ 0.0, # resistor nodes
        vp[3] - vn[3] ~ R * fp[3], fp[3] + fn[3] ~ 0.0,
        vc ~ vp[4] - vn[4], D(vc) ~ fp[4] / C, fp[4] + fn[4] ~ 0.0, # capacitor node
        multiconnect((vn, fn) => (vp, fp), Graphs.cycle_digraph(4))...
    ]
    @named sys = System(
        eqs, t,
        [scalarize(vp); scalarize(vn); scalarize(fp); scalarize(fn); vc],
        [V, R, C]
    )
    sys = mtkcompile(sys)
    prob = ODEProblem(sys, [sys.vc => 0.0], (0.0, 30.0))
    sol = solve(prob, Rodas5P(); reltol = 1.0e-8, abstol = 1.0e-8)
    @test SciMLBase.successful_retcode(sol)
    # the capacitor charges as vc(t) = V * (1 - exp(-t / (2RC))) toward V
    @test sol[sys.vc] ≈ 12.0 .* (1.0 .- exp.(.-sol.t ./ 2.0)) atol = 1.0e-4
    @test sol[sys.vn[4]][end] ≈ 0.0 atol = 1.0e-4
    @test isapprox(sol[sys.fp[2]][1], 6.0, atol = 1.0e-3) # I(0) = V / (2R)
end
