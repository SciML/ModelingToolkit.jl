using ModelingToolkit, ModelingToolkitStandardLibrary, Test
using OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks

function Segment(; name)
    @named R = Resistor(; R = 1)
    @named r = Resistor(; R = 1)
    @named C = Capacitor(; C = 1)

    @named p1 = Pin() # top-left
    @named p2 = Pin() # top-right
    @named n = Pin() # bottom

    eqs = [connect(p1, R.p)
           connect(R.n, p2, r.p)
           connect(r.n, C.p)
           connect(C.n, n)]

    return ODESystem(eqs, t, [], [];
                     name = name,
                     systems = [r, R, C, n, p1, p2])
end

function Strip(; name)
    num_segments = 10
    # construct `num_segments` segments
    segments = [Segment(; name = Symbol(:St_, seg))
                for seg in 1:num_segments]

    @named p1 = Pin() # top-left
    @named p2 = Pin() # top-right
    @named n = Pin() # bottom

    eqs = [connect(p1, segments[1].p1)
           connect(p2, segments[end].p2)
           [connect(n, seg.n) for seg in segments]...
           [connect(segments[i].p2, segments[i + 1].p1) for i in 1:(num_segments - 1)]...]

    return ODESystem(eqs, t, [], []; name,
                     systems = [p1, p2, n, segments...])
end

@variables t
@named source = Voltage()
@named c = Constant(k = 0.01)

@named ground = Ground()
@named strip = Strip()

rc_eqs = [connect(c.output, source.V)
          connect(source.p, strip.p1, strip.p2)
          connect(strip.n, source.n, ground.g)]

@named rc_model = ODESystem(rc_eqs, t, systems = [strip, c, source, ground])
sys = structural_simplify(rc_model)

@show ModelingToolkit.observed(sys)

ref_prob = ODEProblem(sys, [], (0, 10))
prob = ODAEProblem(sys, [], (0, 10))

ref_sol = solve(ref_prob, Tsit5())
@test_nowarn sol = solve(prob, Tsit5())



# test that the observed variables are correct
@test sol[strip₊St_1₊r₊n₊v] ≈ sol[strip₊St_1₊r₊n₊v]
