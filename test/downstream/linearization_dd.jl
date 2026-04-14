## Test that dummy_derivatives can be set to zero
# The call to Link(; m = 0.2, l = 10, I = 1, g = -9.807) hangs forever on Julia v1.6
using LinearAlgebra
using ModelingToolkit
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Mechanical.MultiBody2D
using ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition
using Test
import Symbolics
import NonlinearSolve
using Setfield: @set

using ControlSystemsMTK
using ControlSystemsMTK.ControlSystemsBase: sminreal, minreal, poles
connect = ModelingToolkit.connect

function rm_bindings(sys)
    @set sys.bindings = empty(bindings(sys))
end

@independent_variables t
D = Differential(t)

@named link1 = Link(; m = 0.2, l = 10, I = 1, g = -9.807)
link1 = rm_bindings(link1)
@named cart = TranslationalPosition.Mass(; m = 1)
@named fixed = Fixed()
@named force = Force(use_support = false)

eqs = [
    connect(link1.TX1, cart.flange)
    connect(cart.flange, force.flange)
    connect(link1.TY1, fixed.flange)
]

@named model = System(eqs, t, [], []; systems = [link1, cart, force, fixed])
lin_outputs = [cart.s, cart.v, link1.A, link1.dA]
lin_inputs = [force.f.u]

# => nothing to remove extra defaults
op = Dict(
    cart.v => 0, link1.A => -pi / 2, link1.dA => 0, force.f.u => 0,
    link1.x1 => nothing, link1.y1 => nothing, link1.x2 => nothing, link1.x_cm => nothing
)
guesses = [link1.fx1 => 0]
@info "named_ss"
G = named_ss(
    model, lin_inputs, lin_outputs; allow_symbolic = true, op,
    allow_input_derivatives = true, zero_dummy_der = false, guesses,
    balance = true,
)
G = sminreal(G)
@info "minreal"
G = minreal(G)
@info "poles"
ps = poles(G)

@test minimum(abs, ps) < 1.0e-6
@test minimum(abs, complex(0, 1.3777260367206716) .- ps) < 1.0e-10

lin_fun, ssys = ModelingToolkit.linearization_function(
    model, lin_inputs, lin_outputs, allow_symbolic = true, op = op,
    zero_dummy_der = false, guesses = guesses
);
lsys, = ModelingToolkit.linearize(ssys, lin_fun; op, allow_input_derivatives = true);
lsyss,
    sysss = ModelingToolkit.linearize_symbolic(
    model, lin_inputs, lin_outputs;
    allow_input_derivatives = true, allow_symbolic = true,
)

lsyss2 = (;
    A = lin_fun.prob[lsyss.A], B = lin_fun.prob[lsyss.B],
    C = lin_fun.prob[lsyss.C], D = lin_fun.prob[lsyss.D],
)
@test lsyss2.A ≈ lsys.A
@test lsyss2.B ≈ lsys.B
@test lsyss2.C ≈ lsys.C
@test lsyss2.D ≈ lsys.D
