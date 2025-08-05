## Test that dummy_derivatives can be set to zero
# The call to Link(; m = 0.2, l = 10, I = 1, g = -9.807) hangs forever on Julia v1.6
using LinearAlgebra
using ModelingToolkit
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Mechanical.MultiBody2D
using ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition
using Test
import NonlinearSolve

using ControlSystemsMTK
using ControlSystemsMTK.ControlSystemsBase: sminreal, minreal, poles
connect = ModelingToolkit.connect

@independent_variables t
D = Differential(t)

@named link1 = Link(; m = 0.2, l = 10, I = 1, g = -9.807)
@named cart = TranslationalPosition.Mass(; m = 1, s = 0)
@named fixed = Fixed()
@named force = Force(use_support = false)

eqs = [connect(link1.TX1, cart.flange)
       connect(cart.flange, force.flange)
       connect(link1.TY1, fixed.flange)]

@named model = System(eqs, t, [], []; systems = [link1, cart, force, fixed])
lin_outputs = [cart.s, cart.v, link1.A, link1.dA]
lin_inputs = [force.f.u]

# => nothing to remove extra defaults
op = Dict(cart.s => 10, cart.v => 0, link1.A => -pi / 2, link1.dA => 0, force.f.u => 0,
    link1.x1 => nothing, link1.y1 => nothing, link1.x2 => nothing, link1.x_cm => nothing)
guesses = [link1.fx1 => 0]
@info "named_ss"
G = named_ss(model, lin_inputs, lin_outputs; allow_symbolic = true, op,
    allow_input_derivatives = true, zero_dummy_der = true, guesses)
G = sminreal(G)
@info "minreal"
G = minreal(G)
@info "poles"
ps = poles(G)

@test minimum(abs, ps) < 1e-6
@test minimum(abs, complex(0, 1.3777260367206716) .- ps) < 1e-10

lsys,
syss = linearize(model, lin_inputs, lin_outputs, allow_symbolic = true, op = op,
    allow_input_derivatives = true, zero_dummy_der = true, guesses = guesses)
lsyss,
sysss = ModelingToolkit.linearize_symbolic(model, lin_inputs, lin_outputs;
    allow_input_derivatives = true)

dummyder = setdiff(unknowns(sysss), unknowns(model))
# op2 = merge(ModelingToolkit.guesses(model), op, Dict(x => 0.0 for x in dummyder))
op2 = merge(ModelingToolkit.defaults(syss), op)
op2[link1.fy1] = -op2[link1.g] * op2[link1.m]
op2[cart.f] = 0

@test substitute(lsyss.A, op2) ≈ lsys.A
# We cannot pivot symbolically, so the part where a linear solve is required
# is not reliable.
@test substitute(lsyss.B, op2)[1:6, 1] ≈ lsys.B[1:6, 1]
@test substitute(lsyss.C, op2) ≈ lsys.C
@test substitute(lsyss.D, op2) ≈ lsys.D
