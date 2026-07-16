# [Iterate Limiting: the `limited` operator](@id limiting)

Nonlinear solves over device equations with exponential I-V characteristics — diodes,
BJTs, MOSFETs — are the classic hard case for Newton's method: a volt-sized overshoot of
a junction voltage puts `exp(v/Vt)` astronomically far from the linearization. SPICE-family
simulators handle this with *limiting*: each Newton update of a sensitive quantity is
clipped to a trusted move relative to its previous value. The Predictor/Corrector
Newton-Raphson (PCNR) method of Aadithya, Keiter & Mei recasts limiting consistently by
making each limited quantity an explicit unknown, applying the limiter as a corrector
between Newton steps, and always evaluating residuals at the corrected iterates.

ModelingToolkit exposes this symbolically through the [`limited`](@ref) operator, in the
spirit of the Modelica-style [`homotopy`](@ref homotopy) operator: a component author
annotates the model once, and every nonlinear solve built from the model gets
predictor/corrector limiting automatically.

```julia
limited(actual, limiter)
```

  - `actual` — the expression being limited (e.g. a junction voltage).
  - `limiter` — the correction rule, written in terms of the reserved placeholders
    [`limitnew`](@ref) (the proposed value) and [`limitold`](@ref) (the previously
    accepted value), plus parameters.

## Example: a diode component with `pnjlim`

Register the SPICE3 junction limiting function as an opaque symbolic function (keeping
Julia's short-circuit branch semantics), and annotate the diode's junction voltage:

```@example limiting
using ModelingToolkit, NonlinearSolve

function pnjlim(vnew, vold, vt, vcrit)
    if vnew > vcrit && abs(vnew - vold) > 2vt
        if vold > 0
            arg = 1 + (vnew - vold) / vt
            vnew = arg > 0 ? vold + vt * log(arg) : vcrit
        else
            vnew = vt * log(vnew / vt)
        end
    end
    return vnew
end
@register_symbolic pnjlim(vnew, vold, vt, vcrit)

function DCDiode(; name, Is = 1.0e-14, Vt = 0.025)
    @variables v i
    ps = @parameters begin
        (Is::Float64 = Is)
        (Vt::Float64 = Vt)
        (vcrit::Float64 = Vt * log(Vt / (sqrt(2) * Is)))
    end
    eqs = [i ~ Is * (exp(limited(v, pnjlim(limitnew, limitold, Vt, vcrit)) / Vt) - 1)]
    return System(eqs, [v, i], ps; name)
end

function DCResistor(; name, R = 1.0e3)
    @variables v i
    @parameters R = R
    return System([v ~ i * R], [v, i], [R]; name)
end

@named diode = DCDiode()
@named res = DCResistor()
@parameters Vs = 5.0
connections = [res.i ~ diode.i, Vs ~ res.v + diode.v]
@named circuit = System(connections, [], [Vs]; systems = [diode, res])
csys = mtkcompile(circuit)

prob = NonlinearProblem(csys, [diode.v => 0.0, res.i => 0.0])
sol = solve(prob, NewtonRaphson())
sol[diode.v], sol[diode.i], sol.stats.nsteps
```

Without the annotation, plain Newton needs a couple hundred millivolt-creep iterations on
this circuit; with it, the solve converges in about a dozen. Nothing about the *usage*
changed — the limiting behavior travels with the component definition.

## What `mtkcompile` does with `limited`

For **time-independent** systems, each unique `limited(actual, limiter)` node is lowered
in the PCNR augmented form:

 1. an auxiliary irreducible unknown `limited_k` is introduced and the node is replaced
    by it (irreducible, so structural simplification keeps the limited quantity as the
    surviving representative of its alias class — the reduction of the augmented system
    happens symbolically);
 2. the consistency equation `limited_k ~ actual` is appended, and `limited_k` receives
    the symbolic guess `actual`;
 3. the limiters are compiled into a `postcondition` corrector hook on the generated
    `SciMLBase.NonlinearFunction`, which NonlinearSolve.jl's native solvers apply to every
    accepted iterate before evaluating the residual there.

Solving therefore requires a solver that supports `postcondition` (e.g. `NewtonRaphson`
and the other native NonlinearSolve.jl methods); unsupported solvers throw instead of
silently ignoring the limiter.

For **time-dependent** systems the operator is stripped to `actual` during `mtkcompile`,
so the same component library compiles unchanged for transient simulation. (Limiting
inside the nonlinear solves of implicit time steppers is a separate, future integration.)

## Contracts

  - `limiter` must satisfy `limiter == limitnew` when `limitnew == limitold`, so
    solutions are fixed points of the correction.
  - `limiter` may reference `limitnew`, `limitold`, and parameters (including bound
    parameters); referencing other unknowns is an error.
  - `limited` nodes may not be nested.
  - Jacobians treat `limited(actual, limiter)` as `actual`: the limiter is a corrector
    between iterations, not part of the residual.

## Docstrings

```@docs
ModelingToolkitBase.limited
ModelingToolkitBase.limitnew
ModelingToolkitBase.limitold
```

**Reference:** K. V. Aadithya, E. R. Keiter, T. Mei, *Predictor/Corrector Newton-Raphson
(PCNR): A Simple, Flexible, Scalable, Modular, and Consistent Replacement for Limiting in
Circuit Simulation*, Scientific Computing in Electrical Engineering, 2020.
