"""
Test that `update_initializeprob!` correctly syncs `Initial()` values from
the current integrator state before re-solving algebraic constraints.

Uses a minimal switched DAE (SwitchedPump): one ODE + one algebraic + one
input control parameter σ ∈ {0,1}.

    ẋ = -α(x - σ·(H₀ - k·z²))     ODE: head decays toward pump curve
    0 = z - σ·√(max(0, H₀ - k·z²)) algebraic: flow from operating point

When σ changes mid-simulation, `reeval_internals_due_to_modification!`
triggers `initialize_dae!` → `OverrideInit` → `update_initializeprob!`.
The fix ensures:
  - `Initial(x)` is set from current `u` (not stale guess)
  - `Initial(σ)` is set from current `ps[σ]` (not stale construction value)
so the nonlinear init solve produces the correct algebraic state.
"""

using Test
using ModelingToolkit
using OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as Dt

# =====================================================================
# Minimal switched DAE
# =====================================================================

@component function SwitchedPump(; name)
    pars = @parameters begin
        H_0 = 1.0           # shut-off head
        k = 0.5             # resistance coefficient
        alpha = 5.0          # decay rate
    end
    vars = @variables begin
        x(t) = 0.5,   [guess = 0.5]   # differential: head
        z(t),          [guess = 0.0]   # algebraic: flow
        sigma(t),      [input = true]  # control: 0 or 1
    end
    eqs = [
        Dt(x) ~ -alpha * (x - sigma * (H_0 - k * z^2))
        z ~ sigma * sqrt(max(0, H_0 - k * z^2))
    ]
    return System(eqs, t, vars, pars; name)
end

# Analytic solution: z = σ√(H₀/(1+k))
z_analytic(H_0, k, sigma) = sigma * sqrt(H_0 / (1 + k))

# =====================================================================
# Test: manual stepping with reeval_internals_due_to_modification!
# =====================================================================

@testset "update_initializeprob! syncs Initial() values" begin
    @named sys = SwitchedPump()
    nns = toggle_namespacing(sys, false)
    csys = mtkcompile(sys; inputs=[nns.sigma], allow_parameter=false)

    H_0 = 1.0; k = 0.5
    z_on  = z_analytic(H_0, k, 1)
    z_off = z_analytic(H_0, k, 0)

    prob = ODEProblem(csys, [nns.sigma => 1], (0.0, Inf))
    integ = init(prob, Rodas5P(); dtmax=0.1)

    @test integ.isdae == true  # confirms DAE codepath is active

    # --- Steady state (pump ON) ---
    for _ in 1:200; step!(integ, 0.01, true); end
    @test integ.sol.retcode == ReturnCode.Success
    @test isapprox(integ.sol[nns.z][end], z_on; rtol=1e-3)

    # === Switch ON → OFF ===
    # Use integ[var] for current state (integ.sol[var][end] is the history, not yet updated)
    x_pre = integ[nns.x]

    integ.ps[nns.sigma] = 0
    SciMLBase.reeval_internals_due_to_modification!(integ, false)
    set_proposed_dt!(integ, 0.01)

    x_post = integ[nns.x]
    z_post = integ[nns.z]

    # Differential state x MUST be continuous across switch
    @test isapprox(x_post, x_pre; atol=1e-10)

    # Algebraic state z MUST satisfy new constraint (z=0 when σ=0)
    @test isapprox(z_post, z_off; atol=1e-6)

    # === Step forward in OFF mode ===
    for _ in 1:200; step!(integ, 0.01, true); end

    # === Switch OFF → ON ===
    x_pre2 = integ[nns.x]

    integ.ps[nns.sigma] = 1
    SciMLBase.reeval_internals_due_to_modification!(integ, false)
    set_proposed_dt!(integ, 0.01)

    x_post2 = integ[nns.x]
    z_post2 = integ[nns.z]

    # Differential state x MUST be continuous across switch
    @test isapprox(x_post2, x_pre2; atol=1e-10)

    # Algebraic state z MUST jump to ON operating point
    @test isapprox(z_post2, z_on; atol=1e-2)
end
