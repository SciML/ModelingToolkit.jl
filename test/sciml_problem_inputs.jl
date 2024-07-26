### Prepares Tests ###

# Fetch packages
using ModelingToolkit, JumpProcesses, NonlinearSolve, OrdinaryDiffEq, StaticArrays,
      SteadyStateDiffEq, StochasticDiffEq, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)
seed = rand(rng, 1:100)

### Basic Tests ###

# Prepares a models and initial conditions/parameters (of different forms) to be used as problem inputs.
begin
    # Prepare system components.
    @parameters kp kd k1 k2=0.5 Z0
    @variables X(t) Y(t) Z(t)=Z0
    alg_eqs = [
        0 ~ kp - k1 * X + k2 * Y - kd * X,
        0 ~ -k1 * Y + k1 * X - k2 * Y + k2 * Z,
        0 ~ k1 * Y - k2 * Z
    ]
    diff_eqs = [
        D(X) ~ kp - k1 * X + k2 * Y - kd * X,
        D(Y) ~ -k1 * Y + k1 * X - k2 * Y + k2 * Z,
        D(Z) ~ k1 * Y - k2 * Z
    ]
    noise_eqs = fill(0.01, 3, 6)
    jumps = [
        MassActionJump(kp, Pair{Symbolics.BasicSymbolic{Real}, Int64}[], [X => 1]),
        MassActionJump(kd, [X => 1], [X => -1]),
        MassActionJump(k1, [X => 1], [X => -1, Y => 1]),
        MassActionJump(k2, [Y => 1], [X => 1, Y => -1]),
        MassActionJump(k1, [Y => 1], [Y => -1, Z => 1]),
        MassActionJump(k2, [Z => 1], [Y => 1, Z => -1])
    ]

    # Create systems (without structural_simplify, since that might modify systems to affect intended tests).
    osys = complete(ODESystem(diff_eqs, t; name = :osys))
    ssys = complete(SDESystem(
        diff_eqs, noise_eqs, t, [X, Y, Z], [kp, kd, k1, k2]; name = :ssys))
    jsys = complete(JumpSystem(jumps, t, [X, Y, Z], [kp, kd, k1, k2]; name = :jsys))
    nsys = complete(NonlinearSystem(alg_eqs; name = :nsys))

    u0_alts = [
        # Vectors not providing default values.
        [X => 4, Y => 5],
        [osys.X => 4, osys.Y => 5],
        # Vectors providing default values.
        [X => 4, Y => 5, Z => 10],
        [osys.X => 4, osys.Y => 5, osys.Z => 10],
        # Static vectors not providing default values.
        SA[X => 4, Y => 5],
        SA[osys.X => 4, osys.Y => 5],
        # Static vectors providing default values.
        SA[X => 4, Y => 5, Z => 10],
        SA[osys.X => 4, osys.Y => 5, osys.Z => 10],
        # Dicts not providing default values.
        Dict([X => 4, Y => 5]),
        Dict([osys.X => 4, osys.Y => 5]),
        # Dicts providing default values.
        Dict([X => 4, Y => 5, Z => 10]),
        Dict([osys.X => 4, osys.Y => 5, osys.Z => 10]),
        # Tuples not providing default values.
        (X => 4, Y => 5),
        (osys.X => 4, osys.Y => 5),
        # Tuples providing default values.
        (X => 4, Y => 5, Z => 10),
        (osys.X => 4, osys.Y => 5, osys.Z => 10)
    ]
    tspan = (0.0, 10.0)
    p_alts = [
        # Vectors not providing default values.
        [kp => 1.0, kd => 0.1, k1 => 0.25, Z0 => 10],
        [osys.kp => 1.0, osys.kd => 0.1, osys.k1 => 0.25, osys.Z0 => 10],
        # Vectors providing default values.
        [kp => 1.0, kd => 0.1, k1 => 0.25, k2 => 0.5, Z0 => 10],
        [osys.kp => 1.0, osys.kd => 0.1, osys.k1 => 0.25, osys.k2 => 0.5, osys.Z0 => 10],
        # Static vectors not providing default values.
        SA[kp => 1.0, kd => 0.1, k1 => 0.25, Z0 => 10],
        SA[osys.kp => 1.0, osys.kd => 0.1, osys.k1 => 0.25, osys.Z0 => 10],
        # Static vectors providing default values.
        SA[kp => 1.0, kd => 0.1, k1 => 0.25, k2 => 0.5, Z0 => 10],
        SA[osys.kp => 1.0, osys.kd => 0.1, osys.k1 => 0.25, osys.k2 => 0.5, osys.Z0 => 10],
        # Dicts not providing default values.
        Dict([kp => 1.0, kd => 0.1, k1 => 0.25, Z0 => 10]),
        Dict([osys.kp => 1.0, osys.kd => 0.1, osys.k1 => 0.25, osys.Z0 => 10]),
        # Dicts providing default values.
        Dict([kp => 1.0, kd => 0.1, k1 => 0.25, k2 => 0.5, Z0 => 10]),
        Dict([osys.kp => 1.0, osys.kd => 0.1, osys.k1 => 0.25,
            osys.k2 => 0.5, osys.Z0 => 10]),
        # Tuples not providing default values.
        (kp => 1.0, kd => 0.1, k1 => 0.25, Z0 => 10),
        (osys.kp => 1.0, osys.kd => 0.1, osys.k1 => 0.25, osys.Z0 => 10),
        # Tuples providing default values.
        (kp => 1.0, kd => 0.1, k1 => 0.25, k2 => 0.5, Z0 => 10),
        (osys.kp => 1.0, osys.kd => 0.1, osys.k1 => 0.25, osys.k2 => 0.5, osys.Z0 => 10)
    ]
end

# Perform ODE simulations (singular and ensemble).
let
    # Creates normal and ensemble problems.
    base_oprob = ODEProblem(osys, u0_alts[1], tspan, p_alts[1])
    base_sol = solve(base_oprob, Tsit5(); saveat = 1.0)
    base_eprob = EnsembleProblem(base_oprob)
    base_esol = solve(base_eprob, Tsit5(); trajectories = 2, saveat = 1.0)

    # Simulates problems for all input types, checking that identical solutions are found.
    # test failure.
    for u0 in u0_alts, p in p_alts
        oprob = remake(base_oprob; u0, p)
        @test base_sol == solve(oprob, Tsit5(); saveat = 1.0)
        eprob = remake(base_eprob; u0, p)
        @test base_esol == solve(eprob, Tsit5(); trajectories = 2, saveat = 1.0)
    end
end

# Solves a nonlinear problem (EnsembleProblems are not possible for these).
let
    base_nlprob = NonlinearProblem(nsys, u0_alts[1], p_alts[1])
    base_sol = solve(base_nlprob, NewtonRaphson())
    # Solves problems for all input types, checking that identical solutions are found.
    for u0 in u0_alts, p in p_alts
        nlprob = remake(base_nlprob; u0, p)
        @test base_sol == solve(nlprob, NewtonRaphson())
    end
end

# Perform steady state simulations (singular and ensemble).
let
    # Creates normal and ensemble problems.
    base_ssprob = SteadyStateProblem(osys, u0_alts[1], p_alts[1])
    base_sol = solve(base_ssprob, DynamicSS(Tsit5()))
    base_eprob = EnsembleProblem(base_ssprob)
    base_esol = solve(base_eprob, DynamicSS(Tsit5()); trajectories = 2)

    # Simulates problems for all input types, checking that identical solutions are found.
    # test failure.
    for u0 in u0_alts, p in p_alts
        ssprob = remake(base_ssprob; u0, p)
        @test base_sol == solve(ssprob, DynamicSS(Tsit5()))
        eprob = remake(base_eprob; u0, p)
        @test base_esol == solve(eprob, DynamicSS(Tsit5()); trajectories = 2)
    end
end
