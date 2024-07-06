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

# Perform SDE simulations (singular and ensemble).
let
    # Creates normal and ensemble problems.
    base_sprob = SDEProblem(ssys, u0_alts[1], tspan, p_alts[1])
    base_sol = solve(base_sprob, ImplicitEM(); seed, saveat = 1.0)
    base_eprob = EnsembleProblem(base_sprob)
    base_esol = solve(base_eprob, ImplicitEM(); seed, trajectories = 2, saveat = 1.0)

    # Simulates problems for all input types, checking that identical solutions are found.
    @test_broken false # first remake in subsequent test yields a `ERROR: type Nothing has no field portion`.
    for u0 in u0_alts, p in p_alts
        #    sprob = remake(base_sprob; u0, p)
        #    @test base_sol == solve(sprob, ImplicitEM(); seed, saveat = 1.0)
        #    eprob = remake(base_eprob; u0, p)
        #    @test base_esol == solve(eprob, ImplicitEM(); seed, trajectories = 2, saveat = 1.0)
    end
end

# Perform jump simulations (singular and ensemble).
let
    # Creates normal and ensemble problems.
    base_dprob = DiscreteProblem(jsys, u0_alts[1], tspan, p_alts[1])
    base_jprob = JumpProblem(jsys, base_dprob, Direct(); rng)
    base_sol = solve(base_jprob, SSAStepper(); seed, saveat = 1.0)
    base_eprob = EnsembleProblem(base_jprob)
    base_esol = solve(base_eprob, SSAStepper(); seed, trajectories = 2, saveat = 1.0)

    # Simulates problems for all input types, checking that identical solutions are found.
    @test_broken false # first remake in subsequent test yields a `ERROR: type Nothing has no field portion`.
    for u0 in u0_alts, p in p_alts
        #    jprob = remake(base_jprob; u0, p)
        #    @test base_sol == solve(base_jprob, SSAStepper(); seed, saveat = 1.0)
        #    eprob = remake(base_eprob; u0, p)
        #    @test base_esol == solve(eprob, SSAStepper(); seed, trajectories = 2, saveat = 1.0)
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

### Checks Errors On Faulty Inputs ###

# Checks various erroneous problem inputs, ensuring that these throw errors.
let
    # Prepare system components.
    @parameters k1 k2 k3
    @variables X1(t) X2(t) X3(t)
    alg_eqs = [
        0 ~ -k1 * X1 + k2 * X2,
        0 ~ k1 * X1 - k2 * X2
    ]
    diff_eqs = [
        D(X1) ~ -k1 * X1 + k2 * X2,
        D(X2) ~ k1 * X1 - k2 * X2
    ]
    noise_eqs = fill(0.01, 2, 2)
    jumps = [
        MassActionJump(k1, [X1 => 1], [X1 => -1, X2 => 1]),
        MassActionJump(k2, [X2 => 1], [X1 => 1, X2 => -1])
    ]

    # Create systems (without structural_simplify, since that might modify systems to affect intended tests).
    osys = complete(ODESystem(diff_eqs, t; name = :osys))
    ssys = complete(SDESystem(diff_eqs, noise_eqs, t, [X1, X2], [k1, k2]; name = :ssys))
    jsys = complete(JumpSystem(jumps, t, [X1, X2], [k1, k2]; name = :jsys))
    nsys = complete(NonlinearSystem(alg_eqs; name = :nsys))

    # Declares valid initial conditions and parameter values
    u0_valid = [X1 => 1, X2 => 2]
    ps_valid = [k1 => 0.5, k2 => 0.1]

    # Declares invalid initial conditions and parameters. This includes both cases where values are
    # missing, or additional ones are given. Includes vector/Tuple/Dict forms.
    u0s_invalid = [
        # Missing a value.
        [X1 => 1],
        [osys.X1 => 1],
        SA[X1 => 1],
        SA[osys.X1 => 1],
        Dict([X1 => 1]),
        Dict([osys.X1 => 1]),
        (X1 => 1),
        (osys.X1 => 1),
        # Contain an additional value.
        [X1 => 1, X2 => 2, X3 => 3],
        SA[X1 => 1, X2 => 2, X3 => 3],
        Dict([X1 => 1, X2 => 2, X3 => 3]),
        (X1 => 1, X2 => 2, X3 => 3)
    ]
    ps_invalid = [
        # Missing a value.
        [k1 => 1.0],
        [osys.k1 => 1.0],
        SA[k1 => 1.0],
        SA[osys.k1 => 1.0],
        Dict([k1 => 1.0]),
        Dict([osys.k1 => 1.0]),
        (k1 => 1.0),
        (osys.k1 => 1.0),
        # Contain an additional value.
        [k1 => 1.0, k2 => 2.0, k3 => 3.0],
        SA[k1 => 1.0, k2 => 2.0, k3 => 3.0],
        Dict([k1 => 1.0, k2 => 2.0, k3 => 3.0]),
        (k1 => 1.0, k2 => 2.0, k3 => 3.0)
    ]

    # Loops through all potential parameter sets, checking their inputs yield errors.
    # Broken tests are due to this issue: https://github.com/SciML/ModelingToolkit.jl/issues/2779
    for ps in [[ps_valid]; ps_invalid], u0 in [[u0_valid]; u0s_invalid]
        # Handles problems with/without tspan separately. Special check ensuring that valid inputs passes.
        for (xsys, XProblem) in zip([osys, ssys, jsys],
            [ODEProblem, SDEProblem, DiscreteProblem])
            if isequal(ps, ps_valid) && isequal(u0, u0_valid)
                XProblem(xsys, u0, (0.0, 1.0), ps)
                @test true
            else
                @test_broken false
                continue
                @test_throws Exception XProblem(xsys, u0, (0.0, 1.0), ps)
            end
        end
        for (xsys, XProblem) in zip([nsys, osys], [NonlinearProblem, SteadyStateProblem])
            if isequal(ps, ps_valid) && isequal(u0, u0_valid)
                XProblem(xsys, u0, ps)
                @test true
            else
                @test_broken false
                continue
                @test_throws Exception XProblem(xsys, u0, ps)
            end
        end
    end
end

### Vector Parameter/Variable Inputs ###

begin
    # Declares the model (with vector species/parameters, with/without default values, and observables).
    @variables X(t)[1:2] Y(t)[1:2]=[10.0, 20.0] XY(t)[1:2]
    @parameters p[1:2] d[1:2]=[0.2, 0.5]
    alg_eqs = [
        0 ~ p[1] - d[1] * X[1],
        0 ~ p[2] - d[2] * X[2],
        0 ~ p[1] - d[1] * Y[1],
        0 ~ p[2] - d[2] * Y[2]
    ]
    diff_eqs = [
        D(X[1]) ~ p[1] - d[1] * X[1],
        D(X[2]) ~ p[2] - d[2] * X[2],
        D(Y[1]) ~ p[1] - d[1] * Y[1],
        D(Y[2]) ~ p[2] - d[2] * Y[2]
    ]
    noise_eqs = fill(0.01, 4, 8)
    jumps = [
        MassActionJump(p[1], Pair{Symbolics.BasicSymbolic{Real}, Int64}[], [X[1] => 1]),
        MassActionJump(p[2], Pair{Symbolics.BasicSymbolic{Real}, Int64}[], [X[2] => 1]),
        MassActionJump(d[1], [X[1] => 1], [X[1] => -1]),
        MassActionJump(d[2], [X[2] => 1], [X[2] => -1]),
        MassActionJump(p[1], Pair{Symbolics.BasicSymbolic{Real}, Int64}[], [Y[1] => 1]),
        MassActionJump(p[2], Pair{Symbolics.BasicSymbolic{Real}, Int64}[], [Y[2] => 1]),
        MassActionJump(d[1], [Y[1] => 1], [Y[1] => -1]),
        MassActionJump(d[2], [Y[2] => 1], [Y[2] => -1])
    ]
    observed = [XY[1] ~ X[1] + Y[1], XY[2] ~ X[2] + Y[2]]

    # Create systems (without structural_simplify, since that might modify systems to affect intended tests).
    osys = complete(ODESystem(diff_eqs, t; observed, name = :osys))
    ssys = complete(SDESystem(
        diff_eqs, noise_eqs, t, [X[1], X[2], Y[1], Y[2]], [p, d]; observed, name = :ssys))
    jsys = complete(JumpSystem(
        jumps, t, [X[1], X[2], Y[1], Y[2]], [p, d]; observed, name = :jsys))
    nsys = complete(NonlinearSystem(alg_eqs; observed, name = :nsys))

    # Declares various u0 versions (scalarised and vector forms).
    u0_alts_vec = [
        # Vectors not providing default values.
        [X => [1.0, 2.0]],
        [X[1] => 1.0, X[2] => 2.0],
        [osys.X => [1.0, 2.0]],
        [osys.X[1] => 1.0, osys.X[2] => 2.0],
        # Vectors providing default values.
        [X => [1.0, 2.0], Y => [10.0, 20.0]],
        [X[1] => 1.0, X[2] => 2.0, Y[1] => 10.0, Y[2] => 20.0],
        [osys.X => [1.0, 2.0], osys.Y => [10.0, 20.0]],
        [osys.X[1] => 1.0, osys.X[2] => 2.0, osys.Y[1] => 10.0, osys.Y[2] => 20.0],
        # Static vectors not providing default values.
        SA[X => [1.0, 2.0]],
        SA[X[1] => 1.0, X[2] => 2.0],
        SA[osys.X => [1.0, 2.0]],
        SA[osys.X[1] => 1.0, osys.X[2] => 2.0],
        # Static vectors providing default values.
        SA[X => [1.0, 2.0], Y => [10.0, 20.0]],
        SA[X[1] => 1.0, X[2] => 2.0, Y[1] => 10.0, Y[2] => 20.0],
        SA[osys.X => [1.0, 2.0], osys.Y => [10.0, 20.0]],
        SA[osys.X[1] => 1.0, osys.X[2] => 2.0, osys.Y[1] => 10.0, osys.Y[2] => 20.0],
        # Dicts not providing default values.
        Dict([X => [1.0, 2.0]]),
        Dict([X[1] => 1.0, X[2] => 2.0]),
        Dict([osys.X => [1.0, 2.0]]),
        Dict([osys.X[1] => 1.0, osys.X[2] => 2.0]),
        # Dicts providing default values.
        Dict([X => [1.0, 2.0], Y => [10.0, 20.0]]),
        Dict([X[1] => 1.0, X[2] => 2.0, Y[1] => 10.0, Y[2] => 20.0]),
        Dict([osys.X => [1.0, 2.0], osys.Y => [10.0, 20.0]]),
        Dict([osys.X[1] => 1.0, osys.X[2] => 2.0, osys.Y[1] => 10.0, osys.Y[2] => 20.0]),
        # Tuples not providing default values.
        (X => [1.0, 2.0]),
        (X[1] => 1.0, X[2] => 2.0),
        (osys.X => [1.0, 2.0]),
        (osys.X[1] => 1.0, osys.X[2] => 2.0),
        # Tuples providing default values.
        (X => [1.0, 2.0], Y => [10.0, 20.0]),
        (X[1] => 1.0, X[2] => 2.0, Y[1] => 10.0, Y[2] => 20.0),
        (osys.X => [1.0, 2.0], osys.Y => [10.0, 20.0]),
        (osys.X[1] => 1.0, osys.X[2] => 2.0, osys.Y[1] => 10.0, osys.Y[2] => 20.0)
    ]

    # Declares various ps versions (vector forms only).
    p_alts_vec = [
        # Vectors not providing default values.
        [p => [1.0, 2.0]],
        [osys.p => [1.0, 2.0]],
        # Vectors providing default values.
        [p => [4.0, 5.0], d => [0.2, 0.5]],
        [osys.p => [4.0, 5.0], osys.d => [0.2, 0.5]],
        # Static vectors not providing default values.
        SA[p => [1.0, 2.0]],
        SA[osys.p => [1.0, 2.0]],
        # Static vectors providing default values.
        SA[p => [4.0, 5.0], d => [0.2, 0.5]],
        SA[osys.p => [4.0, 5.0], osys.d => [0.2, 0.5]],
        # Dicts not providing default values.
        Dict([p => [1.0, 2.0]]),
        Dict([osys.p => [1.0, 2.0]]),
        # Dicts providing default values.
        Dict([p => [4.0, 5.0], d => [0.2, 0.5]]),
        Dict([osys.p => [4.0, 5.0], osys.d => [0.2, 0.5]]),
        # Tuples not providing default values.
        (p => [1.0, 2.0]),
        (osys.p => [1.0, 2.0]),
        # Tuples providing default values.
        (p => [4.0, 5.0], d => [0.2, 0.5]),
        (osys.p => [4.0, 5.0], osys.d => [0.2, 0.5])
    ]

    # Declares a timespan.
    tspan = (0.0, 10.0)
end

# Perform ODE simulations (singular and ensemble).
let
    # Creates normal and ensemble problems.
    base_oprob = ODEProblem(osys, u0_alts_vec[1], tspan, p_alts_vec[1])
    base_sol = solve(base_oprob, Tsit5(); saveat = 1.0)
    base_eprob = EnsembleProblem(base_oprob)
    base_esol = solve(base_eprob, Tsit5(); trajectories = 2, saveat = 1.0)

    # Simulates problems for all input types, checking that identical solutions are found.
    @test_broken false # Does not work for certain inputs, likely related to https://github.com/SciML/ModelingToolkit.jl/issues/2804.
    for u0 in u0_alts_vec, p in p_alts_vec
        oprob = remake(base_oprob; u0, p)
        # @test base_sol == solve(oprob, Tsit5(); saveat = 1.0)
        eprob = remake(base_eprob; u0, p)
        # @test base_esol == solve(eprob, Tsit5(); trajectories = 2, saveat = 1.0)
    end
end

# Perform SDE simulations (singular and ensemble).
let
    # Creates normal and ensemble problems.
    base_sprob = SDEProblem(ssys, u0_alts_vec[1], tspan, p_alts_vec[1])
    base_sol = solve(base_sprob, ImplicitEM(); seed, saveat = 1.0)
    base_eprob = EnsembleProblem(base_sprob)
    base_esol = solve(base_eprob, ImplicitEM(); seed, trajectories = 2, saveat = 1.0)

    # Simulates problems for all input types, checking that identical solutions are found.
    @test_broken false # Does not work for certain inputs, likely related to https://github.com/SciML/ModelingToolkit.jl/issues/2804.
    for u0 in u0_alts_vec, p in p_alts_vec
        sprob = remake(base_sprob; u0, p)
        # @test base_sol == solve(sprob, ImplicitEM(); seed, saveat = 1.0)
        eprob = remake(base_eprob; u0, p)
        # @test base_esol == solve(eprob, ImplicitEM(); seed, trajectories = 2, saveat = 1.0)
    end
end

# Perform jump simulations (singular and ensemble).
# Fails. At least partially related to https://github.com/SciML/ModelingToolkit.jl/issues/2804.
@test_broken let
    # Creates normal and ensemble problems.
    base_dprob = DiscreteProblem(jsys, u0_alts_vec[1], tspan, p_alts_vec[1])
    base_jprob = JumpProblem(jsys, base_dprob, Direct(); rng)
    base_sol = solve(base_jprob, SSAStepper(); seed, saveat = 1.0)
    base_eprob = EnsembleProblem(base_jprob)
    base_esol = solve(base_eprob, SSAStepper(); seed, trajectories = 2, saveat = 1.0)

    # Simulates problems for all input types, checking that identical solutions are found.
    @test_broken false # Does not work for certain inputs, likely related to https://github.com/SciML/ModelingToolkit.jl/issues/2804.
    for u0 in u0_alts_vec, p in p_alts_vec
        jprob = remake(base_jprob; u0, p)
        # @test base_sol == solve(base_jprob, SSAStepper(); seed, saveat = 1.0)
        eprob = remake(base_eprob; u0, p)
        # @test base_esol == solve(eprob, SSAStepper(); seed, trajectories = 2, saveat = 1.0)
    end
end

# Solves a nonlinear problem (EnsembleProblems are not possible for these).
let
    base_nlprob = NonlinearProblem(nsys, u0_alts_vec[1], p_alts_vec[1])
    base_sol = solve(base_nlprob, NewtonRaphson())
    @test_broken false # Does not work for certain inputs, likely related to https://github.com/SciML/ModelingToolkit.jl/issues/2804.
    for u0 in u0_alts_vec, p in p_alts_vec
        nlprob = remake(base_nlprob; u0, p)
        # @test base_sol == solve(nlprob, NewtonRaphson())
    end
end

# Perform steady state simulations (singular and ensemble).
# Fails. At least partially related to https://github.com/SciML/ModelingToolkit.jl/issues/2804.
let
    # Creates normal and ensemble problems.
    base_ssprob = SteadyStateProblem(osys, u0_alts_vec[1], p_alts_vec[1])
    base_sol = solve(base_ssprob, DynamicSS(Tsit5()))
    base_eprob = EnsembleProblem(base_ssprob)
    base_esol = solve(base_eprob, DynamicSS(Tsit5()); trajectories = 2)

    # Simulates problems for all input types, checking that identical solutions are found.
    @test_broken false # Does not work for certain inputs, likely related to https://github.com/SciML/ModelingToolkit.jl/issues/2804.
    for u0 in u0_alts_vec, p in p_alts_vec
        ssprob = remake(base_ssprob; u0, p)
        # @test base_sol == solve(ssprob, DynamicSS(Tsit5()))
        eprob = remake(base_eprob; u0, p)
        # @test base_esol == solve(eprob, DynamicSS(Tsit5()); trajectories = 2)
    end
end
