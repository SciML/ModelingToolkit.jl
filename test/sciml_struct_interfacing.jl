### Prepares Tests ###

# Fetch packages
using ModelingToolkit, JumpProcesses, NonlinearSolve, OrdinaryDiffEq, Plots, SteadyStateDiffEq, StochasticDiffEq, Test
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit: getp, getu, setp, setu

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)
seed = rand(rng, 1:100)


### Basic Tests ###

# Prepares a model systems.
begin
    # Prepare system components.
    @parameters kp kd k1 k2
    @variables X(t) Y(t) XY(t)
    alg_eqs = [
        0 ~ kp - kd*X  - k1*X + k2*Y
        0 ~ 1 + k1*X - k2*Y - Y
    ]
    diff_eqs = [
        D(X) ~ kp - kd*X  - k1*X + k2*Y
        D(Y) ~ 1 + k1*X - k2*Y - Y
    ]
    noise_eqs = [
        sqrt(kp + X), 
        sqrt(k1 + Y)
    ]
    jumps = [
        ConstantRateJump(kp, [X ~ X + 1]),
        ConstantRateJump(kd*X, [X ~ X - 1]),
        ConstantRateJump(k1*X, [X ~ X - 1, Y ~ Y + 1]),
        ConstantRateJump(k2*Y, [X ~ X + 1, Y ~ Y - 1]),
        ConstantRateJump(1, [Y ~ Y + 1]),
        ConstantRateJump(Y, [Y ~ Y - 1]),
    ]
    observed = [XY ~ X + Y]

    # Create systems (without structural_simplify, since that might modify systems to affect intended tests).
    osys = complete(ODESystem(diff_eqs, t; observed, name = :osys))
    ssys = complete(SDESystem(diff_eqs, noise_eqs, t, [X, Y], [kp, kd, k1, k2]; observed, name = :ssys))
    jsys = complete(JumpSystem(jumps, t, [X, Y], [kp, kd, k1, k2]; observed, name = :jsys))
    nsys = complete(NonlinearSystem(alg_eqs; observed, name = :nsys))
end


# Prepares problems, integrators, and solutions.
begin
    # Sets problem inputs (to be used for all problem creations).
    u0_vals = [X => 4, Y => 5]
    tspan = (0.0, 10.0)
    p_vals = [kp => 1.0, kd => 0.1, k1 => 0.25, k2 => 0.5]

    # Creates problems.
    oprob = ODEProblem(osys, u0_vals, tspan, p_vals)
    sprob = SDEProblem(ssys,u0_vals, tspan, p_vals)
    dprob = DiscreteProblem(jsys, u0_vals, tspan, p_vals)
    jprob = JumpProblem(jsys, deepcopy(dprob), Direct(); rng)
    nprob = NonlinearProblem(nsys, u0_vals, p_vals)
    ssprob = SteadyStateProblem(osys, u0_vals, p_vals)
    problems = [oprob, sprob, dprob, jprob, nprob, ssprob]
    systems = [osys, ssys, jsys, jsys, nsys, osys]

    # Creates an `EnsembleProblem` for each problem.
    eoprob = EnsembleProblem(oprob)
    esprob = EnsembleProblem(sprob)
    edprob = EnsembleProblem(dprob)
    ejprob = EnsembleProblem(jprob)
    enprob = EnsembleProblem(nprob)
    essprob = EnsembleProblem(ssprob)
    eproblems = [eoprob, esprob, edprob, ejprob, enprob, essprob]
    esystems = [osys, ssys, jsys, jsys, nsys, osys]

    # Creates integrators.
    oint = init(oprob, Tsit5(); save_everystep = false)
    sint = init(sprob, ImplicitEM(); save_everystep = false)
    jint = init(jprob, SSAStepper())
    nint = init(nprob, NewtonRaphson(); save_everystep = false)
    @test_broken ssint = init(ssprob, DynamicSS(Tsit5()); save_everystep = false) # https://github.com/SciML/SteadyStateDiffEq.jl/issues/79
    integrators = [oint, sint, jint, nint]
    
    # Creates solutions.
    osol = solve(oprob, Tsit5())
    ssol = solve(sprob, ImplicitEM(); seed)
    jsol = solve(jprob, SSAStepper(); seed)
    nsol = solve(nprob, NewtonRaphson())
    sssol = solve(ssprob, DynamicSS(Tsit5()))
    sols = [osol, ssol, jsol, nsol, sssol]
end

# Tests problem indexing and updating.
let 
    @test_broken false # Currently does not work for nonlinearproblems and their ensemble problems (https://github.com/SciML/SciMLBase.jl/issues/720).
    # for (prob, sys) in zip([deepcopy(problems); deepcopy(eproblems)], [deepcopy(systems); deepcopy(esystems)])
    for (prob, sys) in zip([deepcopy([oprob, sprob, dprob, jprob, ssprob]); deepcopy([eoprob, esprob, edprob, ejprob, essprob])], [deepcopy([osys, ssys, jsys, jsys, osys]); deepcopy([osys, ssys, jsys, jsys, osys])])
        # Get u values (including observables).
        @test prob[X] == prob[sys.X] == prob[:X] == 4
        @test prob[XY] == prob[sys.XY] == prob[:XY] == 9
        @test prob[[XY,Y]] == prob[[sys.XY,sys.Y]] == prob[[:XY,:Y]] == [9, 5]
        @test_broken prob[(XY,Y)] == prob[(sys.XY,sys.Y)] == prob[(:XY,:Y)] == (9, 5) # https://github.com/SciML/SciMLBase.jl/issues/709
        @test getu(prob, X)(prob) == getu(prob, sys.X)(prob) == getu(prob, :X)(prob) == 4
        @test getu(prob, XY)(prob) == getu(prob, sys.XY)(prob) == getu(prob, :XY)(prob) == 9 
        @test getu(prob, [XY,Y])(prob) == getu(prob, [sys.XY,sys.Y])(prob) == getu(prob, [:XY,:Y])(prob) == [9, 5]  
        @test getu(prob, (XY,Y))(prob) == getu(prob, (sys.XY,sys.Y))(prob) == getu(prob, (:XY,:Y))(prob) == (9, 5)

        # Set u values.
        prob[X] = 20
        @test prob[X] == 20
        prob[sys.X] = 30
        @test prob[X] == 30
        prob[:X] = 40
        @test prob[X] == 40
        setu(prob, X)(prob, 50)
        @test prob[X] == 50
        setu(prob, sys.X)(prob, 60)
        @test prob[X] == 60
        setu(prob, :X)(prob, 70)
        @test prob[X] == 70

        # Get p values.
        @test prob.ps[kp] == prob.ps[sys.kp] == prob.ps[:kp] == 1.0    
        @test prob.ps[[k1,k2]] == prob.ps[[sys.k1,sys.k2]] == prob.ps[[:k1,:k2]] == [0.25, 0.5]
        @test prob.ps[(k1,k2)] == prob.ps[(sys.k1,sys.k2)] == prob.ps[(:k1,:k2)] == (0.25, 0.5)
        @test getp(prob, kp)(prob) == getp(prob, sys.kp)(prob) == getp(prob, :kp)(prob) == 1.0
        @test getp(prob, [k1,k2])(prob) == getp(prob, [sys.k1,sys.k2])(prob) == getp(prob, [:k1,:k2])(prob) == [0.25, 0.5]
        @test getp(prob, (k1,k2))(prob) == getp(prob, (sys.k1,sys.k2))(prob) == getp(prob, (:k1,:k2))(prob) == (0.25, 0.5)
        
        # Set p values.
        prob.ps[kp] = 2.0
        @test prob.ps[kp] == 2.0
        prob.ps[sys.kp] = 3.0
        @test prob.ps[kp] == 3.0
        prob.ps[:kp] = 4.0
        @test prob.ps[kp] == 4.0
        setp(prob, kp)(prob, 5.0)
        @test prob.ps[kp] == 5.0
        setp(prob, sys.kp)(prob, 6.0)
        @test prob.ps[kp] == 6.0
        setp(prob, :kp)(prob, 7.0)
        @test prob.ps[kp] == 7.0
    end
end

# Test remake function.
let 
    for (prob, sys) in zip([deepcopy(problems); deepcopy(eproblems)], [deepcopy(systems); deepcopy(esystems)])
        # Remake for all u0s.
        rp = remake(prob; u0 = [X => 1, Y => 2])
        @test rp[[X, Y]] == [1, 2]
        rp = remake(prob; u0 = [sys.X => 3, sys.Y => 4])
        @test rp[[X, Y]] == [3, 4]
        rp = remake(prob; u0 = [:X => 5, :Y => 6])
        @test rp[[X, Y]] == [5, 6]

        # Remake for a single u0.
        rp = remake(prob; u0 = [Y => 7])
        @test rp[[X, Y]] == [4, 7]
        rp = remake(prob; u0 = [sys.Y => 8])
        @test rp[[X, Y]] == [4, 8]
        rp = remake(prob; u0 = [:Y => 9])
        @test rp[[X, Y]] == [4, 9]

        # Remake for all ps.
        rp = remake(prob; p = [kp => 1.0, kd => 2.0, k1 => 3.0, k2 => 4.0])
        @test rp.ps[[kp, kd, k1, k2]] == [1.0, 2.0, 3.0, 4.0]
        rp = remake(prob; p = [sys.kp => 5.0, sys.kd => 6.0, sys.k1 => 7.0, sys.k2 => 8.0])
        @test rp.ps[[kp, kd, k1, k2]] == [5.0, 6.0, 7.0, 8.0]
        rp = remake(prob; p = [:kp => 9.0, :kd => 10.0, :k1 => 11.0, :k2 => 12.0])
        @test rp.ps[[kp, kd, k1, k2]] == [9.0, 10.0, 11.0, 12.0]

        # Remake for a single p.
        rp = remake(prob; p = [k2 => 13.0])
        @test rp.ps[[kp, kd, k1, k2]] == [1.0, 0.1, 0.25, 13.0]
        rp = remake(prob; p = [sys.k2 => 14.0])
        @test rp.ps[[kp, kd, k1, k2]] == [1.0, 0.1, 0.25, 14.0]
        rp = remake(prob; p = [:k2 => 15.0])
        @test rp.ps[[kp, kd, k1, k2]] == [1.0, 0.1, 0.25, 15.0]
    end
end

# Test integrator indexing.
let 
    @test_broken false # NOTE: Multiple problems for `nint` (https://github.com/SciML/SciMLBase.jl/issues/662).
    for (int, sys) in zip(deepcopy([oint, sint, jint]), [osys, ssys, jsys])
        # Get u values.
        @test int[X] == int[sys.X] == int[:X] == 4
        @test int[XY] == int[sys.XY] == int[:XY] == 9
        @test int[[XY,Y]] == int[[sys.XY,sys.Y]] == int[[:XY,:Y]] == [9, 5]
        @test int[(XY,Y)] == int[(sys.XY,sys.Y)] == int[(:XY,:Y)] == (9, 5)
        @test getu(int, X)(int) == getu(int, sys.X)(int) == getu(int, :X)(int) == 4
        @test getu(int, XY)(int) == getu(int, sys.XY)(int) == getu(int, :XY)(int) == 9 
        @test getu(int, [XY,Y])(int) == getu(int, [sys.XY,sys.Y])(int) == getu(int, [:XY,:Y])(int) == [9, 5]  
        @test getu(int, (XY,Y))(int) == getu(int, (sys.XY,sys.Y))(int) == getu(int, (:XY,:Y))(int) == (9, 5)

        # Set u values.
        int[X] = 20
        @test int[X] == 20
        int[sys.X] = 30
        @test int[X] == 30
        int[:X] = 40
        @test int[X] == 40
        setu(int, X)(int, 50)
        @test int[X] == 50
        setu(int, sys.X)(int, 60)
        @test int[X] == 60
        setu(int, :X)(int, 70)
        @test int[X] == 70

        # Get p values.
        @test int.ps[kp] == int.ps[sys.kp] == int.ps[:kp] == 1.0    
        @test int.ps[[k1,k2]] == int.ps[[sys.k1,sys.k2]] == int.ps[[:k1,:k2]] == [0.25, 0.5]
        @test int.ps[(k1,k2)] == int.ps[(sys.k1,sys.k2)] == int.ps[(:k1,:k2)] == (0.25, 0.5)
        @test getp(int, kp)(int) == getp(int, sys.kp)(int) == getp(int, :kp)(int) == 1.0
        @test getp(int, [k1,k2])(int) == getp(int, [sys.k1,sys.k2])(int) == getp(int, [:k1,:k2])(int) == [0.25, 0.5]
        @test getp(int, (k1,k2))(int) == getp(int, (sys.k1,sys.k2))(int) == getp(int, (:k1,:k2))(int) == (0.25, 0.5)
        
        # Set p values.
        int.ps[kp] = 2.0
        @test int.ps[kp] == 2.0
        int.ps[sys.kp] = 3.0
        @test int.ps[kp] == 3.0
        int.ps[:kp] = 4.0
        @test int.ps[kp] == 4.0
        setp(int, kp)(int, 5.0)
        @test int.ps[kp] == 5.0
        setp(int, sys.kp)(int, 6.0)
        @test int.ps[kp] == 6.0
        setp(int, :kp)(int, 7.0)
        @test int.ps[kp] == 7.0
    end
end

# Test solve's save_idxs argument.
# Currently, `save_idxs` is broken with symbolic stuff (https://github.com/SciML/ModelingToolkit.jl/issues/1761).
let 
    for (prob, sys, solver) in zip(deepcopy([oprob, sprob, jprob]), [osys, ssys, jsys], [Tsit5(), ImplicitEM(), SSAStepper()])
        # Save single variable
        @test_broken solve(prob, solver; seed, save_idxs=X)[X][1] == 4
        @test_broken solve(prob, solver; seed, save_idxs=sys.X)[X][1] == 4
        @test_broken solve(prob, solver; seed, save_idxs=:X)[X][1] == 4

        # Save observable.
        @test_broken solve(prob, solver; seed, save_idxs=XY)[XY][1] == 9
        @test_broken solve(prob, solver; seed, save_idxs=sys.XY)[XY][1] == 9
        @test_broken solve(prob, solver; seed, save_idxs=:XY)[XY][1] == 9

        # Save vector of stuff.
        @test_broken solve(prob, solver; seed, save_idxs=[XY,Y])[[XY,Y]][1] == [9, 5]
        @test_broken solve(prob, solver; seed, save_idxs=[sys.XY,sys.Y])[[sys.XY,sys.Y]][1] == [9, 5]
        @test_broken solve(prob, solver; seed, save_idxs=[:XY,:Y])[[:XY,:Y]][1] == [9, 5]
    end
end

# Tests solution indexing.
let 
    for (sol, sys) in zip(deepcopy([osol, ssol, jsol]), [osys, ssys, jsys])
        # Get u values.
        @test sol[X][1] == sol[sys.X][1] == sol[:X][1] == 4
        @test sol[XY][1] == sol[sys.XY][1] == sol[:XY][1] == 9
        @test sol[[XY,Y]][1] == sol[[sys.XY,sys.Y]][1] == sol[[:XY,:Y]][1] == [9, 5]
        @test sol[(XY,Y)][1] == sol[(sys.XY,sys.Y)][1] == sol[(:XY,:Y)][1] == (9, 5)
        @test getu(sol, X)(sol)[1] == getu(sol, sys.X)(sol)[1] == getu(sol, :X)(sol)[1] == 4
        @test getu(sol, XY)(sol)[1] == getu(sol, sys.XY)(sol)[1] == getu(sol, :XY)(sol)[1] == 9 
        @test getu(sol, [XY,Y])(sol)[1] == getu(sol, [sys.XY,sys.Y])(sol)[1] == getu(sol, [:XY,:Y])(sol)[1] == [9, 5]  
        @test getu(sol, (XY,Y))(sol)[1] == getu(sol, (sys.XY,sys.Y))(sol)[1] == getu(sol, (:XY,:Y))(sol)[1] == (9, 5)       

        # Get u values via idxs and functional call.
        @test osol(0.0; idxs=X) == osol(0.0; idxs=sys.X) == osol(0.0; idxs=:X) == 4
        @test osol(0.0; idxs=XY) == osol(0.0; idxs=sys.XY) == osol(0.0; idxs=:XY) == 9
        @test osol(0.0; idxs = [XY,Y]) == osol(0.0; idxs = [sys.XY,sys.Y]) == osol(0.0; idxs = [:XY,:Y]) == [9, 5]
        @test_broken osol(0.0; idxs = (XY,Y)) == osol(0.0; idxs = (sys.XY,sys.Y)) == osol(0.0; idxs = (:XY,:Y)) == (9, 5) # https://github.com/SciML/SciMLBase.jl/issues/711

        # Get p values.
        @test sol.ps[kp] == sol.ps[sys.kp] == sol.ps[:kp] == 1.0    
        @test sol.ps[[k1,k2]] == sol.ps[[sys.k1,sys.k2]] == sol.ps[[:k1,:k2]] == [0.25, 0.5]
        @test sol.ps[(k1,k2)] == sol.ps[(sys.k1,sys.k2)] == sol.ps[(:k1,:k2)] == (0.25, 0.5)
        @test getp(sol, kp)(sol) == getp(sol, sys.kp)(sol) == getp(sol, :kp)(sol) == 1.0
        @test getp(sol, [k1,k2])(sol) == getp(sol, [sys.k1,sys.k2])(sol) == getp(sol, [:k1,:k2])(sol) == [0.25, 0.5]
        @test getp(sol, (k1,k2))(sol) == getp(sol, (sys.k1,sys.k2))(sol) == getp(sol, (:k1,:k2))(sol) == (0.25, 0.5)
    end

    # Handles nonlinear and steady state solutions differently.
    let
        @test_broken false # Currently a problem for nonlinear solutions and steady state solutions (https://github.com/SciML/SciMLBase.jl/issues/720).
        for (sol, sys) in zip(deepcopy([]), []) # zip(deepcopy([nsol, sssol]), [nsys, osys])
            # Get u values.
            @test sol[X] == sol[sys.X] == sol[:X]
            @test sol[XY] == sol[sys.XY][1] == sol[:XY]
            @test sol[[XY,Y]] == sol[[sys.XY,sys.Y]] == sol[[:XY,:Y]]
            @test_broken sol[(XY,Y)] == sol[(sys.XY,sys.Y)] == sol[(:XY,:Y)] # https://github.com/SciML/SciMLBase.jl/issues/710
            @test getu(sol, X)(sol) == getu(sol, sys.X)(sol)[1] == getu(sol, :X)(sol)
            @test getu(sol, XY)(sol) == getu(sol, sys.XY)(sol)[1] == getu(sol, :XY)(sol)
            @test getu(sol, [XY,Y])(sol) == getu(sol, [sys.XY,sys.Y])(sol) == getu(sol, [:XY,:Y])(sol)
            @test_broken getu(sol, (XY,Y))(sol) == getu(sol, (sys.XY,sys.Y))(sol) == getu(sol, (:XY,:Y))(sol)[1] # https://github.com/SciML/SciMLBase.jl/issues/710

            # Get p values.
            @test sol.ps[kp] == sol.ps[sys.kp] == sol.ps[:kp]
            @test sol.ps[[k1,k2]] == sol.ps[[sys.k1,sys.k2]] == sol.ps[[:k1,:k2]]
            @test sol.ps[(k1,k2)] == sol.ps[(sys.k1,sys.k2)] == sol.ps[(:k1,:k2)]
            @test getp(sol, kp)(sol) == getp(sol, sys.kp)(sol) == getp(sol, :kp)(sol)
            @test getp(sol, [k1,k2])(sol) == getp(sol, [sys.k1,sys.k2])(sol) == getp(sol, [:k1,:k2])(sol)
            @test getp(sol, (k1,k2))(sol) == getp(sol, (sys.k1,sys.k2))(sol) == getp(sol, (:k1,:k2))(sol)
        end
    end
end

# Tests plotting.
let 
    for (sol, sys) in zip(deepcopy([osol, ssol, jsol]), [osys, ssys, jsys])
        # Single variable.
        @test length(plot(sol; idxs = X).series_list) == 1
        @test length(plot(sol; idxs = XY).series_list) == 1
        @test length(plot(sol; idxs = sys.X).series_list) == 1
        @test length(plot(sol; idxs = sys.XY).series_list) == 1
        @test length(plot(sol; idxs = :X).series_list) == 1
        @test length(plot(sol; idxs = :XY).series_list) == 1

        # As vector.
        @test length(plot(sol; idxs = [X,Y]).series_list) == 2
        @test length(plot(sol; idxs = [XY,Y]).series_list) == 2
        @test length(plot(sol; idxs = [sys.X,sys.Y]).series_list) == 2
        @test length(plot(sol; idxs = [sys.XY,sys.Y]).series_list) == 2
        @test length(plot(sol; idxs = [:X,:Y]).series_list) == 2
        @test length(plot(sol; idxs = [:XY,:Y]).series_list) == 2

        # As tuple.
        @test length(plot(sol; idxs = (X, Y)).series_list) == 1
        @test length(plot(sol; idxs = (XY, Y)).series_list) == 1
        @test length(plot(sol; idxs = (sys.X, sys.Y)).series_list) == 1
        @test length(plot(sol; idxs = (sys.XY, sys.Y)).series_list) == 1
        @test length(plot(sol; idxs = (:X, :Y)).series_list) == 1
        @test length(plot(sol; idxs = (:XY, :Y)).series_list) == 1
    end     
end


### Mass Action Jump Rate Updating Correctness ###

# Checks that the rates of mass action jumps are correctly updated after parameter values are changed.
let
    # Creates the model.
    @parameters p1 p2
    @variables A(t) B(t) C(t)
    maj = MassActionJump(p1*p2, [A => 1, B => 1], [A => -1, B => -1, C => 1])
    @mtkbuild majsys = JumpSystem([maj], t, [A, B, C], [p1, p2])

    # Creates a JumpProblem and integrator. Checks that the initial mass action rate is correct.
    u0 = [A => 1, B => 2, C => 3]
    ps = [p1 => 3.0, p2 => 2.0]
    dprob = DiscreteProblem(majsys, u0, (0.0, 1.0), ps)
    jprob = JumpProblem(majsys, dprob, Direct())
    jint = init(jprob, SSAStepper())
    @test jprob.massaction_jump.scaled_rates[1] == 6.0

    # Checks that the mass action rate is correctly updated after normal indexing.
    jprob.ps[p1] = 4.0
    @test jprob.massaction_jump.scaled_rates[1] == 8.0
    jprob.ps[majsys.p1] = 5.0
    @test jprob.massaction_jump.scaled_rates[1] == 10.0
    jprob.ps[:p1] = 6.0
    @test jprob.massaction_jump.scaled_rates[1] == 12.0
    setp(jprob, p1)(jprob, 7.0)
    @test jprob.massaction_jump.scaled_rates[1] == 14.0
    setp(jprob, majsys.p1)(jprob, 8.0)
    @test jprob.massaction_jump.scaled_rates[1] == 16.0
    setp(jprob, :p1)(jprob, 3.0)
    @test jprob.massaction_jump.scaled_rates[1] == 6.0

    # Check that the mass action rate is correctly updated when `remake` is used.
    # Checks both when partial and full parameter vectors are provided to `remake`.
    @test remake(jprob; p = [p1 => 4.0]).massaction_jump.scaled_rates[1] == 8.0
    @test remake(jprob; p = [majsys.p1 => 5.0]).massaction_jump.scaled_rates[1] == 10.0
    @test remake(jprob; p = [:p1 => 6.0]).massaction_jump.scaled_rates[1] == 12.0
    @test remake(jprob; p = [p1 => 4.0, p2 => 3.0]).massaction_jump.scaled_rates[1] == 12.0
    @test remake(jprob; p = [majsys.p1 => 5.0, majsys.p2 => 4.0]).massaction_jump.scaled_rates[1] == 20.0
    @test remake(jprob; p = [:p1 => 6.0, :p2 => 5.0]).massaction_jump.scaled_rates[1] == 30.0

    # Checks that updating an integrators parameter values does not affect mass action rate until after
    # `reset_aggregated_jumps!` have been applied as well (wt which point the correct rate is achieved).
    jint.ps[p1] = 4.0
    @test jint.cb.condition.ma_jumps.scaled_rates[1] == 30.0
    reset_aggregated_jumps!(jint)
    @test jint.cb.condition.ma_jumps.scaled_rates[1] == 8.0

    jint.ps[majsys.p1] = 5.0
    @test jint.cb.condition.ma_jumps.scaled_rates[1] == 8.0
    reset_aggregated_jumps!(jint)
    @test jint.cb.condition.ma_jumps.scaled_rates[1] == 10.0

    jint.ps[:p1] = 6.0
    @test jint.cb.condition.ma_jumps.scaled_rates[1] == 10.0
    reset_aggregated_jumps!(jint)
    @test jint.cb.condition.ma_jumps.scaled_rates[1] == 12.0

    setp(jint, p1)(jint, 7.0)
    @test jint.cb.condition.ma_jumps.scaled_rates[1] == 12.0
    reset_aggregated_jumps!(jint)
    @test jint.cb.condition.ma_jumps.scaled_rates[1] == 14.0

    setp(jint, majsys.p1)(jint, 8.0)
    @test jint.cb.condition.ma_jumps.scaled_rates[1] == 14.0
    reset_aggregated_jumps!(jint)
    @test jint.cb.condition.ma_jumps.scaled_rates[1] == 16.0

    setp(jint, :p1)(jint, 3.0)
    @test jint.cb.condition.ma_jumps.scaled_rates[1] == 16.0
    reset_aggregated_jumps!(jint)
    @test jint.cb.condition.ma_jumps.scaled_rates[1] == 6.0
end

