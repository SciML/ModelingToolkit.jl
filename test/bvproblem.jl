using BoundaryValueDiffEq, OrdinaryDiffEq, BoundaryValueDiffEqAscher
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

### Test Collocation solvers on simple problems 
solvers = [MIRK4, RadauIIa5, LobattoIIIa3]
daesolvers = [Ascher2, Ascher4, Ascher6]

let
     @parameters α=7.5 β=4.0 γ=8.0 δ=5.0
     @variables x(t)=1.0 y(t)=2.0
     
     eqs = [D(x) ~ α * x - β * x * y,
         D(y) ~ -γ * y + δ * x * y]
     
     u0map = [x => 1.0, y => 2.0]
     parammap = [α => 7.5, β => 4, γ => 8.0, δ => 5.0]
     tspan = (0.0, 10.0)
     
     @mtkbuild lotkavolterra = ODESystem(eqs, t)
     op = ODEProblem(lotkavolterra, u0map, tspan, parammap)
     osol = solve(op, Vern9())
     
     bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(
         lotkavolterra, u0map, tspan, parammap; eval_expression = true)
     
     for solver in solvers
         sol = solve(bvp, solver(), dt = 0.01)
         @test isapprox(sol.u[end], osol.u[end]; atol = 0.01)
         @test sol.u[1] == [1.0, 2.0]
     end
     
     # Test out of place
     bvp2 = SciMLBase.BVProblem{false, SciMLBase.AutoSpecialize}(
         lotkavolterra, u0map, tspan, parammap; eval_expression = true)
     
     for solver in solvers
         sol = solve(bvp2, solver(), dt = 0.01)
         @test isapprox(sol.u[end], osol.u[end]; atol = 0.01)
         @test sol.u[1] == [1.0, 2.0]
     end
end

### Testing on pendulum
let
     @parameters g=9.81 L=1.0
     @variables θ(t) = π / 2
     
     eqs = [D(D(θ)) ~ -(g / L) * sin(θ)]
     
     @mtkbuild pend = ODESystem(eqs, t)
     
     u0map = [θ => π / 2, D(θ) => π / 2]
     parammap = [:L => 1.0, :g => 9.81]
     tspan = (0.0, 6.0)
     
     op = ODEProblem(pend, u0map, tspan, parammap)
     osol = solve(op, Vern9())
     
     bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap)
     for solver in solvers
         sol = solve(bvp, solver(), dt = 0.01)
         @test isapprox(sol.u[end], osol.u[end]; atol = 0.01)
         @test sol.u[1] == [π / 2, π / 2]
     end
     
     # Test out-of-place
     bvp2 = SciMLBase.BVProblem{false, SciMLBase.FullSpecialize}(pend, u0map, tspan, parammap)
     
     for solver in solvers
         sol = solve(bvp2, solver(), dt = 0.01)
         @test isapprox(sol.u[end], osol.u[end]; atol = 0.01)
         @test sol.u[1] == [π / 2, π / 2]
     end
end

###################################################
### ODESystem with Constraint Equations, DAEs with constraints ###
###################################################

# Cartesian pendulum from the docs.
# DAE IVP solved using BoundaryValueDiffEq solvers.
let
    @parameters g
    @variables x(t) y(t) [state_priority = 10] λ(t)
    eqs = [D(D(x)) ~ λ * x
           D(D(y)) ~ λ * y - g
           x^2 + y^2 ~ 1]
    @mtkbuild pend = ODESystem(eqs, t)

    tspan = (0.0, 1.5)
    u0map = [x => 1, y => 0]
    parammap = [g => 1]
    guesses = [λ => 1]

    prob = ODEProblem(pend, u0map, tspan, pmap; guesses)
    sol = solve(prob, Rodas5P())

    bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; guesses)
    
    for solver in solvers
        sol = solve(bvp, solver(), dt = 0.01)
        @test isapprox(sol.u[end], osol.u[end]; atol = 0.01)
        conditions = getfield.(equations(pend)[3:end], :rhs)
        @test [sol[conditions][1]; sol[x][1] - 1; sol[y][1]] ≈ 0 
    end

    bvp2 = SciMLBase.BVProblem{false, SciMLBase.FullSpecialize}(pend, u0map, tspan, parammap)
    for solver in solvers
        sol = solve(bvp, solver(), dt = 0.01)
        @test isapprox(sol.u[end], osol.u[end]; atol = 0.01)
        conditions = getfield.(equations(pend)[3:end], :rhs)
        @test [sol[conditions][1]; sol[x][1] - 1; sol[y][1]] ≈ 0 
    end
end

function test_solvers(solvers, prob, u0map, constraints, equations = []; dt = 0.01)
    for solver in solvers
        sol = solve(bvp, solver(); dt)

        for (k, v) in u0map
            @test sol[k][1] == v
        end
         
        for cons in constraints
            @test sol[cons.rhs - cons.lhs] ≈ 0
        end

        for eq in equations
            @test sol[eq] ≈ 0
        end
    end
end

# Simple ODESystem with BVP constraints.
let
    @parameters α=7.5 β=4.0 γ=8.0 δ=5.0
    @variables x(..) y(t)
    
    eqs = [D(x) ~ α * x - β * x * y,
        D(y) ~ -γ * y + δ * x * y]
    
    u0map = [y => 2.0]
    parammap = [α => 7.5, β => 4, γ => 8.0, δ => 5.0]
    tspan = (0.0, 10.0)
    guesses = [x => 1.0]

    @mtkbuild lotkavolterra = ODESystem(eqs, t)
    op = ODEProblem(lotkavolterra, u0map, tspan, parammap, guesses = guesses)

    constraints = [x(6.) ~ 3]
    bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; guesses, constraints)
    test_solvers(solvers, bvp, u0map, constraints)

    # Testing that more complicated constraints give correct solutions.
    constraints = [y(2.) + x(8.) ~ 12]
    bvp = SciMLBase.BVProblem{false, SciMLBase.FullSpecialize}(pend, u0map, tspan, parammap; guesses, constraints)
    test_solvers(solvers, bvp, u0map, constraints)

    constraints = [α * β - x(6.) ~ 24]
    bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; guesses, constraints)
    test_solvers(solvers, bvp, u0map, constraints)

    # Testing that errors are properly thrown when malformed constraints are given.
    @variables bad(..)
    constraints = [x(1.) + bad(3.) ~ 10]
    @test_throws Exception bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; guesses, constraints)

    constraints = [x(t) + y(t) ~ 3]
    @test_throws Exception bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; guesses, constraints)

    @parameters bad2
    constraints = [bad2 + x(0.) ~ 3]
    @test_throws Exception bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; guesses, constraints)
end

# Adding a midpoint boundary constraint.
# Solve using BVDAE solvers.
let 
    @parameters g
    @variables x(..) y(t) [state_priority = 10] λ(t)
    eqs = [D(D(x(t))) ~ λ * x(t)
           D(D(y)) ~ λ * y - g
           x(t)^2 + y^2 ~ 1]
    @mtkbuild pend = ODESystem(eqs, t)

    tspan = (0.0, 1.5)
    u0map = [x(t) => 0.6, y => 0.8]
    parammap = [g => 1]
    guesses = [λ => 1]

    constraints = [x(0.5) ~ 1]
    bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; constraints, guesses, check_length = false)
    test_solvers(daesolvers, bvp, u0map, constraints)

    bvp2 = SciMLBase.BVProblem{false, SciMLBase.FullSpecialize}(pend, u0map, tspan, parammap)
    test_solvers(daesolvers, bvp2, u0map, constraints, get_alg_eqs(pend))

    # More complicated constraints.
    u0map = [x(t) => 0.6]
    guesses = [λ => 1, y(t) => 0.8]

    constraints = [x(0.5) ~ 1, 
                   x(0.3)^3 + y(0.6)^2 ~ 0.5]
    bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; constraints, guesses, check_length = false)
    test_solvers(daesolvers, bvp, u0map, constraints, get_alg_eqs(pend))

    constraints = [x(0.4) * g ~ y(0.2),
                   y(0.7) ~ 0.3]
    bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; constraints, guesses, check_length = false)
    test_solvers(daesolvers, bvp, u0map, constraints, get_alg_eqs(pend))
end
