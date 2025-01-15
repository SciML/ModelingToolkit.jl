### TODO: update when BoundaryValueDiffEqAscher is updated to use the normal boundary condition conventions 

using BoundaryValueDiffEq, OrdinaryDiffEq, BoundaryValueDiffEqAscher
using ModelingToolkit
using SciMLBase
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit: process_constraints

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
     @variables θ(t) = π / 2 θ_t(t)
     
     eqs = [D(θ) ~ θ_t
            D(θ_t) ~ -(g / L) * sin(θ)]
     
     @mtkbuild pend = ODESystem(eqs, t)
     
     u0map = [θ => π / 2, θ_t => π / 2]
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

##################################################################
### ODESystem with constraint equations, DAEs with constraints ###
##################################################################

# Test generation of boundary condition function using `process_constraints`. Compare solutions to manually written boundary conditions
let
    @parameters α=1.5 β=1.0 γ=3.0 δ=1.0
    @variables x(..) y(..)
    eqs = [D(x(t)) ~ α * x(t) - β * x(t) * y(t),
           D(y(t)) ~ -γ * y(t) + δ * x(t) * y(t)]
    
    tspan = (0., 1.)
    @mtkbuild lksys = ODESystem(eqs, t)

    function lotkavolterra!(du, u, p, t) 
        du[1] = p[1]*u[1] - p[2]*u[1]*u[2]
        du[2] = -p[4]*u[2] + p[3]*u[1]*u[2]
    end

    function lotkavolterra(u, p, t) 
        [p[1]*u[1] - p[2]*u[1]*u[2], -p[4]*u[2] + p[3]*u[1]*u[2]]
    end
    # Compare the built bc function to the actual constructed one.
    function bc!(resid, u, p, t) 
        resid[1] = u[1][1] - 1.
        resid[2] = u[1][2] - 2.
        nothing
    end
    function bc(u, p, t)
        [u[1][1] - 1., u[1][2] - 2.]
    end

    u0 = [1., 2.]; p = [1.5, 1., 1., 3.]
    genbc_iip = ModelingToolkit.process_constraints(lksys, nothing, u0, [1, 2], tspan, true)
    genbc_oop = ModelingToolkit.process_constraints(lksys, nothing, u0, [1, 2], tspan, false)

    bvpi1 = SciMLBase.BVProblem(lotkavolterra!, bc!, [1.,2.], tspan, p)
    bvpi2 = SciMLBase.BVProblem(lotkavolterra!, genbc_iip, [1.,2.], tspan, p)

    sol1 = solve(bvpi1, MIRK4(), dt = 0.05)
    sol2 = solve(bvpi2, MIRK4(), dt = 0.05)
    @test sol1 ≈ sol2

    bvpo1 = BVProblem(lotkavolterra, bc, [1,2], tspan, p)
    bvpo2 = BVProblem(lotkavolterra, genbc_oop, [1,2], tspan, p)

    sol1 = solve(bvpo1, MIRK4(), dt = 0.05)
    sol2 = solve(bvpo2, MIRK4(), dt = 0.05)
    @test sol1 ≈ sol2

    # Test with a constraint.
    constraints = [y(0.5) ~ 2.]

    function bc!(resid, u, p, t) 
        resid[1] = u(0.0)[1] - 1.
        resid[2] = u(0.5)[2] - 2.
    end
    function bc(u, p, t)
        [u(0.0)[1] - 1., u(0.5)[2] - 2.]
    end

    u0 = [1, 1.]
    genbc_iip = ModelingToolkit.process_constraints(lksys, constraints, u0, [1], tspan, true)
    genbc_oop = ModelingToolkit.process_constraints(lksys, constraints, u0, [1], tspan, false)

    bvpi1 = SciMLBase.BVProblem(lotkavolterra!, bc!, u0, tspan, p)
    bvpi2 = SciMLBase.BVProblem(lotkavolterra!, genbc_iip, u0, tspan, p)
    bvpi3 = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(lksys, [x(t) => 1.], tspan; guesses = [y(t) => 1.], constraints)
    bvpi4 = SciMLBase.BVProblem{true, SciMLBase.FullSpecialize}(lksys, [x(t) => 1.], tspan; guesses = [y(t) => 1.], constraints)
    
    @btime sol1 = solve(bvpi1, MIRK4(), dt = 0.01)
    @btime sol2 = solve(bvpi2, MIRK4(), dt = 0.01)
    @btime sol3 = solve(bvpi3, MIRK4(), dt = 0.01)
    @btime sol4 = solve(bvpi4, MIRK4(), dt = 0.01)
    @test sol1 ≈ sol2 ≈ sol3 ≈ sol4 # don't get true equality here, not sure why

    bvpo1 = BVProblem(lotkavolterra, bc, [1,2], tspan, p)
    bvpo2 = BVProblem(lotkavolterra, genbc_oop, [1,2], tspan, p)
    bvpo3 = SciMLBase.BVProblem{false, SciMLBase.FullSpecialize}(lksys, [x(t) => 1.], tspan, parammap; guesses = [y(t) => 1.], constraints)

    @btime sol1 = solve(bvpo1, MIRK4(), dt = 0.05)
    @btime sol2 = solve(bvpo2, MIRK4(), dt = 0.05)
    @btime sol3 = solve(bvpo3, MIRK4(), dt = 0.05)
    @test sol1 ≈ sol2 ≈ sol3
end

function test_solvers(solvers, prob, u0map, constraints, equations = []; dt = 0.05, atol = 1e-4)
    for solver in solvers
        println("Solver: $solver")
        @btime sol = solve(prob, solver(), dt = dt, abstol = atol)
        @test SciMLBase.successful_retcode(sol.retcode)
        p = prob.p; t = sol.t; bc = prob.f.bc
        ns = length(prob.u0)
        if isinplace(bvp.f)
            resid = zeros(ns)
            bc!(resid, sol, p, t)
            @test isapprox(zeros(ns), resid; atol)
        else
            @test isapprox(zeros(ns), bc(sol, p, t); atol)
        end

        for (k, v) in u0map
            @test sol[k][1] == v
        end
         
        # for cons in constraints
        #     @test sol[cons.rhs - cons.lhs] ≈ 0
        # end

        for eq in equations
            @test sol[eq] ≈ 0
        end
    end
end

# Simple ODESystem with BVP constraints.
let
    @parameters α=1.5 β=1.0 γ=3.0 δ=1.0
    @variables x(..) y(..)
    
    eqs = [D(x(t)) ~ α * x(t) - β * x(t) * y(t),
           D(y(t)) ~ -γ * y(t) + δ * x(t) * y(t)]
    
    u0map = []
    tspan = (0.0, 10.0)
    guesses = [x(t) => 4.0, y(t) => 2.]

    @mtkbuild lksys = ODESystem(eqs, t)

    constraints = [x(6.) ~ 3.5, x(3.) ~ 7.]
    bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(lksys, u0map, tspan; guesses, constraints)
    test_solvers(solvers, bvp, u0map, constraints; dt = 0.1)

    # Testing that more complicated constraints give correct solutions.
    constraints = [y(2.) + x(8.) ~ 2.]
    bvp = SciMLBase.BVProblem{false, SciMLBase.FullSpecialize}(lotkavolterra, u0map, tspan, parammap; guesses, constraints)
    test_solvers(solvers, bvp, u0map, constraints)

    constraints = [α * β - x(6.) ~ 0.5]
    bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(lotkavolterra, u0map, tspan, parammap; guesses, constraints)
    test_solvers(solvers, bvp, u0map, constraints)

    # Testing that errors are properly thrown when malformed constraints are given.
    @variables bad(..)
    constraints = [x(1.) + bad(3.) ~ 10]
    @test_throws Exception bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(lotkavolterra, u0map, tspan, parammap; guesses, constraints)

    constraints = [x(t) + y(t) ~ 3]
    @test_throws Exception bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(lotkavolterra, u0map, tspan, parammap; guesses, constraints)

    @parameters bad2
    constraints = [bad2 + x(0.) ~ 3]
    @test_throws Exception bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(lotkavolterra, u0map, tspan, parammap; guesses, constraints)
end

# Cartesian pendulum from the docs.
# DAE IVP solved using BoundaryValueDiffEq solvers.
# let
#     @parameters g
#     @variables x(t) y(t) [state_priority = 10] λ(t)
#     eqs = [D(D(x)) ~ λ * x
#            D(D(y)) ~ λ * y - g
#            x^2 + y^2 ~ 1]
#     @mtkbuild pend = ODESystem(eqs, t)
# 
#     tspan = (0.0, 1.5)
#     u0map = [x => 1, y => 0]
#     pmap = [g => 1]
#     guess = [λ => 1]
# 
#     prob = ODEProblem(pend, u0map, tspan, pmap; guesses = guess)
#     osol = solve(prob, Rodas5P())
# 
#     zeta = [0., 0., 0., 0., 0.]
#     bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; guesses = guess)
#     
#     for solver in solvers
#         sol = solve(bvp, solver(zeta), dt = 0.001)
#         @test isapprox(sol.u[end], osol.u[end]; atol = 0.01)
#         conditions = getfield.(equations(pend)[3:end], :rhs)
#         @test isapprox([sol[conditions][1]; sol[x][1] - 1; sol[y][1]], zeros(5), atol = 0.001)
#     end
# 
#     bvp2 = SciMLBase.BVProblem{false, SciMLBase.FullSpecialize}(pend, u0map, tspan, parammap)
#     for solver in solvers
#         sol = solve(bvp, solver(zeta), dt = 0.01)
#         @test isapprox(sol.u[end], osol.u[end]; atol = 0.01)
#         conditions = getfield.(equations(pend)[3:end], :rhs)
#         @test [sol[conditions][1]; sol[x][1] - 1; sol[y][1]] ≈ 0 
#     end
# end

# Adding a midpoint boundary constraint.
# Solve using BVDAE solvers.
# let 
#     @parameters g
#     @variables x(..) y(t) [state_priority = 10] λ(t)
#     eqs = [D(D(x(t))) ~ λ * x(t)
#            D(D(y)) ~ λ * y - g
#            x(t)^2 + y^2 ~ 1]
#     @mtkbuild pend = ODESystem(eqs, t)
# 
#     tspan = (0.0, 1.5)
#     u0map = [x(t) => 0.6, y => 0.8]
#     parammap = [g => 1]
#     guesses = [λ => 1]
# 
#     constraints = [x(0.5) ~ 1]
#     bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; constraints, guesses, check_length = false)
#     test_solvers(daesolvers, bvp, u0map, constraints)
# 
#     bvp2 = SciMLBase.BVProblem{false, SciMLBase.FullSpecialize}(pend, u0map, tspan, parammap)
#     test_solvers(daesolvers, bvp2, u0map, constraints, get_alg_eqs(pend))
# 
#     # More complicated constraints.
#     u0map = [x(t) => 0.6]
#     guesses = [λ => 1, y(t) => 0.8]
# 
#     constraints = [x(0.5) ~ 1, 
#                    x(0.3)^3 + y(0.6)^2 ~ 0.5]
#     bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; constraints, guesses, check_length = false)
#     test_solvers(daesolvers, bvp, u0map, constraints, get_alg_eqs(pend))
# 
#     constraints = [x(0.4) * g ~ y(0.2),
#                    y(0.7) ~ 0.3]
#     bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; constraints, guesses, check_length = false)
#     test_solvers(daesolvers, bvp, u0map, constraints, get_alg_eqs(pend))
# end
