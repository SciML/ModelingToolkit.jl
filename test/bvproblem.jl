### TODO: update when BoundaryValueDiffEqAscher is updated to use the normal boundary condition conventions 
using OrdinaryDiffEq
using BoundaryValueDiffEqMIRK, BoundaryValueDiffEqAscher
using ModelingToolkit
using SciMLBase
using ModelingToolkit: t_nounits as t, D_nounits as D

### Test Collocation solvers on simple problems 
solvers = [MIRK4]
daesolvers = [Ascher2, Ascher4, Ascher6]

@testset "Lotka-Volterra" begin
    @parameters α=7.5 β=4.0 γ=8.0 δ=5.0
    @variables x(t)=1.0 y(t)=2.0

    eqs = [D(x) ~ α * x - β * x * y,
        D(y) ~ -γ * y + δ * x * y]

    u0map = [x => 1.0, y => 2.0]
    parammap = [α => 7.5, β => 4, γ => 8.0, δ => 5.0]
    tspan = (0.0, 10.0)

    @mtkcompile lotkavolterra = System(eqs, t)
    op = ODEProblem(lotkavolterra, [u0map; parammap], tspan)
    osol = solve(op, Vern9())

    bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(
        lotkavolterra, [u0map; parammap], tspan)

    for solver in solvers
        sol = solve(bvp, solver(), dt = 0.01)
        @test isapprox(sol.u[end], osol.u[end]; atol = 0.01)
        @test sol[[x, y], 1] == [1.0, 2.0]
    end

    # Test out of place
    bvp2 = SciMLBase.BVProblem{false, SciMLBase.AutoSpecialize}(
        lotkavolterra, [u0map; parammap], tspan)

    for solver in solvers
        sol = solve(bvp2, solver(), dt = 0.01)
        @test isapprox(sol.u[end], osol.u[end]; atol = 0.01)
        @test sol[[x, y], 1] == [1.0, 2.0]
    end
end

### Testing on pendulum
@testset "Pendulum" begin
    @parameters g=9.81 L=1.0
    @variables θ(t)=π / 2 θ_t(t)

    eqs = [D(θ) ~ θ_t
           D(θ_t) ~ -(g / L) * sin(θ)]

    @mtkcompile pend = System(eqs, t)

    u0map = [θ => π / 2, θ_t => π / 2]
    parammap = [:L => 1.0, :g => 9.81]
    tspan = (0.0, 6.0)

    op = ODEProblem(pend, [u0map; parammap], tspan)
    osol = solve(op, Vern9())

    bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(
        pend, [u0map; parammap], tspan)
    for solver in solvers
        sol = solve(bvp, solver(), dt = 0.01)
        @test isapprox(sol.u[end], osol.u[end]; atol = 0.01)
        @test sol.u[1] == [π / 2, π / 2]
    end

    # Test out-of-place
    bvp2 = SciMLBase.BVProblem{false, SciMLBase.FullSpecialize}(
        pend, [u0map; parammap], tspan)

    for solver in solvers
        sol = solve(bvp2, solver(), dt = 0.01)
        @test isapprox(sol.u[end], osol.u[end]; atol = 0.01)
        @test sol.u[1] == [π / 2, π / 2]
    end
end

##################################################################
### System with constraint equations, DAEs with constraints ###
##################################################################

# Test generation of boundary condition function using `generate_function_bc`. Compare solutions to manually written boundary conditions
@testset "Boundary Condition Compilation" begin
    @parameters α=1.5 β=1.0 γ=3.0 δ=1.0
    @variables x(..) y(..)
    eqs = [D(x(t)) ~ α * x(t) - β * x(t) * y(t),
        D(y(t)) ~ -γ * y(t) + δ * x(t) * y(t)]

    tspan = (0.0, 1.0)
    @mtkcompile lksys = System(eqs, t)

    function lotkavolterra!(du, u, p, t)
        du[1] = p[1] * u[1] - p[2] * u[1] * u[2]
        du[2] = -p[4] * u[2] + p[3] * u[1] * u[2]
    end

    function lotkavolterra(u, p, t)
        [p[1] * u[1] - p[2] * u[1] * u[2], -p[4] * u[2] + p[3] * u[1] * u[2]]
    end

    # Test with a constraint.
    constr = [y(0.5) ~ 2.0]
    @mtkcompile lksys = System(eqs, t; constraints = constr)

    function bc!(resid, u, p, t)
        resid[1] = u(0.0)[1] - 1.0
        resid[2] = u(0.5)[2] - 2.0
    end
    function bc(u, p, t)
        [u(0.0)[1] - 1.0, u(0.5)[2] - 2.0]
    end

    u0 = [1.0, 1.0]
    tspan = (0.0, 1.0)
    p = [1.5, 1.0, 1.0, 3.0]
    bvpi1 = SciMLBase.BVProblem(lotkavolterra!, bc!, u0, tspan, p)
    bvpi2 = SciMLBase.BVProblem(lotkavolterra, bc, u0, tspan, p)
    bvpi3 = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(
        lksys, [x(t) => 1.0], tspan; guesses = [y(t) => 1.0])
    bvpi4 = SciMLBase.BVProblem{false, SciMLBase.FullSpecialize}(
        lksys, [x(t) => 1.0], tspan; guesses = [y(t) => 1.0])

    sol1 = solve(bvpi1, MIRK4(), dt = 0.01)
    sol2 = solve(bvpi2, MIRK4(), dt = 0.01)
    sol3 = solve(bvpi3, MIRK4(), dt = 0.01)
    sol4 = solve(bvpi4, MIRK4(), dt = 0.01)
    @test sol1.t ≈ sol2.t ≈ sol3.t ≈ sol4.t
    @test sol1.u ≈ sol2.u ≈ sol3[[x(t), y(t)]] ≈ sol4[[x(t), y(t)]]
    # @test sol1 ≈ sol2 ≈ sol3 ≈ sol4 # don't get true equality here, not sure why
end

function test_solvers(
        solvers, prob, u0map, constraints, equations = []; dt = 0.05, atol = 1e-2)
    for solver in solvers
        println("Solver: $solver")
        sol = solve(prob, solver(), dt = dt, abstol = atol)
        @test SciMLBase.successful_retcode(sol.retcode)
        p = prob.p
        t = sol.t
        bc = prob.f.bc
        ns = length(prob.u0)
        if isinplace(prob.f)
            resid = zeros(ns)
            bc(resid, sol, p, t)
            @test isapprox(zeros(ns), resid; atol)
            @show resid
        else
            @test isapprox(zeros(ns), bc(sol, p, t); atol)
            @show bc(sol, p, t)
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

# Simple System with BVP constraints.
@testset "ODE with constraints" begin
    @parameters α=1.5 β=1.0 γ=3.0 δ=1.0
    @variables x(..) y(..)

    eqs = [D(x(t)) ~ α * x(t) - β * x(t) * y(t),
        D(y(t)) ~ -γ * y(t) + δ * x(t) * y(t)]

    u0map = []
    tspan = (0.0, 1.0)
    guess = [x(t) => 4.0, y(t) => 2.0]
    constr = [x(0.6) ~ 3.5, x(0.3) ~ 7.0]
    @mtkcompile lksys = System(eqs, t; constraints = constr)

    bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(
        lksys, u0map, tspan; guesses = guess)
    test_solvers(solvers, bvp, u0map, constr; dt = 0.05)

    # Testing that more complicated constraints give correct solutions.
    constr = [y(0.2) + x(0.8) ~ 3.0, y(0.3) ~ 2.0]
    @mtkcompile lksys = System(eqs, t; constraints = constr)
    bvp = SciMLBase.BVProblem{false, SciMLBase.FullSpecialize}(
        lksys, u0map, tspan; guesses = guess)
    test_solvers(solvers, bvp, u0map, constr; dt = 0.05)

    constr = [α * β - x(0.6) ~ 0.0, y(0.2) ~ 3.0]
    @mtkcompile lksys = System(eqs, t; constraints = constr)
    bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(
        lksys, u0map, tspan; guesses = guess)
    test_solvers(solvers, bvp, u0map, constr)
end

# Cartesian pendulum from the docs.
# DAE IVP solved using BoundaryValueDiffEq solvers.
# let
#     @parameters g
#     @variables x(t) y(t) [state_priority = 10] λ(t)
#     eqs = [D(D(x)) ~ λ * x
#            D(D(y)) ~ λ * y - g
#            x^2 + y^2 ~ 1]
#     @mtkcompile pend = System(eqs, t)
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
#     constr = [x(0.5) ~ 1]
#     @mtkcompile pend = System(eqs, t; constr)
# 
#     tspan = (0.0, 1.5)
#     u0map = [x(t) => 0.6, y => 0.8]
#     parammap = [g => 1]
#     guesses = [λ => 1]
# 
#     bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; guesses, check_length = false)
#     test_solvers(daesolvers, bvp, u0map, constr)
# 
#     bvp2 = SciMLBase.BVProblem{false, SciMLBase.FullSpecialize}(pend, u0map, tspan, parammap)
#     test_solvers(daesolvers, bvp2, u0map, constr, get_alg_eqs(pend))
# 
#     # More complicated constr.
#     u0map = [x(t) => 0.6]
#     guesses = [λ => 1, y(t) => 0.8]
# 
#     constr = [x(0.5) ~ 1, 
#                    x(0.3)^3 + y(0.6)^2 ~ 0.5]
#     @mtkcompile pend = System(eqs, t; constr)
#     bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; guesses, check_length = false)
#     test_solvers(daesolvers, bvp, u0map, constr, get_alg_eqs(pend))
# 
#     constr = [x(0.4) * g ~ y(0.2),
#                    y(0.7) ~ 0.3]
#     @mtkcompile pend = System(eqs, t; constr)
#     bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; guesses, check_length = false)
#     test_solvers(daesolvers, bvp, u0map, constr, get_alg_eqs(pend))
# end

@testset "Cost function compilation" begin
    @parameters α=1.5 β=1.0 γ=3.0 δ=1.0
    @variables x(..) y(..)
    t = ModelingToolkit.t_nounits

    eqs = [D(x(t)) ~ α * x(t) - β * x(t) * y(t),
        D(y(t)) ~ -γ * y(t) + δ * x(t) * y(t)]

    tspan = (0.0, 1.0)
    u0map = [x(t) => 4.0, y(t) => 2.0]
    parammap = [α => 7.5, β => 4, γ => 8.0, δ => 5.0]
    costs = [x(0.6), x(0.3)^2]
    consolidate(u, sub) = (u[1] + 3)^2 + u[2] + sum(sub; init = 0)
    @mtkcompile lksys = System(eqs, t; costs, consolidate)

    @test_throws ModelingToolkit.SystemCompatibilityError ODEProblem(
        lksys, [u0map; parammap], tspan)
    prob = ODEProblem(lksys, [u0map; parammap], tspan; check_compatibility = false)
    sol = solve(prob, Tsit5())
    costfn = ModelingToolkit.generate_cost(
        lksys; expression = Val{false}, wrap_gfw = Val{true})
    _t = tspan[2]
    @test costfn(sol, prob.p, _t) ≈ (sol(0.6; idxs = x(t)) + 3)^2 + sol(0.3; idxs = x(t))^2

    ### With a parameter
    @parameters t_c
    costs = [y(t_c) + x(0.0), x(0.4)^2]
    consolidate(u, sub) = log(u[1]) - u[2] + sum(sub; init = 0)
    @mtkcompile lksys = System(eqs, t; costs, consolidate)
    @test t_c ∈ Set(parameters(lksys))
    push!(parammap, t_c => 0.56)
    prob = ODEProblem(lksys, [u0map; parammap], tspan; check_compatibility = false)
    sol = solve(prob, Tsit5())
    costfn = ModelingToolkit.generate_cost(
        lksys; expression = Val{false}, wrap_gfw = Val{true})
    @test costfn(sol, prob.p, _t) ≈
          log(sol(0.56; idxs = y(t)) + sol(0.0; idxs = x(t))) - sol(0.4; idxs = x(t))^2
end
