using ModelingToolkit
using NonlinearSolve, SCCNonlinearSolve
using OrdinaryDiffEq
using SciMLBase, Symbolics
using LinearAlgebra, Test
using ModelingToolkit: t_nounits as t, D_nounits as D

@testset "Trivial case" begin
    function f!(du, u, p)
        du[1] = cos(u[2]) - u[1]
        du[2] = sin(u[1] + u[2]) + u[2]
        du[3] = 2u[4] + u[3] + 1.0
        du[4] = u[5]^2 + u[4]
        du[5] = u[3]^2 + u[5]
        du[6] = u[1] + u[2] + u[3] + u[4] + u[5] + 2.0u[6] + 2.5u[7] + 1.5u[8]
        du[7] = u[1] + u[2] + u[3] + 2.0u[4] + u[5] + 4.0u[6] - 1.5u[7] + 1.5u[8]
        du[8] = u[1] + 2.0u[2] + 3.0u[3] + 5.0u[4] + 6.0u[5] + u[6] - u[7] - u[8]
    end
    @variables u[1:8] [irreducible = true]
    eqs = Any[0 for _ in 1:8]
    f!(eqs, u, nothing)
    eqs = 0 .~ eqs
    @named model = NonlinearSystem(eqs)
    @test_throws ["simplified", "required"] SCCNonlinearProblem(model, [])
    _model = structural_simplify(model; split = false)
    @test_throws ["not compatible"] SCCNonlinearProblem(_model, [])
    model = structural_simplify(model)
    prob = NonlinearProblem(model, [u => zeros(8)])
    sccprob = SCCNonlinearProblem(model, [u => zeros(8)])
    sol1 = solve(prob, NewtonRaphson())
    sol2 = solve(sccprob, NewtonRaphson())
    @test SciMLBase.successful_retcode(sol1)
    @test SciMLBase.successful_retcode(sol2)
    @test sol1[u] ≈ sol2[u]
end

@testset "With parameters" begin
    function f!(du, u, (p1, p2), t)
        x = (*)(p1[4], u[1])
        y = (*)(p1[4], (+)(0.1016, (*)(-1, u[1])))
        z1 = ifelse((<)(p2[1], 0),
            (*)((*)(457896.07999999996, p1[2]), sqrt((*)(1.1686468413521012e-5, p1[3]))),
            0)
        z2 = ifelse((>)(p2[1], 0),
            (*)((*)((*)(0.58, p1[2]), sqrt((*)(1 // 86100, p1[3]))), u[4]),
            0)
        z3 = ifelse((>)(p2[1], 0),
            (*)((*)(457896.07999999996, p1[2]), sqrt((*)(1.1686468413521012e-5, p1[3]))),
            0)
        z4 = ifelse((<)(p2[1], 0),
            (*)((*)((*)(0.58, p1[2]), sqrt((*)(1 // 86100, p1[3]))), u[5]),
            0)
        du[1] = p2[1]
        du[2] = (+)(z1, (*)(-1, z2))
        du[3] = (+)(z3, (*)(-1, z4))
        du[4] = (+)((*)(-1, u[2]), (*)((*)(1 // 86100, y), u[4]))
        du[5] = (+)((*)(-1, u[3]), (*)((*)(1 // 86100, x), u[5]))
    end
    p = (
        [0.04864391799335977, 7.853981633974484e-5, 1.4034843205574914,
            0.018241469247509915, 300237.05, 9.226186337232914],
        [0.0508])
    u0 = [0.0, 0.0, 0.0, 789476.0, 101325.0]
    tspan = (0.0, 1.0)
    mass_matrix = [1.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0;
                   0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0]
    dt = 1e-3
    function nlf(u1, (u0, p))
        resid = Any[0 for _ in u0]
        f!(resid, u1, p, 0.0)
        return mass_matrix * (u1 - u0) - dt * resid
    end

    prob = NonlinearProblem(nlf, u0, (u0, p))
    @test_throws Exception solve(prob, SimpleNewtonRaphson(), abstol = 1e-9)
    sol = solve(prob, TrustRegion(); abstol = 1e-9)

    @variables u[1:5] [irreducible = true]
    @parameters p1[1:6] p2
    eqs = 0 .~ collect(nlf(u, (u0, (p1, p2))))
    @mtkbuild sys = NonlinearSystem(eqs, [u], [p1, p2])
    sccprob = SCCNonlinearProblem(sys, [u => u0], [p1 => p[1], p2 => p[2][]])
    sccsol = solve(sccprob, SimpleNewtonRaphson(); abstol = 1e-9)
    @test SciMLBase.successful_retcode(sccsol)
    @test norm(sccsol.resid) < norm(sol.resid)
end

@testset "Transistor amplifier" begin
    C = [k * 1e-6 for k in 1:5]
    Ub = 6
    UF = 0.026
    α = 0.99
    β = 1e-6
    R0 = 1000
    R = 9000
    Ue(t) = 0.1 * sin(200 * π * t)

    function transamp(out, du, u, p, t)
        g(x) = 1e-6 * (exp(x / 0.026) - 1)
        y1, y2, y3, y4, y5, y6, y7, y8 = u
        out[1] = -Ue(t) / R0 + y1 / R0 + C[1] * du[1] - C[1] * du[2]
        out[2] = -Ub / R + y2 * 2 / R - (α - 1) * g(y2 - y3) - C[1] * du[1] + C[1] * du[2]
        out[3] = -g(y2 - y3) + y3 / R + C[2] * du[3]
        out[4] = -Ub / R + y4 / R + α * g(y2 - y3) + C[3] * du[4] - C[3] * du[5]
        out[5] = -Ub / R + y5 * 2 / R - (α - 1) * g(y5 - y6) - C[3] * du[4] + C[3] * du[5]
        out[6] = -g(y5 - y6) + y6 / R + C[4] * du[6]
        out[7] = -Ub / R + y7 / R + α * g(y5 - y6) + C[5] * du[7] - C[5] * du[8]
        out[8] = y8 / R - C[5] * du[7] + C[5] * du[8]
    end

    u0 = [0, Ub / 2, Ub / 2, Ub, Ub / 2, Ub / 2, Ub, 0]
    du0 = [
        51.338775,
        51.338775,
        -Ub / (2 * (C[2] * R)),
        -24.9757667,
        -24.9757667,
        -Ub / (2 * (C[4] * R)),
        -10.00564453,
        -10.00564453
    ]
    daeprob = DAEProblem(transamp, du0, u0, (0.0, 0.1))
    daesol = solve(daeprob, DImplicitEuler())

    t0 = daesol.t[5]
    t1 = daesol.t[6]
    u0 = daesol.u[5]
    u1 = daesol.u[6]
    dt = t1 - t0

    @variables y(t)[1:8]
    eqs = Any[0 for _ in 1:8]
    transamp(eqs, collect(D(y)), y, nothing, t)
    eqs = 0 .~ eqs
    subrules = Dict(Symbolics.unwrap(D(y[i])) => ((y[i] - u0[i]) / dt) for i in 1:8)
    eqs = substitute.(eqs, (subrules,))
    @mtkbuild sys = NonlinearSystem(eqs)
    prob = NonlinearProblem(sys, [y => u0], [t => t0])
    sol = solve(prob, NewtonRaphson(); abstol = 1e-12)

    sccprob = SCCNonlinearProblem(sys, [y => u0], [t => t0])
    sccsol = solve(sccprob, NewtonRaphson(); abstol = 1e-12)

    @test sol.u≈sccsol.u atol=1e-10
end

@testset "Expression caching" begin
    @variables x[1:4] = rand(4)
    val = Ref(0)
    function func(x, y)
        val[] += 1
        x + y
    end
    @register_symbolic func(x, y)
    @mtkbuild sys = NonlinearSystem([0 ~ x[1]^3 + x[2]^3 - 5
                                     0 ~ sin(x[1] - x[2]) - 0.5
                                     0 ~ func(x[1], x[2]) * exp(x[3]) - x[4]^3 - 5
                                     0 ~ func(x[1], x[2]) * exp(x[4]) - x[3]^3 - 4])
    sccprob = SCCNonlinearProblem(sys, [])
    sccsol = solve(sccprob, NewtonRaphson())
    @test SciMLBase.successful_retcode(sccsol)
    @test val[] == 1
end

import ModelingToolkitStandardLibrary.Blocks as B
import ModelingToolkitStandardLibrary.Mechanical.Translational as T
import ModelingToolkitStandardLibrary.Hydraulic.IsothermalCompressible as IC

@testset "Caching of subexpressions of different types" begin
    liquid_pressure(rho, rho_0, bulk) = (rho / rho_0 - 1) * bulk
    gas_pressure(rho, rho_0, p_gas, rho_gas) = rho * ((0 - p_gas) / (rho_0 - rho_gas))
    full_pressure(rho, rho_0, bulk, p_gas, rho_gas) = ifelse(
        rho >= rho_0, liquid_pressure(rho, rho_0, bulk),
        gas_pressure(rho, rho_0, p_gas, rho_gas))

    @component function Volume(;
            #parameters
            area,
            direction = +1,
            x_int,
            name)
        pars = @parameters begin
            area = area
            x_int = x_int
            rho_0 = 1000
            bulk = 1e9
            p_gas = -1000
            rho_gas = 1
        end

        vars = @variables begin
            x(t) = x_int
            dx(t), [guess = 0]
            p(t), [guess = 0]
            f(t), [guess = 0]
            rho(t), [guess = 0]
            m(t), [guess = 0]
            dm(t), [guess = 0]
        end

        systems = @named begin
            port = IC.HydraulicPort()
            flange = T.MechanicalPort()
        end

        eqs = [
               # connectors
               port.p ~ p
               port.dm ~ dm
               flange.v * direction ~ dx
               flange.f * direction ~ -f

               # differentials
               D(x) ~ dx
               D(m) ~ dm

               # physics
               p ~ full_pressure(rho, rho_0, bulk, p_gas, rho_gas)
               f ~ p * area
               m ~ rho * x * area]

        return ODESystem(eqs, t, vars, pars; name, systems)
    end

    systems = @named begin
        fluid = IC.HydraulicFluid(; bulk_modulus = 1e9)

        src1 = IC.Pressure(;)
        src2 = IC.Pressure(;)

        vol1 = Volume(; area = 0.01, direction = +1, x_int = 0.1)
        vol2 = Volume(; area = 0.01, direction = +1, x_int = 0.1)

        mass = T.Mass(; m = 10)

        sin1 = B.Sine(; frequency = 0.5, amplitude = +0.5e5, offset = 10e5)
        sin2 = B.Sine(; frequency = 0.5, amplitude = -0.5e5, offset = 10e5)
    end

    eqs = [connect(fluid, src1.port)
           connect(fluid, src2.port)
           connect(src1.port, vol1.port)
           connect(src2.port, vol2.port)
           connect(vol1.flange, mass.flange, vol2.flange)
           connect(src1.p, sin1.output)
           connect(src2.p, sin2.output)]

    initialization_eqs = [mass.s ~ 0.0
                          mass.v ~ 0.0]

    @mtkbuild sys = ODESystem(eqs, t, [], []; systems, initialization_eqs)
    prob = ODEProblem(sys, [], (0, 5))
    sol = solve(prob)
    @test SciMLBase.successful_retcode(sol)
end
