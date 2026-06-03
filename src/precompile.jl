include(pkgdir(ModelingToolkitBase, "src", "precompile.jl"))

PrecompileTools.@compile_workload begin
    t = ModelingToolkitBase.t_nounits
    D = ModelingToolkitBase.D_nounits

    function f!(du, u, p)
        du[1] = cos(u[2]) - u[1]
        du[2] = sin(u[1] + u[2]) + u[2]
        du[3] = 2u[4] + u[3] + 1.0
        du[4] = u[5]^2 + u[4]
        du[5] = u[3]^2 + u[5]
        du[6] = u[1] + u[2] + u[3] + u[4] + u[5] + 2.0u[6] + 2.5u[7] + 1.5u[8]
        du[7] = u[1] + u[2] + u[3] + 2.0u[4] + u[5] + 4.0u[6] - 1.5u[7] + 1.5u[8]
        du[8] = u[1] + 2.0u[2] + 3.0u[3] + 5.0u[4] + 6.0u[5] + u[6] - u[7] - u[8]

        du[9] = u[4] + u[5] + u[6] + u[7] + u[8] + 5.3u[9] + 5.8u[10] + 4.8u[11]
        du[10] = u[4] + u[5] + u[6] + 5.3u[7] + u[8] + 7.3u[9] - 4.8u[10] + 4.8u[11]
        du[11] = u[4] + 5.3u[5] + 6.3u[6] + 8.3u[7] + 9.3u[8] + u[9] - u[10] - u[11]

        du[12] = u[7] + u[8] + u[9] + u[10] + u[11] + 8.6u[12] + 8.11u[13] + 7.11u[14]
        du[13] = u[7] + u[8] + u[9] + 8.6u[10] + u[11] + 10.6u[12] - 7.11u[13] + 7.11u[14]
        du[14] = u[7] + 8.6u[8] + 9.6u[9] + 11.6u[10] + 12.6u[11] + u[12] - u[13] - u[14]

        du[15] = u[10] + u[11] + u[12] + u[13] + u[14] + 11.9u[15] + 11.14u[16] + 10.14u[17]
        du[16] = u[10] + u[11] + u[12] + 11.9u[13] + u[14] + 13.9u[15] - 10.14u[16] + 10.14u[17]
        du[17] = u[10] + 11.9u[11] + 12.9u[12] + 14.9u[13] + 15.9u[14] + u[15] - u[16] - u[17]

        du[18] = u[13] + u[14] + u[15] + u[16] + u[17] + 14.12u[18] + 14.17u[19] + 13.17u[20]
        du[19] = u[13] + u[14] + u[15] + 14.12u[16] + u[17] + 16.12u[18] - 13.17u[19] + 13.17u[20]
        du[20] = u[13] + 14.12u[14] + 15.12u[15] + 17.12u[16] + 18.12u[17] + u[18] - u[19] - u[20]
    end
    @variables u[1:20] = rand(20) [irreducible = true]
    eqs = Num[0 for _ in 1:20]
    f!(eqs, u, nothing)
    eqs = 0 .~ eqs
    @mtkcompile model = System(eqs)
    sccprob = SCCNonlinearProblem(model, nothing)
    solve(sccprob, SimpleNonlinearSolve.SimpleNewtonRaphson())
end

