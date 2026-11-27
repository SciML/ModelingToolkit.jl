using Test, FMI, FMIZoo
using OrdinaryDiffEq: Tsit5

# Later releases name their archives differently than FMIZoo expects.
const REFERENCE_FMUS = "0.0.39"

function reference_fmu(model, version; type)
    return loadFMU(model, "ModelicaReferenceFMUs", REFERENCE_FMUS, "$version.0"; type)
end

@testset "Reference FMU Model Exchange simulation" begin
    @testset "Dahlquist v$version at t = $time" for version in (2, 3), time in (0.5, 1.0)
        fmu = reference_fmu("Dahlquist", version; type = :ME)
        sol = simulateME(
            fmu, (0.0, time); solver = Tsit5(), recordValues = ["x"],
            saveat = [time], showProgress = false
        )
        @test sol.success
    end
end
