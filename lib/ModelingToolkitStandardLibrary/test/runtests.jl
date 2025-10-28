using SafeTestsets, Test

const GROUP = get(ENV, "GROUP", "All")

@time begin
    if GROUP == "QA" || GROUP == "All"
        @time @safetestset "Aqua" begin
            include("aqua.jl")
        end
    end

    if GROUP == "Core" || GROUP == "All"
        @testset "Core" begin
            # Blocks
            @safetestset "Blocks: utils" begin
                include("Blocks/utils.jl")
            end
            @safetestset "Blocks: math" begin
                include("Blocks/math.jl")
            end
            @safetestset "Blocks: nonlinear" begin
                include("Blocks/nonlinear.jl")
            end
            @safetestset "Blocks: continuous" begin
                include("Blocks/continuous.jl")
            end
            @safetestset "Blocks: sources" begin
                include("Blocks/sources.jl")
            end
            @safetestset "Blocks: analysis points" begin
                include("Blocks/test_analysis_points.jl")
            end

            # Electrical
            @safetestset "Analog Circuits" begin
                include("Electrical/analog.jl")
            end

            @safetestset "Digital Circuits" begin
                include("Electrical/digital.jl")
            end
            @safetestset "Chua Circuit Demo" begin
                include("chua_circuit.jl")
            end

            # Thermal
            @safetestset "Thermal Circuits" begin
                include("Thermal/thermal.jl")
            end
            @safetestset "Thermal Demo" begin
                include("Thermal/demo.jl")
                include("Thermal/motor.jl")
                include("Thermal/piston.jl")
            end

            # Magnetic
            @safetestset "Magnetic" begin
                include("Magnetic/magnetic.jl")
            end

            # Mechanical
            @safetestset "Mechanical Rotation" begin
                include("Mechanical/rotational.jl")
            end
            @safetestset "Mechanical Translation" begin
                include("Mechanical/translational.jl")
            end
            @safetestset "Mechanical Translation Modelica" begin
                include("Mechanical/translational_modelica.jl")
            end
            @safetestset "Multi-Domain" begin
                include("multi_domain.jl")
            end

            # Hydraulic
            @safetestset "Hydraulic IsothermalCompressible" begin
                include("Hydraulic/isothermal_compressible.jl")
            end
        end
    end
end
