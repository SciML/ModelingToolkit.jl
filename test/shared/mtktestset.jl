using SafeTestsets

const MTKBasePath = joinpath(dirname(dirname(@__DIR__)), "lib", "ModelingToolkitBase")

macro mtktestset(name, file)
    return quote
        @safetestset $name begin
            using ModelingToolkit
            import ModelingToolkitBase
            include(joinpath(pkgdir(ModelingToolkitBase), "test", $file))
        end
    end
end
