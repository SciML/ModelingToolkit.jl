include("shared/mtktestset.jl")

@safetestset "Linearization Dummy Derivative Tests" include("downstream/linearization_dd.jl")
@safetestset "Inverse Models Test" include("downstream/inversemodel.jl")
@safetestset "Disturbance model Test" include("downstream/test_disturbance_model.jl")
