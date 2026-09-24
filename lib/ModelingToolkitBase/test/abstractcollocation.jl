using ModelingToolkitBase, Test

const PUBLIC_COLLOCATION_SELECTORS = (
    :JuMPCollocation,
    :InfiniteOptCollocation,
    :CasADiCollocation,
    :PyomoCollocation,
)
const INTERNAL_COLLOCATION_HOOKS = (
    :prepare_and_optimize!,
    :get_t_values,
    :get_U_values,
    :get_V_values,
    :get_P_values,
    :successful_solve,
)

@testset "AbstractCollocation public boundary" begin
    public_names = names(ModelingToolkitBase)
    @test :AbstractCollocation ∈ public_names
    @test isabstracttype(ModelingToolkitBase.AbstractCollocation)

    abstract_collocation_docs = sprint(
        show, MIME"text/plain"(), Base.Docs.doc(ModelingToolkitBase.AbstractCollocation)
    )
    @test occursin("not a third-party extension interface", abstract_collocation_docs)
    @test occursin("Do not subtype AbstractCollocation", abstract_collocation_docs)
    @test occursin("solve(prob, solver)", abstract_collocation_docs)

    for selector in PUBLIC_COLLOCATION_SELECTORS
        @test selector ∈ public_names
        @test Base.Docs.doc(getproperty(ModelingToolkitBase, selector)) !== nothing
    end

    for hook in INTERNAL_COLLOCATION_HOOKS
        @test hook ∉ public_names
    end
end
