using SciCompDSL, Test
using ModelingToolkitBase: getmetadata, t_nounits as t, D_nounits as D

@mtkmodel EmptyModel begin end

@test EmptyModel(; name = :empty) !== nothing

struct Author end

@mtkmodel DecayModel begin
    @metadata begin
        Author = "Test Author"
    end
    @parameters begin
        k = 1.0
    end
    @variables begin
        x(t) = 1.0
    end
    @equations begin
        D(x) ~ -k * x
    end
end

decay = DecayModel(; name = :decay)
@test decay !== nothing
@test getmetadata(decay, Author, nothing) == "Test Author"
