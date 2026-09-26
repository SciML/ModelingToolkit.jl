using ModelingToolkitBase, Symbolics, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D

@testset "Broadcasts over array slices namespace through subsystems" begin
    @variables U(t)[1:2, 1:2, 1:3] I(t)[1:2, 1:3]
    @parameters R[1:2, 1:2, 1:3]
    eqs = [D(U[k, :, :]) ~ -U[k, :, :] ./ R[k, :, :] .+ I for k in 1:2]
    push!(eqs, I ~ U[1, :, :] .+ U[2, :, :])
    @named cell = System(eqs, t, [U, I], [R])
    @named outer = System(Equation[], t; systems = [cell])
    flat = flatten(expand_connections(outer))
    @test length(equations(flat)) == 3
    namespaced = Set(Symbolics.getname.(unknowns(flat)))
    @test namespaced == Set([:cell₊U, :cell₊I])
    rhs = equations(flat)[end].rhs
    variables = Set(
        Symbolics.getname(first(Symbolics.arguments(Symbolics.value(v))))
            for v in Symbolics.get_variables(rhs)
    )
    @test variables == Set([:cell₊U])
end
