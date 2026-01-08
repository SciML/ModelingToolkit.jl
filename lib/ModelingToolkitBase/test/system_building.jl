using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D

@testset "`state_priorities` and `irreducibles` kwargs take priority over metadata" begin
    @variables x(t) [state_priority = 3] y(t) [irreducible = false]
    @named sys = System([D(x) ~ 2x + y], t; state_priorities = [x => -4], irreducibles = [y])
    @test state_priorities(sys)[x] == -4
    @test y in irreducibles(sys)
end
