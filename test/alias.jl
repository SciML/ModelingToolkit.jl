using Test
using ModelingToolkit: AliasGraph

ag = AliasGraph(10)
ag[1] = 1 => 2
ag[2] = -1 => 3
ag[4] = -1 => 1
ag[5] = -1 => 4
for _ in 1:5 # check ag is robust
    @test ag[1] == (-1, 3)
    @test ag[2] == (-1, 3)
    @test ag[4] == (1, 3)
    @test ag[5] == (-1, 3)
end

@test 1 in keys(ag)
@test 2 in keys(ag)
@test !(3 in keys(ag))
@test 4 in keys(ag)
@test 5 in keys(ag)
