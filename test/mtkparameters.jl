using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D, MTKParameters
using SymbolicIndexingInterface
using SciMLStructures: SciMLStructures, canonicalize, Tunable, Discrete, Constants

@parameters a b c d::Integer e[1:3] f[1:3, 1:3]::Int g::Vector{AbstractFloat} h::String
@named sys = ODESystem(
    Equation[], t, [], [a, c, d, e, f, g, h], parameter_dependencies = [b => 2a],
    continuous_events = [[a ~ 0] => [c ~ 0]], defaults = Dict(a => 0.0))
sys = complete(sys)

ivs = Dict(c => 3a, d => 4, e => [5.0, 6.0, 7.0],
    f => ones(Int, 3, 3), g => [0.1, 0.2, 0.3], h => "foo")

ps = MTKParameters(sys, ivs)
@test_nowarn copy(ps)
# dependent initialization, also using defaults
@test getp(sys, a)(ps) == getp(sys, b)(ps) == getp(sys, c)(ps) == 0.0
@test getp(sys, d)(ps) isa Int

ivs[a] = 1.0
ps = MTKParameters(sys, ivs)
@test_broken getp(sys, g) # SII bug
for (p, val) in ivs
    isequal(p, g) && continue # broken
    if isequal(p, c)
        val = 3ivs[a]
    end
    idx = parameter_index(sys, p)
    # ensure getindex with `ParameterIndex` works
    @test ps[idx] == getp(sys, p)(ps) == val
end

# ensure setindex! with `ParameterIndex` works
ps[parameter_index(sys, a)] = 3.0
@test getp(sys, a)(ps) == 3.0
setp(sys, a)(ps, 1.0)

@test getp(sys, a)(ps) == getp(sys, b)(ps) / 2 == getp(sys, c)(ps) / 3 == 1.0

for (portion, values) in [(Tunable(), vcat(ones(9), [1.0, 4.0, 5.0, 6.0, 7.0]))
                          (Discrete(), [3.0])
                          (Constants(), [0.1, 0.2, 0.3])]
    buffer, repack, alias = canonicalize(portion, ps)
    @test alias
    @test sort(collect(buffer)) == values
    @test all(isone,
        canonicalize(portion, SciMLStructures.replace(portion, ps, ones(length(buffer))))[1])
    # make sure it is out-of-place
    @test sort(collect(buffer)) == values
    SciMLStructures.replace!(portion, ps, ones(length(buffer)))
    # make sure it is in-place
    @test all(isone, canonicalize(portion, ps)[1])
    repack(zeros(length(buffer)))
    @test all(iszero, canonicalize(portion, ps)[1])
end

setp(sys, a)(ps, 2.0) # test set_parameter!
@test getp(sys, a)(ps) == 2.0

setp(sys, e)(ps, 5ones(3)) # with an array
@test getp(sys, e)(ps) == 5ones(3)

setp(sys, f[2, 2])(ps, 42) # with a sub-index
@test getp(sys, f[2, 2])(ps) == 42

# SII bug
@test_broken setp(sys, g)(ps, ones(100)) # with non-fixed-length array
@test_broken getp(sys, g)(ps) == ones(100)

setp(sys, h)(ps, "bar") # with a non-numeric
@test getp(sys, h)(ps) == "bar"
