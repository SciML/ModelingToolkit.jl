using ModelingToolkitBase, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using Symbolics, SciMLBase
import ModelingToolkitBase as MTKB

replace_at(v, pairs...) = (w = copy(v); foreach(((i, x),) -> w[i] = x, pairs); w)

@testset "array unknowns resolve identically from array, per-element and mixed values" begin
    n = 6
    @variables u(t)[1:n]
    @parameters p q b[1:n]
    sys = complete(System([D(u) ~ -u + b], t, [u], [p, q, b]; name = :arrsys))
    bval = collect(1.0:n)
    uval = collect(10.0:(9 + n))
    duval = collect(100.0:(99 + n))

    @testset "array value" begin
        prob = DAEProblem(
            sys, [u => uval, D(u) => duval, b => bval, p => 1.0, q => 2.0], (0.0, 1.0);
            build_initializeprob = false
        )
        @test prob.u0 == uval
        @test prob.du0 == duval
    end

    @testset "per-element values" begin
        op = Pair[b => bval, p => 1.0, q => 2.0]
        for i in 1:n
            push!(op, u[i] => uval[i])
            push!(op, D(u)[i] => duval[i])
        end
        prob = DAEProblem(sys, op, (0.0, 1.0); build_initializeprob = false)
        @test prob.u0 == uval
        @test prob.du0 == duval
    end

    @testset "mixed numeric, parameter and element-referencing values" begin
        op = [
            u[1] => p, u[2] => 2.0, u[3] => u[2] + q, u[4] => 4.0, u[5] => 2 * u[1], u[6] => 6.0,
            D(u)[1] => 1.0, D(u)[2] => q, D(u)[3] => 3.0, D(u)[4] => D(u)[2] + 1, D(u)[5] => 5.0,
            D(u)[6] => p + q, b => bval, p => 3.0, q => 2.0,
        ]
        prob = DAEProblem(sys, op, (0.0, 1.0); build_initializeprob = false)
        @test prob.u0 == [3.0, 2.0, 4.0, 4.0, 6.0, 6.0]
        @test prob.du0 == [1.0, 2.0, 3.0, 3.0, 5.0, 5.0]
    end

    @testset "symbolic array value" begin
        prob = ODEProblem(
            sys, [u => [p * i for i in 1:n], p => 3.0, q => 2.0, b => bval], (0.0, 1.0);
            build_initializeprob = false
        )
        @test prob.u0 == 3.0 .* (1:n)
    end
end

@testset "varmap_to_vars resolves symbolic array values on IR and non-IR paths" begin
    n = 5
    @variables x(t)[1:n]
    @parameters p
    sys = complete(System([D(x) ~ -x], t, [x], [p]; name = :xsys))
    ir = MTKB.get_irstructure(sys)
    vars = Symbolics.unwrap.(collect(x))
    expected = [3.0 * i for i in 1:n]
    chain = Dict{Any, Any}(x[i] => x[i - 1] + 3.0 for i in 2:n)
    chain[x[1]] = p
    chain[p] = 3.0
    for options in ((; ir), (;))
        @test MTKB.varmap_to_vars(
            Dict{Any, Any}(x => [p * i for i in 1:n], p => 3.0), vars; options...
        ) == expected
        @test MTKB.varmap_to_vars(copy(chain), vars; options...) == expected
        # symbolic values that cannot be resolved are left as they are
        unresolved = MTKB.varmap_to_vars(
            Dict{Any, Any}(x => fill(p, n)), vars; options..., allow_symbolic = true
        )
        @test all(isequal(Symbolics.unwrap(p)), unresolved)
    end
end

@testset "symbolic values of a matrix unknown keep their element order" begin
    @variables m(t)[1:2, 1:3]
    @parameters p
    sys = complete(System([D(m) ~ -m], t, [m], [p]; name = :msys))
    ir = MTKB.get_irstructure(sys)
    vars = Symbolics.unwrap.(vec(collect(m)))
    value = [p * (i + 10j) for i in 1:2, j in 1:3]
    for options in ((; ir), (;))
        @test MTKB.varmap_to_vars(Dict{Any, Any}(m => value, p => 2.0), vars; options...) ==
            vec([2.0 * (i + 10j) for i in 1:2, j in 1:3])
    end
end
