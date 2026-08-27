using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using Symbolics
using OrdinaryDiffEq
using SymbolicIndexingInterface
using SciMLBase
using LinearAlgebra
using Test

# A record of scalars, and a record containing a record. Both concrete, as struct
# variables require.
struct Pair2
    x::Real
    y::Real
end
struct Nested
    x::Real
    f::Pair2
end
@symstruct Pair2
@symstruct Nested

@testset "declaration and field access" begin
    @variables s(t)::Pair2 n(t)::Nested

    @test s isa Symbolics.SymStruct{Pair2}
    @test Symbolics.symtype(Symbolics.unwrap(s)) === Pair2
    @test Symbolics.symtype(Symbolics.unwrap(s.x)) === Real

    # A field access is a projection of the one variable, not a free-standing variable,
    # so variable search reports the leaves rather than inventing new roots.
    @test issetequal(
        Symbolics.unwrap.([s.x, s.y]),
        collect(Symbolics.SymStruct{Pair2}(Symbolics.unwrap(s)))
    )

    # Nesting composes, and leaf enumeration flattens through it in declaration order.
    @test isequal(Symbolics.unwrap(n.f.x), collect(Symbolics.SymStruct{Nested}(Symbolics.unwrap(n)))[2])
    @test length(Symbolics.SymStruct{Nested}(Symbolics.unwrap(n))) == 3
end

@testset "flat record integrates" begin
    @variables s(t)::Pair2
    # s.x' = -s.x, s.y' = s.x  =>  s.x = e^-t, s.y = 1 - e^-t
    sys = mtkcompile(System([D(s.x) ~ -s.x, D(s.y) ~ s.x], t, [s], []; name = :flat))

    # The record is one unknown before compilation; its leaves are the states after.
    @test issetequal(unknowns(sys), Symbolics.unwrap.([s.x, s.y]))

    prob = ODEProblem(sys, [s.x => 1.0, s.y => 0.0], (0.0, 1.0))
    # No algebraic equations, so this must stay a plain ODE.
    @test prob.f.mass_matrix === I

    sol = solve(prob, Tsit5())
    @test SciMLBase.successful_retcode(sol)
    @test sol[s.x][end]≈exp(-1) rtol=1e-5
    @test sol[s.y][end]≈1 - exp(-1) rtol=1e-5
end

@testset "nested record integrates" begin
    @variables n(t)::Nested
    # n.x'   = n.x + t         => n.x   = 2e^t - t - 1
    # n.f.x' = n.x + n.f.x + 1 => n.f.x = e^t + 2t*e^t + t + 1
    # n.f.y' = n.f.y           => n.f.y = 3e^t
    eqs = [D(n.x) ~ n.x + t
           D(n.f.x) ~ n.x + n.f.x + 1
           D(n.f.y) ~ n.f.y]
    sys = mtkcompile(System(eqs, t, [n], []; name = :nested))

    @test issetequal(unknowns(sys), Symbolics.unwrap.([n.x, n.f.x, n.f.y]))

    prob = ODEProblem(sys, [n.x => 1.0, n.f.x => 2.0, n.f.y => 3.0], (0.0, 1.0))
    sol = solve(prob, Tsit5())
    @test SciMLBase.successful_retcode(sol)
    @test sol[n.x][end]≈2exp(1) - 2 rtol=1e-5
    @test sol[n.f.x][end]≈3exp(1) + 2 rtol=1e-5
    @test sol[n.f.y][end]≈3exp(1) rtol=1e-5
end

@testset "unknowns are inferred from field accesses" begin
    @variables n(t)::Nested
    eqs = [D(n.x) ~ -n.x, D(n.f.x) ~ n.x, D(n.f.y) ~ n.f.y]
    # No unknowns passed: they must be discovered from the equations, at leaf granularity.
    sys = mtkcompile(System(eqs, t; name = :inferred))
    @test issetequal(unknowns(sys), Symbolics.unwrap.([n.x, n.f.x, n.f.y]))
end

@testset "initial conditions" begin
    @variables n(t)::Nested
    eqs = [D(n.x) ~ n.x, D(n.f.x) ~ n.f.x, D(n.f.y) ~ n.f.y]
    sys = mtkcompile(System(eqs, t, [n], []; name = :ics))
    idx = Base.Fix1(variable_index, sys)

    # A whole-record value is split field-wise, as an array value is split element-wise.
    whole = ODEProblem(sys, [n => Nested(1, Pair2(2, 3))], (0.0, 1.0))
    @test whole.u0[idx(n.x)] == 1
    @test whole.u0[idx(n.f.x)] == 2
    @test whole.u0[idx(n.f.y)] == 3

    # Leaf-wise entries are equivalent.
    leafwise = ODEProblem(
        sys, [n.x => 1.0, n.f.x => 2.0, n.f.y => 3.0], (0.0, 1.0))
    @test leafwise.u0 == whole.u0

    # An explicit leaf entry overrides the record it came from.
    override = ODEProblem(
        sys, [n => Nested(1, Pair2(2, 3)), n.f.y => 99.0], (0.0, 1.0))
    @test override.u0[idx(n.f.y)] == 99
    @test override.u0[idx(n.x)] == 1
end

@testset "symbolic indexing" begin
    @variables n(t)::Nested
    sys = mtkcompile(System([D(n.x) ~ -n.x, D(n.f.x) ~ n.x, D(n.f.y) ~ n.f.y], t, [n], []; name = :sii))
    prob = ODEProblem(sys, [n.x => 1.0, n.f.x => 2.0, n.f.y => 3.0], (0.0, 1.0))

    for (leaf, val) in ((n.x, 1.0), (n.f.x, 2.0), (n.f.y, 3.0))
        @test is_variable(sys, leaf)
        @test variable_index(sys, leaf) isa Int
        @test prob[leaf] == val
    end
    # Distinct leaves must occupy distinct slots.
    @test allunique(variable_index.((sys,), [n.x, n.f.x, n.f.y]))
    # Derivatives of a leaf resolve to the state introduced for them.
    @test is_variable(sys, n.f.y) && !is_variable(sys, D(n.f.y))
end

@testset "namespacing" begin
    @variables s(t)::Pair2
    @named inner = System([D(s.x) ~ -s.x, D(s.y) ~ s.x], t, [s], [])
    @named outer = System(Equation[], t, [], []; systems = [inner])
    # The name belongs to the record; the field access is structural.
    @test isequal(Symbolics.unwrap(outer.inner.s.x), Symbolics.unwrap(getproperty(outer.inner, :s).x))
    @test occursin("inner", string(Symbolics.unwrap(outer.inner.s.x)))
end

# A record holding an array of records, so an access path interleaves fields and indices
# (`h.x[1].y`) rather than grouping them.
struct Leaf
    y::Real
    z::Real
end
struct Holder
    x::Vector{Leaf}
end
@symstruct Leaf
@symstruct Holder begin
    shape(:x) = [1:2]
end

@testset "access paths interleave fields and indices" begin
    @variables h(t)::Holder n(t)::Nested u(t)[1:2] w(t)
    SI = Symbolics.SymbolicUtils.StableIndex

    # The whole chain resolves to the record, in either order.
    @test isequal(
        ModelingToolkitBase.split_field_access(Symbolics.unwrap(n.f.x)),
        (Symbolics.unwrap(n), true))
    @test isequal(
        ModelingToolkitBase.split_field_access(Symbolics.unwrap(h.x[1].y)),
        (Symbolics.unwrap(h), true))

    # Paths are recorded outermost last, so the index sits where it was written.
    @test isequal(ModelingToolkitBase.get_struct_access(Symbolics.unwrap(n.f.x)), [:f, :x])
    @test isequal(
        ModelingToolkitBase.get_struct_access(Symbolics.unwrap(h.x[1].y)),
        [:x, SI([1]), :y])

    # Commutes with operators, as `split_indexed_var` does.
    @test isequal(
        ModelingToolkitBase.split_field_access(Symbolics.unwrap(D(h.x[2].z)))[1],
        D(Symbolics.unwrap(h)))
    @test isequal(
        ModelingToolkitBase.get_struct_access(Symbolics.unwrap(D(h.x[2].z))),
        [:x, SI([2]), :z])

    # Plain arrays still resolve to the array; plain scalars are not access paths.
    @test isequal(
        ModelingToolkitBase.split_field_access(Symbolics.unwrap(u[1])),
        (Symbolics.unwrap(u), true))
    @test isequal(
        ModelingToolkitBase.split_field_access(Symbolics.unwrap(w)),
        (Symbolics.unwrap(w), false))
end

@testset "record of records-array integrates" begin
    @variables h(t)::Holder
    eqs = [D(h.x[1].y) ~ -h.x[1].y, D(h.x[1].z) ~ h.x[1].y,
           D(h.x[2].y) ~ -h.x[2].y, D(h.x[2].z) ~ h.x[2].y]
    sys = mtkcompile(System(eqs, t, [h], []; name = :holder))
    @test issetequal(
        unknowns(sys),
        Symbolics.unwrap.([h.x[1].y, h.x[1].z, h.x[2].y, h.x[2].z])
    )

    prob = ODEProblem(
        sys,
        [h.x[1].y => 1.0, h.x[1].z => 0.0, h.x[2].y => 2.0, h.x[2].z => 0.0],
        (0.0, 1.0)
    )
    sol = solve(prob, Tsit5())
    @test SciMLBase.successful_retcode(sol)
    @test sol[h.x[1].y][end]≈exp(-1) rtol=1e-5
    @test sol[h.x[2].y][end]≈2exp(-1) rtol=1e-5
end

@testset "array variables are unaffected" begin
    @variables u(t)[1:2]
    @parameters q[1:2]
    sys = mtkcompile(System([D(u[1]) ~ q[1] * u[2], D(u[2]) ~ -q[2] * u[1]], t, [u], [q]; name = :arr))
    prob = ODEProblem(sys, [u => [1.0, 0.0], q => [2.0, 3.0]], (0.0, 1.0))
    @test prob[u] == [1.0, 0.0]
    @test prob[u[1]] == 1.0
    @test prob.ps[q] == [2.0, 3.0]
    @test SciMLBase.successful_retcode(solve(prob, Tsit5()))
end
