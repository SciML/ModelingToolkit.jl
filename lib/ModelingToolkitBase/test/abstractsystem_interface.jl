using ModelingToolkitBase
using SymbolicIndexingInterface: SymbolicIndexingInterface as SII
using Test
MT = ModelingToolkitBase

# Every call below is part of the documented `AbstractSystem` interface (see
# `docs/src/API/abstract_system_interface.md`). Nothing here reads a field of a system
# directly or calls an entry point that is specific to the `System` implementation.

@independent_variables t
@variables x(t) y(t)
@parameters p q
D = Differential(t)

"A subtype carrying only the two required fields."
struct MinimalSystem <: MT.AbstractSystem
    name::Symbol
    systems::Vector
end

"A subtype carrying the optional fields backing the equation/unknown/parameter accessors."
struct ComponentSystem <: MT.AbstractSystem
    eqs::Vector{Equation}
    unknowns::Vector
    ps::Vector
    iv::Any
    name::Symbol
    systems::Vector
end

"A subtype which supplies its independent variable through the scalar extension point."
struct RenamedIVSystem <: MT.AbstractSystem
    time::Any
    name::Symbol
    systems::Vector
end
MT.independent_variable(sys::RenamedIVSystem) = getfield(sys, :time)

leaf = ComponentSystem([D(x) ~ p * x], [x], [p], t, :leaf, [])
branch = ComponentSystem([D(y) ~ q * y], [y], [q], t, :branch, [leaf])
root = ComponentSystem(Equation[], [], [], t, :root, [branch])

@testset "Required surface" begin
    sys = MinimalSystem(:sys, [])
    @test nameof(sys) === :sys
    @test isempty(MT.get_systems(sys))

    parent = MinimalSystem(:parent, [sys])
    @test map(nameof, MT.get_systems(parent)) == [:sys]
end

@testset "Accessors defaulted for absent optional fields" begin
    sys = MinimalSystem(:sys, [])
    @test MT.independent_variable(sys) === nothing
    @test isempty(independent_variables(sys))
    @test MT.description(sys) == ""
    @test MT.iscomplete(sys) == false
    @test MT.does_namespacing(sys) == true
    @test isempty(MT.assertions(sys))
    @test isempty(MT.continuous_events(sys))
    @test isempty(MT.discrete_events(sys))
    @test isempty(SII.independent_variable_symbols(sys))
    @test SII.constant_structure(sys)
end

@testset "`has_x` reports optional field availability" begin
    minimal = MinimalSystem(:sys, [])

    @test MT.has_name(minimal)
    @test MT.has_name(leaf)

    for hasfn in (MT.has_eqs, MT.has_unknowns, MT.has_ps, MT.has_iv)
        @test hasfn(minimal) == false
        @test hasfn(leaf) == true
    end

    # Neither subtype stores observed equations, so neither may be asked for them.
    @test MT.has_observed(minimal) == false
    @test MT.has_observed(leaf) == false

    @test MT.get_eqs(leaf) == [D(x) ~ p * x]
    @test isequal(MT.get_unknowns(leaf), [x])
    @test isequal(MT.get_ps(leaf), [p])
    @test isequal(MT.get_iv(leaf), t)

    # `get_x` on absent storage throws, so generic code has to guard it with `has_x`.
    @test_throws ErrorException MT.get_eqs(minimal)
    @test_throws ErrorException MT.get_unknowns(minimal)
    @test_throws ErrorException MT.get_ps(minimal)

    # The event accessors are the documented exceptions: they default to empty.
    @test MT.has_continuous_events(minimal) == false
    @test MT.has_discrete_events(minimal) == false
    @test isempty(MT.get_continuous_events(minimal))
    @test isempty(MT.get_discrete_events(minimal))
end

@testset "Scalar independent variable extension point" begin
    sys = RenamedIVSystem(t, :sys, [])
    @test isequal(MT.independent_variable(sys), t)
    @test isequal(independent_variables(sys), [t])
    @test isequal(SII.independent_variable_symbols(sys), [t])
end

@testset "Hierarchical accessors namespace subsystem contents" begin
    branch_y = MT.renamespace(branch, y)
    branch_q = MT.renamespace(branch, q)
    branch_leaf_x = MT.renamespace(branch, MT.renamespace(leaf, x))
    branch_leaf_p = MT.renamespace(branch, MT.renamespace(leaf, p))

    @test isequal(unknowns(root), [branch_y, branch_leaf_x])
    @test isequal(parameters(root), [branch_q, branch_leaf_p])
    @test equations(root) == [
        D(branch_y) ~ branch_q * branch_y,
        D(branch_leaf_x) ~ branch_leaf_p * branch_leaf_x,
    ]

    @test isequal(unknowns(leaf), [x])
    @test isequal(parameters(leaf), [p])
    @test equations(leaf) == [D(x) ~ p * x]

    # Locally stored entries come first, then the subsystems in `get_systems` order.
    leaf_x = MT.renamespace(leaf, x)
    leaf_p = MT.renamespace(leaf, p)
    @test isequal(unknowns(branch), [y, leaf_x])
    @test isequal(parameters(branch), [q, leaf_p])
    @test equations(branch) == [D(y) ~ q * y, D(leaf_x) ~ leaf_p * leaf_x]
end

@testset "Toplevel accessors ignore subsystems" begin
    @test isempty(MT.equations_toplevel(root))
    @test isempty(MT.unknowns_toplevel(root))
    @test isempty(MT.parameters_toplevel(root))
    @test isempty(MT.continuous_events_toplevel(root))
    @test isempty(MT.discrete_events_toplevel(root))

    @test MT.equations_toplevel(branch) == [D(y) ~ q * y]
    @test isequal(MT.unknowns_toplevel(branch), [y])
    @test isequal(MT.parameters_toplevel(branch), [q])
end

@testset "Equation classification" begin
    algebraic = ComponentSystem([0 ~ p - x], [x], [p], t, :algebraic, [])

    @test has_diff_equations(leaf)
    @test !has_alg_equations(leaf)
    @test diff_equations(leaf) == [D(x) ~ p * x]
    @test isempty(alg_equations(leaf))

    @test has_alg_equations(algebraic)
    @test !has_diff_equations(algebraic)
    @test alg_equations(algebraic) == [0 ~ p - x]
    @test isempty(diff_equations(algebraic))
end

@testset "Namespacing helpers" begin
    leaf_x = MT.renamespace(leaf, x)
    leaf_p = MT.renamespace(leaf, p)
    @test MT.namespace_equations(leaf) == [D(leaf_x) ~ leaf_p * leaf_x]
end

@testset "SymbolicIndexingInterface queries backed by optional fields" begin
    @test isequal(SII.variable_symbols(leaf), [x])
    @test isequal(SII.parameter_symbols(leaf), [p])
    @test SII.is_variable(leaf, x)
    @test !SII.is_variable(leaf, p)
    @test SII.variable_index(leaf, x) == 1
    @test SII.is_parameter(leaf, p)
    @test !SII.is_parameter(leaf, x)
    @test SII.parameter_index(leaf, p) == 1
end
