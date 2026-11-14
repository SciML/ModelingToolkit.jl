using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D, MTKParameters
using ModelingToolkitBase: SciMLBase
using ModelingToolkitBase.StaticArrays: SArray, SVector
using Enzyme
using Enzyme: EnzymeCore
using ForwardDiff
using StableRNGs
using Test

# `make_zero!`/`remake_zero!` on an `MTKParameters` shadow must leave it exactly as
# Enzyme's generic methods do (what is zeroed, what throws, the `seen` bookkeeping),
# never write to the primal (parts of a `make_zero` shadow are the primal's own objects:
# integer buffers, inactive nonnumeric parameters), and skip the generic `IdSet` walk
# for the float buffers of plain models.

const LEAF = Union{Number, Symbol, String, Module, Type, Function}

# every leaf value reachable from `x`, in a fixed order, with undefined slots and the
# aliasing of mutable objects (as the index of their first visit)
collect_leaves(x) = collect_leaves!(Any[], IdDict{Any, Int}(), x)
collect_leaves!(out, ids, x::LEAF) = (push!(out, x); out)
function collect_leaves!(out, ids, x::ForwardDiff.Dual)
    push!(out, ForwardDiff.value(x))
    append!(out, ForwardDiff.partials(x))
    return out
end
function collect_leaves!(out, ids, x)
    if ismutable(x)
        if haskey(ids, x)
            push!(out, (:alias, ids[x]))
            return out
        end
        ids[x] = length(ids) + 1
    end
    if x isa Union{Array, Tuple}
        for i in eachindex(x)
            if !(x isa Array) || isassigned(x, i)
                collect_leaves!(out, ids, x[i])
            else
                push!(out, :undef)
            end
        end
    else
        for i in 1:nfields(x)
            isdefined(x, i) ? collect_leaves!(out, ids, getfield(x, i)) : push!(out, :undef)
        end
    end
    return out
end
same_leaves(a, b) = isequal(collect_leaves(a), collect_leaves(b))

# mutable objects reachable from `x`
reachable(x) = reachable!(Base.IdSet{Any}(), x)
reachable!(set, ::LEAF) = set
function reachable!(set, x)
    if ismutable(x)
        x in set && return set
        push!(set, x)
    end
    if x isa Array
        for i in eachindex(x)
            isassigned(x, i) && reachable!(set, x[i])
        end
    else
        for i in 1:nfields(x)
            isdefined(x, i) && reachable!(set, getfield(x, i))
        end
    end
    return set
end

# isbits types made only of one float type (`SVector`s, `Dual`s, complex numbers) can be
# dirtied through a view reinterpreted as that float type
function single_float(T::Type)
    T <: Union{Float16, Float32, Float64} && return T
    (isbitstype(T) && fieldcount(T) > 0) || return nothing
    Fs = unique(map(single_float, fieldtypes(T)))
    return length(Fs) == 1 ? only(Fs) : nothing
end

# overwrite every float the shadow owns (not the parts it shares with the primal); the
# shadow is a fresh `make_zero`, so this stands in for an accumulated gradient
dirty!(rng, x, primal) = x in primal ? x : _dirty!(rng, x, primal)
_dirty!(rng, x::LEAF, primal) = x
function _dirty!(rng, x::Array{T}, primal) where {T}
    if T === Union{}
    elseif T <: Union{Float16, Float32, Float64}
        x .= 1 .+ rand(rng, T, size(x))
    elseif single_float(T) !== nothing
        F = single_float(T)
        v = reinterpret(F, x)
        v .= 1 .+ rand(rng, F, length(v))
    elseif T <: Union{AbstractFloat, Complex}
        for i in eachindex(x)
            x[i] = x[i] isa Complex ? complex(1 + rand(rng), 1 + rand(rng)) : 1 + rand(rng)
        end
    else
        for i in eachindex(x)
            isassigned(x, i) && dirty!(rng, x[i], primal)
        end
    end
    return x
end
function _dirty!(rng, x::Tuple, primal)
    foreach(xi -> dirty!(rng, xi, primal), x)
    return x
end
function _dirty!(rng, x, primal)
    for i in 1:nfields(x)
        isdefined(x, i) || continue
        xi = getfield(x, i)
        if ismutable(x) && xi isa AbstractFloat && !isconst(typeof(x), i)
            setfield!(x, i, one(xi) + rand(rng, typeof(xi)))
        else
            dirty!(rng, xi, primal)
        end
    end
    return x
end

# callable nonnumeric parameter holding differentiable data
struct ActiveScale
    w::Vector{Float64}
end
(s::ActiveScale)(x) = s.w[1] * x
# callable nonnumeric parameter declared constant: `make_zero` shares it with the primal
struct InactiveTable
    v::Vector{Float64}
end
(s::InactiveTable)(x) = s.v[1] * x
EnzymeCore.EnzymeRules.inactive_type(::Type{<:InactiveTable}) = true

# a mutable buffer type with its own `make_zero!` and no `remake_zero!` (as GPU arrays):
# the generic walk zeroes such a field through `make_zero!`, also from `remake_zero!`
const OWN_MAKE_ZERO_CALLS = Ref(0)
mutable struct OwnZeroVector <: AbstractVector{Float64}
    data::Vector{Float64}
end
Base.size(x::OwnZeroVector) = size(x.data)
Base.getindex(x::OwnZeroVector, i::Int) = x.data[i]
Base.setindex!(x::OwnZeroVector, v, i::Int) = (x.data[i] = v)
function EnzymeCore.make_zero!(x::OwnZeroVector, seen)
    OWN_MAKE_ZERO_CALLS[] += 1
    fill!(x.data, 0.0)
    return nothing
end
EnzymeCore.make_zero!(x::OwnZeroVector) = EnzymeCore.make_zero!(x, nothing)

# the same, immutable: the generic walk still zeroes it through its `make_zero!`
struct ImmutableOwnZeroVector <: AbstractVector{Float64}
    data::Vector{Float64}
end
Base.size(x::ImmutableOwnZeroVector) = size(x.data)
Base.getindex(x::ImmutableOwnZeroVector, i::Int) = x.data[i]
function EnzymeCore.make_zero!(x::ImmutableOwnZeroVector, seen)
    OWN_MAKE_ZERO_CALLS[] += 1
    fill!(x.data, 0.0)
    return nothing
end

# an array type without elements that still holds differentiable data
struct EmptyWithData <: AbstractVector{Union{}}
    data::Vector{Float64}
end
Base.size(::EmptyWithData) = (0,)

# types that contain themselves, inactive and active
struct RecursiveNode
    children::Vector{RecursiveNode}
    w::Vector{Float64}
end
EnzymeCore.EnzymeRules.inactive_type(::Type{RecursiveNode}) = true
struct RecursiveVector <: AbstractVector{Float64}
    next::Vector{RecursiveVector}
    data::Vector{Float64}
end
Base.size(x::RecursiveVector) = size(x.data)
Base.getindex(x::RecursiveVector, i::Int) = x.data[i]

# a mutable object with its own `make_zero!` that counts the objects it zeroes, for an
# object that is both a nonnumeric parameter and a cache
const TRACKED_ZEROED = Ref(0)
mutable struct TrackedZero
    w::Vector{Float64}
end
function EnzymeCore.make_zero!(x::TrackedZero, seen)
    x in seen && return nothing
    push!(seen, x)
    TRACKED_ZEROED[] += 1
    fill!(x.w, 0.0)
    return nothing
end

# types declared inactive (unique to this file): buffer tuples, and a whole
# `MTKParameters`, whose `make_zero` shadows are the primal's own objects
struct InactiveTag end
struct InactiveTupleTag end
const InactiveTuple = Tuple{Vector{InactiveTupleTag}, Vector{Float64}}
EnzymeCore.EnzymeRules.inactive_type(::Type{InactiveTuple}) = true
const P_INACTIVE = MTKParameters([7.0], Float64[], (), (), ([InactiveTag()],), ())
EnzymeCore.EnzymeRules.inactive_type(::Type{typeof(P_INACTIVE)}) = true
EnzymeCore.EnzymeRules.inactive_type(::Type{Tuple{Vector{Float32}, Vector{ComplexF32}}}) = true

# an immutable buffer whose float data Enzyme sees only through its `is_mutable_array`
# hook and the buffer's own `make_zero`/`make_zero!`
struct EncodedFloats <: AbstractVector{Float64}
    bits::Vector{UInt64}
end
Base.size(x::EncodedFloats) = size(x.bits)
Base.getindex(x::EncodedFloats, i::Int) = reinterpret(Float64, x.bits[i])
Enzyme.Compiler.is_mutable_array(::Type{EncodedFloats}) = true
function EnzymeCore.make_zero(
        ::Type{EncodedFloats}, seen::IdDict, x::EncodedFloats, ::Val{C} = Val(false)
    ) where {C}
    haskey(seen, x) && return seen[x]
    return seen[x] = EncodedFloats(zeros(UInt64, length(x)))
end
EnzymeCore.make_zero!(x::EncodedFloats, seen) = (fill!(x.bits, 0); nothing)

# a nonnumeric parameter sharing the tunable buffer, whose own `make_zero!` checks
# whether the generic walk has already visited that buffer
const SAW_SHARED_TUNABLE = Ref(0)
mutable struct SharesTunable
    data::Vector{Float64}
end
function EnzymeCore.make_zero!(x::SharesTunable, seen)
    x in seen && return nothing
    push!(seen, x)
    x.data in seen && (SAW_SHARED_TUNABLE[] += 1)
    fill!(x.data, 0.0)
    return nothing
end

# array-valued parameter nested beyond one level, sharing its subtrees
nested_tree(x, n) = n == 0 ? x : (t = nested_tree(x, n - 1); [[t], [t]])

# a nonnumeric tag declared inactive
struct ConstTag end
EnzymeCore.EnzymeRules.inactive_type(::Type{ConstTag}) = true
# a float type declared inactive, which Enzyme's activity analysis still treats as active
# (floats are classified before `inactive_type` is consulted)
struct InactiveDeclaredFloat <: AbstractFloat
    x::Float64
end
Base.zero(::Type{InactiveDeclaredFloat}) = InactiveDeclaredFloat(0.0)
EnzymeCore.EnzymeRules.inactive_type(::Type{InactiveDeclaredFloat}) = true
collect_leaves!(out, ids, x::InactiveDeclaredFloat) = (push!(out, x.x); out)

allocations(f, x) = (f(x); @allocated f(x))

# the generic method on the `MTKParameters` itself (inside the despecialization wrapper,
# which only forwards to it)
function generic!(zero!, s, seen)
    target = s isa SciMLBase.DespecializedParameters ? getfield(s, :params) : s
    return invoke(zero!, Tuple{Any, Any}, target, seen)
end

# whether the one-argument methods skip the generic walk for `p`
fast_zero(p) = ModelingToolkitBase.enzyme_fast_zero(typeof(p))
fast_zero(p::SciMLBase.DespecializedParameters) = false

function check_zeroing(p; label, fast = false)
    return @testset "$label" begin
        @test fast_zero(p) == fast
        rng = StableRNG(1)
        primal = deepcopy(p)
        shadow = Enzyme.make_zero(p)
        dirty!(rng, shadow, reachable(p))
        @test same_leaves(p, primal)  # dirtying touched only shadow-owned data
        for zero! in (Enzyme.make_zero!, Enzyme.remake_zero!)
            expected = deepcopy(shadow)
            calls = OWN_MAKE_ZERO_CALLS[]
            generic_error = try
                generic!(zero!, expected, Base.IdSet{Any}())
                nothing
            catch err
                err isa ArgumentError || rethrow()
                sprint(showerror, err)
            end
            generic_calls = OWN_MAKE_ZERO_CALLS[] - calls
            for args in ((), (Base.IdSet{Any}(),))
                s = deepcopy(shadow)
                calls = OWN_MAKE_ZERO_CALLS[]
                our_error = try
                    zero!(s, args...)
                    nothing
                catch err
                    err isa ArgumentError || rethrow()
                    sprint(showerror, err)
                end
                @test our_error == generic_error
                # also the state left behind by an error
                @test same_leaves(s, expected)
                @test OWN_MAKE_ZERO_CALLS[] - calls == generic_calls
            end
        end
        # in place, on the shadow that shares objects with the primal
        try
            Enzyme.make_zero!(shadow)
        catch err
            err isa ArgumentError || rethrow()
        end
        Enzyme.remake_zero!(shadow)
        @test same_leaves(p, primal)
    end
end

# `seen` bookkeeping as in the generic methods: visited objects are registered, and
# objects already in the set are skipped
function check_seen(p; label)
    return @testset "$label: seen set" begin
        rng = StableRNG(2)
        shadow = dirty!(rng, Enzyme.make_zero(p), reachable(p))
        for zero! in (Enzyme.make_zero!, Enzyme.remake_zero!)
            for preseen in (s -> (), s -> (s,), s -> (s.tunable,))
                s, e = deepcopy(shadow), deepcopy(shadow)
                seen_s = Base.IdSet{Any}(preseen(s))
                seen_e = Base.IdSet{Any}(preseen(e))
                zero!(s, seen_s)
                invoke(zero!, Tuple{Any, Any}, e, seen_e)
                @test same_leaves(s, e)
                @test (s in seen_s) == (e in seen_e)
                @test (s.tunable in seen_s) == (e.tunable in seen_e)
            end
        end
    end
end

# `MTKParameters` built directly from buffers, for buffer types the constructor from a
# system does not produce
mtkp(tunable; constant = ()) = MTKParameters(tunable, Float64[], (), constant, (), ())

@testset "make_zero!/remake_zero! of an MTKParameters shadow" begin
    # a plain model: float tunables, initial values, a non-tunable float and a discrete
    # float parameter, as `MTKParameters(sys, op)` builds them
    @variables y(t)
    @discretes q(t)
    @parameters α β [tunable = false]
    @named plain = System(
        [D(y) ~ -α * y + β * q], t;
        continuous_events = [
            ModelingToolkitBase.SymbolicContinuousCallback(
                [y ~ 0.5] => [q ~ 2.0], discrete_parameters = q
            ),
        ]
    )
    plain = mtkcompile(plain)
    for T in (Float64, Float32)
        op = [y => T(1), α => T(1), β => T(2), q => T(0.5)]
        p_plain = SciMLBase.unwrap_parameters(
            ODEProblem{true, SciMLBase.FullSpecialize}(plain, op, (T(0), T(1))).p
        )
        @test p_plain.tunable isa Vector{T} && p_plain.initials isa Vector{T}
        @test p_plain.discrete isa Tuple{ModelingToolkitBase.BlockedArray{T}}
        @test p_plain.constant isa Tuple{Vector{T}}
        @test all(!isempty, (p_plain.tunable, p_plain.initials, p_plain.discrete[1], p_plain.constant[1]))
        @test fast_zero(p_plain)
        check_zeroing(p_plain; label = "plain $T model", fast = true)
        check_seen(p_plain; label = "plain $T model")
        # zeroed without the generic `IdSet` walk
        @test allocations(Enzyme.remake_zero!, Enzyme.make_zero(p_plain)) == 0
        @test allocations(Enzyme.make_zero!, Enzyme.make_zero(p_plain)) == 0
    end

    # every other layout is left to the generic walk; these compare it against itself
    # through the one-argument methods
    @discretes c(t)
    @parameters a b d::Integer e[1:3] f[1:3, 1:3]::Int g::Vector{AbstractFloat} h::String
    @named sys = System(
        Equation[], t, [], [a, b, c, d, e, f, g, h];
        continuous_events = [
            ModelingToolkitBase.SymbolicContinuousCallback(
                [a ~ 0] => [c ~ 0], discrete_parameters = c
            ),
        ], bindings = [b => 2a],
        initial_conditions = Dict(a => 0.0)
    )
    sys = complete(sys)
    ivs = Dict(
        a => 1.5, c => 3a, d => 4, e => [5.0, 6.0, 7.0],
        f => ones(Int, 3, 3), g => [0.1, 0.2, 0.3], h => "foo"
    )
    ps = MTKParameters(sys, ivs)
    @test ps.discrete isa Tuple{ModelingToolkitBase.BlockedArray{Float64}}
    @test !fast_zero(ps)
    check_zeroing(ps; label = "float, integer, array and string buffers")
    check_seen(ps; label = "float, integer, array and string buffers")
    check_zeroing(MTKParameters(sys, ivs; floatT = Float32); label = "Float32 buffers")
    ps_static = MTKParameters(sys, ivs; p_constructor = x -> SArray{Tuple{size(x)...}}(x))
    @test ps_static.tunable isa SVector && ps_static.initials isa SVector
    check_zeroing(ps_static; label = "static buffers")
    check_zeroing(SciMLBase.DespecializedParameters(ps); label = "despecialized wrapper")

    @variables x(t)
    @parameters k (scale::ActiveScale)(..) (table::InactiveTable)(..)
    @named sys2 = System([D(x) ~ -k * scale(x) + table(x)], t)
    sys2 = mtkcompile(sys2)
    prob = ODEProblem(
        sys2,
        [x => 1.0, k => 2.0, scale => ActiveScale([0.5]), table => InactiveTable([0.1])],
        (0.0, 1.0)
    )
    p2 = SciMLBase.unwrap_parameters(prob.p)
    @test Enzyme.make_zero(p2).nonnumeric !== p2.nonnumeric
    @test !fast_zero(p2)
    check_zeroing(p2; label = "active and inactive nonnumeric parameters")

    @named sys3 = System([D(x) ~ 2x + t], t)
    sys3, _ = ModelingToolkitBase.add_diffcache(sys3, 8)
    sys3 = complete(sys3)
    p3 = SciMLBase.unwrap_parameters(ODEProblem(sys3, [x => 1.0], (0.0, 1.0)).p)
    @test !isempty(p3.caches)
    @test !fast_zero(p3)
    check_zeroing(p3; label = "DiffCache caches")

    check_zeroing(mtkp(ComplexF64[1 + 2im, 3]); label = "complex buffer", fast = true)
    check_zeroing(mtkp(Union{Float64, ComplexF64}[1.0, 2.0im]); label = "union eltype")
    check_zeroing(
        mtkp(Union{}[]; constant = (Vector{Union{}}[Union{}[]],)); label = "bottom eltype"
    )
    check_zeroing(
        mtkp(ForwardDiff.Dual{Nothing}.([1.0, 2.0], [3.0, 4.0])); label = "Dual buffer"
    )
    check_zeroing(
        mtkp(ForwardDiff.Dual{Nothing}.(Float32[1, 2], Float32[3, 4]));
        label = "Float32 Dual buffer"
    )
    check_zeroing(
        mtkp([RecursiveNode(RecursiveNode[], [7.0])]); label = "inactive recursive type"
    )
    check_zeroing(
        mtkp(RecursiveVector(RecursiveVector[], [1.0])); label = "active recursive type"
    )
    V = Vector{Float64}
    check_zeroing(
        MTKParameters{V, V, Tuple{Vararg{V}}, Tuple{}, Tuple{}, Tuple{}}(
            [1.0], Float64[], ([2.0],), (), (), ()
        );
        label = "Vararg buffer tuple type"
    )
    for n in (3, 40)
        check_zeroing(
            mtkp([1.0]; constant = (fill([1.0, 2.0], n),));
            label = "inner array shared $n times"
        )
    end
    tracked = TrackedZero([1.0])
    p_tracked = MTKParameters([1.0], Float64[], (), (), ([tracked],), ([tracked],))
    for zero! in (Enzyme.make_zero!, Enzyme.remake_zero!)
        s, e = deepcopy(Enzyme.make_zero(p_tracked)), deepcopy(Enzyme.make_zero(p_tracked))
        TRACKED_ZEROED[] = 0
        zero!(s)
        ours = TRACKED_ZEROED[]
        TRACKED_ZEROED[] = 0
        invoke(zero!, Tuple{Any, Any}, e, Base.IdSet{Any}())
        @test ours == TRACKED_ZEROED[] == 1
        @test same_leaves(s, e)
    end
    check_zeroing(
        mtkp([1.0]; constant = ([SVector(1.0, 2.0), SVector(3.0, 4.0)],));
        label = "array of static vectors"
    )
    check_zeroing(mtkp(OwnZeroVector([1.0, 2.0])); label = "buffer with its own make_zero!")
    check_seen(mtkp(OwnZeroVector([1.0, 2.0])); label = "buffer with its own make_zero!")
    check_zeroing(
        mtkp(ImmutableOwnZeroVector([1.0, 2.0]));
        label = "immutable buffer with its own make_zero!"
    )
    check_zeroing(mtkp(EmptyWithData([2.0])); label = "eltype-less array holding data")
    z0 = SVector{0, Float64}()
    check_zeroing(
        MTKParameters(Complex{Union{}}[], z0, (), (), (), ()); label = "Complex{Union{}} eltype"
    )
    check_zeroing(
        mtkp([1.0]; constant = ([2.0], SVector(0.0))); label = "mixed constant buffers"
    )
    check_zeroing(
        mtkp([1.0]; constant = ([InactiveTupleTag()], [7.0])); label = "inactive buffer tuple"
    )
    s_inactive = Enzyme.make_zero(P_INACTIVE)
    @test s_inactive.tunable === P_INACTIVE.tunable
    Enzyme.remake_zero!(s_inactive)
    Enzyme.make_zero!(s_inactive)
    @test P_INACTIVE.tunable == [7.0]
    # a buffer tuple of float vectors declared inactive: the generic walk leaves it alone
    p_inactive_floats = mtkp([1.0]; constant = (Float32[7], ComplexF32[7]))
    s_floats = Enzyme.make_zero(p_inactive_floats)
    s_floats.constant[1][1] = 7
    invoke(Enzyme.remake_zero!, Tuple{Any, Any}, s_floats, Base.IdSet{Any}())
    @test s_floats.constant[1][1] == 7
    @test !fast_zero(p_inactive_floats)
    check_zeroing(p_inactive_floats; label = "inactive float buffer tuple")

    # buffers that share arrays: aliased top-level float buffers are zeroed twice, which
    # leaves the same zeros as the generic walk's single visit
    shared_x = [1.0, 2.0]
    p_aliased = MTKParameters(shared_x, shared_x, (), (shared_x,), (), ())
    @test fast_zero(p_aliased)
    check_zeroing(p_aliased; label = "aliased top-level buffers", fast = true)
    @test allocations(Enzyme.remake_zero!, Enzyme.make_zero(p_aliased)) == 0
    check_zeroing(
        MTKParameters(shared_x, Float64[], (), ([shared_x, shared_x],), (), ());
        label = "array-valued parameter sharing the tunable buffer"
    )
    # a buffer type outside the standard layout, active only through Enzyme hooks
    p_encoded = mtkp(EncodedFloats(collect(reinterpret(UInt64, [3.0]))))
    @test !fast_zero(p_encoded)
    for zero! in (Enzyme.make_zero!, Enzyme.remake_zero!)
        s = Enzyme.make_zero(p_encoded)
        s.tunable.bits[1] = reinterpret(UInt64, 7.0)
        e = deepcopy(s)
        zero!(s)
        invoke(zero!, Tuple{Any, Any}, e, Base.IdSet{Any}())
        @test s.tunable[1] == e.tunable[1] == 0.0
    end
    # a nonnumeric parameter that shares a numeric buffer sees the generic walk's set
    shared = [1.0, 2.0]
    p_shares = MTKParameters(shared, Float64[], (), (), ([SharesTunable(shared)],), ())
    for zero! in (Enzyme.make_zero!, Enzyme.remake_zero!)
        s = Enzyme.make_zero(p_shares)
        @test s.nonnumeric[1][1].data === s.tunable
        SAW_SHARED_TUNABLE[] = 0
        zero!(s)
        ours = SAW_SHARED_TUNABLE[]
        SAW_SHARED_TUNABLE[] = 0
        invoke(zero!, Tuple{Any, Any}, Enzyme.make_zero(p_shares), Base.IdSet{Any}())
        @test ours == SAW_SHARED_TUNABLE[] == 1
    end
    # an abstractly typed buffer field
    check_zeroing(
        MTKParameters{V, V, Any, Tuple{}, Tuple{}, Tuple{}}([1.0], Float64[], ([2.0],), (), (), ());
        label = "abstract buffer field type"
    )
    # array-valued parameters, also nested beyond one level with shared subtrees
    p_nested = mtkp([1.0]; constant = (nested_tree([1.0, 2.0], 4),))
    @test !fast_zero(p_nested)
    check_zeroing(p_nested; label = "nested shared array-valued parameter")
    # nested tunable and initials buffers
    for p in (
            MTKParameters([[7.0]], Float64[], (), (), (), ()),
            MTKParameters(Float64[], [[7.0]], (), (), (), ()),
        )
        @test !fast_zero(p)
        check_zeroing(p; label = "nested tunable or initials buffer")
    end
    # an array-valued parameter whose element type has no dimension
    p_dimless = mtkp([7.0]; constant = (Vector{Array{Float64, Any}}(undef, 1),))
    @test !fast_zero(p_dimless)
    check_zeroing(p_dimless; label = "inner array type without a dimension")
    # the despecialization wrapper reaches `MTKParameters` through the generic walk with a
    # set, which these methods leave to Enzyme (values are checked above)
    @test parentmodule(
        which(Enzyme.remake_zero!, Tuple{SciMLBase.DespecializedParameters, Base.IdSet{Any}})
    ) !== ModelingToolkitBase
    p_declared_float = MTKParameters([1.0], Float64[], (), (), ([InactiveDeclaredFloat(7.0)],), ())
    @test !fast_zero(p_declared_float)
    check_zeroing(p_declared_float; label = "nonnumeric float type declared inactive")
    p_const_nonnumeric = MTKParameters(
        [1.0, 2.0], [0.5], (), ([4.0, 5.0],), ([ConstTag()], [InactiveTable([0.1])]), ()
    )
    @test !fast_zero(p_const_nonnumeric)
    check_zeroing(p_const_nonnumeric; label = "nonnumeric declared inactive")
    # with a set, these methods are not involved: the generic walk runs unchanged
    for zero! in (Enzyme.make_zero!, Enzyme.remake_zero!)
        @test parentmodule(which(zero!, Tuple{typeof(ps), Base.IdSet{Any}})) !==
            ModelingToolkitBase
    end

    # the fast path takes exactly the layout above: float vectors, blocked float vectors
    # as `MTKParameters(sys, op)` builds them, the empty static buffer, empty tuples
    blocked(x) = ModelingToolkitBase.BlockedArray(x, [length(x)])
    for T in (Float16, Float32, Float64, ComplexF16, ComplexF32, ComplexF64)
        p = MTKParameters(T[1], z0, (blocked(T[2]),), (T[3],), (), ())
        @test fast_zero(p)
        check_zeroing(p; label = "$T buffers", fast = true)
    end
    @test fast_zero(MTKParameters([1.0], Float32[2], (blocked([3.0]),), (Float32[4],), (), ()))
    for p in (
            mtkp([1]), mtkp(BigFloat[1]), mtkp(SVector(1.0)), mtkp(view([1.0], 1:1)),
            mtkp([1.0]; constant = ([1],)), mtkp([1.0]; constant = (BitVector([true]),)),
            mtkp([1.0]; constant = ([[2.0]],)), mtkp([1.0]; constant = (blocked([2.0]),)),
            MTKParameters([1.0], Float64[], ([2.0],), (), (), ()),
            MTKParameters([1.0], Float64[], (blocked([2]),), (), (), ()),
            MTKParameters([1.0], Float64[], (blocked([[2.0]]),), (), (), ()),
            MTKParameters([1.0], Float64[], (ModelingToolkitBase.BlockedArray(ones(1, 1), [1], [1]),), (), (), ()),
            MTKParameters([1.0], Float64[], (), (), ([1.0],), ()),
            MTKParameters([1.0], Float64[], (), (), (), ([1.0],)),
        )
        @test !fast_zero(p)
    end
    # at most 16 discrete and constant buffers (MTK builds one per parameter type),
    # counted before the buffers are classified
    p16 = MTKParameters([1.0], Float64[], ntuple(_ -> blocked([1.0]), 8), ntuple(_ -> [1.0, 2.0], 8), (), ())
    @test fast_zero(p16)
    check_zeroing(p16; label = "16 buffers", fast = true)
    @test allocations(Enzyme.remake_zero!, Enzyme.make_zero(p16)) == 0
    @test !fast_zero(mtkp([1.0]; constant = ntuple(_ -> [1.0, 2.0], 17)))
    @test !fast_zero(mtkp([1.0]; constant = ntuple(_ -> [1.0], 1024)))
end

# an `inactive_type` rule declared after `remake_zero!` ran for that type takes effect:
# the type leaves the fast path and is then zeroed exactly as by Enzyme's generic method
# (which may still use its own earlier activity result for the type)
const P_LATE = MTKParameters(Float16[7], Float16[], (), (), (), ())
Enzyme.remake_zero!(Enzyme.make_zero(P_LATE))
const KNOWN_BEFORE = ModelingToolkitBase.enzyme_fast_zero(typeof(P_LATE))
EnzymeCore.EnzymeRules.inactive_type(::Type{typeof(P_LATE)}) = true
@testset "inactive_type declared after the first call" begin
    @test KNOWN_BEFORE
    @test !ModelingToolkitBase.enzyme_fast_zero(typeof(P_LATE))
    for zero! in (Enzyme.make_zero!, Enzyme.remake_zero!)
        s, e = deepcopy(P_LATE), deepcopy(P_LATE)
        zero!(s)
        invoke(zero!, Tuple{Any, Any}, e, Base.IdSet{Any}())
        @test same_leaves(s, e)
    end
end

# a rule declared after the fast path was compiled for a type invalidates that compiled code
# (the fast path allocates nothing; the generic walk, which still visits the rest of the
# object, does)
const P_LATE2 = MTKParameters(ComplexF16[1], ComplexF16[], (), (ComplexF16[2],), (), ())
const S_LATE2 = Enzyme.make_zero(P_LATE2)
const ALLOCS_BEFORE = allocations(Enzyme.remake_zero!, S_LATE2)
EnzymeCore.EnzymeRules.inactive_type(::Type{Tuple{Vector{ComplexF16}}}) = true
@testset "compiled fast path invalidated by a later rule" begin
    @test ALLOCS_BEFORE == 0
    @test allocations(Enzyme.remake_zero!, S_LATE2) > 0
end

# rules that Enzyme's activity analysis does not consult are not queried either: those of
# array types (classified by their element type) and of ghost types (an empty tuple, a
# tuple of singletons, an `MTKParameters` without data); rejected layouts query no rules
struct GhostMarker end
const P_GHOST = MTKParameters(SVector{0, Float64}(), SVector{0, Float64}(), (), (), (), ())
const QUERIED_RULES = (Vector{ComplexF32}, Tuple{}, Tuple{GhostMarker}, typeof(P_GHOST))
for T in QUERIED_RULES
    @eval EnzymeCore.EnzymeRules.inactive_type(::Type{$T}) = error("queried")
end
@testset "no rule queries Enzyme does not make" begin
    for p in (
            MTKParameters(ComplexF32[7], SVector{0, Float64}(), (), (), (), ()),
            MTKParameters([7.0], Float64[], (GhostMarker(),), (), (), ()),
        )
        for zero! in (Enzyme.make_zero!, Enzyme.remake_zero!)
            s, e = deepcopy(Enzyme.make_zero(p)), deepcopy(Enzyme.make_zero(p))
            s.tunable[1] = e.tunable[1] = 7
            zero!(s)
            invoke(zero!, Tuple{Any, Any}, e, Base.IdSet{Any}())
            @test s.tunable == e.tunable && iszero(s.tunable)
        end
    end
    @test ModelingToolkitBase.enzyme_fast_zero(typeof(P_GHOST))
    for zero! in (Enzyme.make_zero!, Enzyme.remake_zero!)
        @test zero!(P_GHOST) === invoke(zero!, Tuple{Any, Any}, P_GHOST, Base.IdSet{Any}())
    end
end
# remove the rules for types this file does not own, which later tests in this process
# would otherwise see
for T in (
        QUERIED_RULES..., Tuple{Vector{Float32}, Vector{ComplexF32}}, typeof(P_LATE),
        Tuple{Vector{ComplexF16}},
    )
    Base.delete_method(which(EnzymeCore.EnzymeRules.inactive_type, Tuple{Type{T}}))
end
