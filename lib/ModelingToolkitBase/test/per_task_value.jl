using ModelingToolkitBase, Test, Serialization
using ModelingToolkitBase: PerTaskValue, MutableCacheKey
using ModelingToolkitBase: t_nounits as t, D_nounits as D

__ptv_test_init() = Int[]

# `sweep` is called between collections because `WeakKeyDict` drops dead entries
# lazily: the task's finalizer only marks the dict dirty, and the next dict
# operation (`length` here) removes them.
function wait_for_collection(ref::WeakRef; sweep = () -> nothing, attempts = 200)
    for _ in 1:attempts
        GC.gc(true)
        yield()
        sweep()
        ref.value === nothing && return true
    end
    return false
end

function _worker_access(ptv, payload)
    t = Threads.@spawn begin
        v = ptv[]
        push!(v, 1)
        payload[] = WeakRef(v)
    end
    wait(t)
    return nothing
end

function _own_payload(payload)
    ptv = PerTaskValue{Vector{Int}}(__ptv_test_init)
    payload[] = WeakRef(ptv[])
    return nothing
end

function _worker_touch(ptv, payload)
    t = Threads.@spawn begin
        payload[] = WeakRef(ptv[])
    end
    wait(t)
    return nothing
end

@testset "PerTaskValue" begin
    @testset "per-task values" begin
        ptv = PerTaskValue{Vector{Int}}(__ptv_test_init)
        val = ptv[]
        @test val == Int[]
        @test ptv[] === val
        push!(ptv[], 1)
        @test ptv[] == [1]

        child = fetch(Threads.@spawn ptv[])
        @test child !== val
        @test child == Int[]

        ptv[] = [7]
        @test ptv[] == [7]

        empty!(ptv)
        @test ptv[] == Int[]
        @test ptv[] !== val
    end

    @testset "concurrent task isolation" begin
        ptv = PerTaskValue{Vector{Int}}(__ptv_test_init)
        mainval = ptv[]
        n = 64
        outs = fill(false, n)
        @sync for i in 1:n
            Threads.@spawn begin
                v = ptv[]
                push!(v, i)
                yield()
                outs[i] = v == [i]
            end
        end
        @test all(outs)
        @test ptv[] === mainval
    end

    @testset "collected tasks release their values" begin
        ptv = PerTaskValue{Vector{Int}}(__ptv_test_init)
        mainval = ptv[]
        payload = Ref{WeakRef}()
        _worker_access(ptv, payload)
        # dead entries are swept by the next dict operation after the task's
        # finalizer marks the dict dirty
        @test wait_for_collection(payload[]; sweep = () -> length(ptv.values))
        @test ptv[] === mainval
        @test length(ptv.values) == 1
    end

    @testset "dropped owner releases values" begin
        payload = Ref{WeakRef}()
        _own_payload(payload)
        @test wait_for_collection(payload[])
    end

    @testset "`System` per-task caches" begin
        @variables x(t)
        @named sys = System([D(x) ~ x], t)

        cache_tlv = getmetadata(sys, MutableCacheKey, nothing)
        for tlv in (getfield(sys, :irstructure_tlv), cache_tlv)
            @test tlv isa PerTaskValue
            payload = Ref{WeakRef}()
            _worker_touch(tlv, payload)
            # the surviving `sys` retains the value only until the dead entry is
            # swept by a subsequent dict operation
            @test wait_for_collection(payload[]; sweep = () -> length(tlv.values))
        end
    end

    @testset "`System` deepcopy" begin
        @variables x(t)
        @named sys = System([D(x) ~ x], t)
        tlv = getfield(sys, :irstructure_tlv)
        tlv[]

        cp = deepcopy(sys)
        cp_tlv = getfield(cp, :irstructure_tlv)
        # per-task values are not copied; the copy re-initializes on demand
        @test cp_tlv !== tlv
        @test length(cp_tlv.values) == 0
        cp_tlv[]
        @test length(cp_tlv.values) == 1
        @test length(tlv.values) == 1

        cache_tlv = getmetadata(sys, MutableCacheKey, nothing)
        cp_cache_tlv = getmetadata(cp, MutableCacheKey, nothing)
        @test cp_cache_tlv isa PerTaskValue
        @test cp_cache_tlv !== cache_tlv
        @test length(cp_cache_tlv.values) == 0
    end

    @testset "`System` serialization" begin
        @variables x(t)
        @named sys = System([D(x) ~ x], t)
        tlv = getfield(sys, :irstructure_tlv)
        tlv[]

        io = IOBuffer()
        serialize(io, sys)
        cp = deserialize(seekstart(io))
        @test cp isa System
        cp_tlv = getfield(cp, :irstructure_tlv)
        @test cp_tlv isa PerTaskValue
        @test cp_tlv !== tlv
        @test length(cp_tlv.values) == 0
        cp_tlv[]
        @test length(cp_tlv.values) == 1
    end
end
