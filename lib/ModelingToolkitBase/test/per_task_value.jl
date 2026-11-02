using ModelingToolkitBase, Test, Serialization
using ModelingToolkitBase: PerTaskValue, MutableCacheKey
using ModelingToolkitBase: t_nounits as t, D_nounits as D

__ptv_test_init() = Int[]

function wait_for_collection(ref::WeakRef; attempts = 200)
    for _ in 1:attempts
        GC.gc(true)
        yield()
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
        # Read-only accesses by the surviving main task cannot sweep the dead
        # entry; the task's own finalizer must have released its value.
        @test wait_for_collection(payload[])
        @test ptv[] === mainval
        @test length(ptv.dict) == 1
    end

    @testset "dropped owner releases values" begin
        payload = Ref{WeakRef}()
        _own_payload(payload)
        # The cleanup finalizer registered on this task only holds a `WeakRef`
        # to `ptv`, so it must not keep the owner alive.
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
            # the surviving `sys` retains the value only for live tasks
            @test wait_for_collection(payload[])
        end
    end

    @testset "deepcopy" begin
        ptv = PerTaskValue{Vector{Int}}(__ptv_test_init)
        push!(ptv[], 1)

        cp = deepcopy(ptv)
        @test cp !== ptv
        @test cp[] == Int[]
        push!(cp[], 9)
        @test ptv[] == [1]
        @test cp[] == [9]

        # shared references stay shared in the copy
        cps = deepcopy([ptv, ptv])
        @test cps[1] === cps[2]
        @test cps[1] !== ptv
    end

    @testset "serialization" begin
        ptv = PerTaskValue{Vector{Int}}(__ptv_test_init)
        push!(ptv[], 1)

        io = IOBuffer()
        serialize(io, ptv)
        cp = deserialize(seekstart(io))
        @test cp isa PerTaskValue{Vector{Int}}
        @test cp !== ptv
        @test cp[] == Int[]
        push!(cp[], 9)
        @test ptv[] == [1]

        # shared references stay shared in the result
        io = IOBuffer()
        serialize(io, [ptv, ptv])
        cps = deserialize(seekstart(io))
        @test cps[1] === cps[2]
        @test cps[1] !== ptv
    end
end
