using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D, IfLifting, get_index_cache
using Test
import Symbolics

# Regression tests for hash/objectid-order non-determinism in simplification.
#
# Symbolic hashing (and `objectid` of freshly-created-per-compile types) is not stable
# across sessions / SymbolicUtils versions / GC of the weak hash-cons cache. Several
# pipeline stages used to order observable output (parameter buffer layout, generated
# callbacks/parameters) by iterating `Dict`/`Set` collections keyed by such values, so the
# compiled structure could differ run-to-run within a single session. These tests compile
# the same model many times (with `GC.gc(true)` between runs to provoke hash-cons eviction)
# and assert the compiled structure is identical every time.

const NRUNS = 150

"""Compile `build()` `NRUNS` times with GC stress; return the number of distinct
`dump`-strings observed (should be 1)."""
function ndistinct(build, dump)
    seen = Set{UInt}()
    for _ in 1:NRUNS
        push!(seen, hash(dump(build())))
        GC.gc(true)
    end
    return length(seen)
end

@testset "Deterministic nonnumeric/constant parameter buffer layout (index_cache)" begin
    # Parameters whose symtypes are freshly created on every build mimic the
    # FunctionWrapper/closure parameter types that multibody-style models route into the
    # nonnumeric buffers. The `Dict{Type,Set}` that lays out those buffers must not be
    # iterated in `objectid` order.
    function build()
        @variables x(t) = 1.0
        @parameters k = 1.0
        ps = ModelingToolkit.SymbolicParam[k]
        for i in 1:6
            # a genuinely fresh, empty struct type on every build
            nm = gensym(:NN)
            FT = @eval Main (struct $nm end; $nm)
            p = ModelingToolkit.toparam(Symbolics.unwrap(Symbolics.variable(Symbol(:nn, i); T = FT)))
            push!(ps, ModelingToolkit.setdefault(p, FT.instance))
        end
        return mtkcompile(System([D(x) ~ -k * x], t, [x], ps; name = :sys))
    end
    function dump(ss)
        ic = get_index_cache(ss)
        io = IOBuffer()
        for p in sort(parameters(ss); by = string)
            println(io, p, " => ", ModelingToolkit.parameter_index(ss, p))
        end
        # normalize the fresh gensym type names so only ORDER differences remain
        replace(String(take!(io)), r"##NN#\d+" => "NN")
    end
    @test ndistinct(build, dump) == 1
end

@testset "Deterministic if-lifting (condition vars, parameters, callbacks)" begin
    function build()
        @variables x1(t)=0.1 x2(t)=0.2 x3(t)=0.3 x4(t)=0.4
        @parameters a=1.0 b=2.0 c=3.0
        eqs = [
            D(x1) ~ ifelse(x2 < 0.5, a * x2, -a * x2) + ifelse(x3 > 0.1, 0.1, -0.1)
            D(x2) ~ ifelse(x1 < 0.0, b, -b) + ifelse(x4 < x3, 1.0, 2.0)
            D(x3) ~ ifelse(x2 < x1, c * x1, c * x4) - ifelse(x4 > 0.2, x3, -x3)
            D(x4) ~ ifelse(x3 < 0.3, x1, x2) + ifelse(x1 > x4, 0.5, -0.5)
        ]
        return mtkcompile(System(eqs, t; name = :top); additional_passes = [IfLifting])
    end
    function dump(ss)
        io = IOBuffer()
        for sec in (equations, unknowns, parameters, observed)
            for e in sec(ss)
                println(io, string(e))
            end
        end
        String(take!(io))
    end
    @test ndistinct(build, dump) == 1
end
