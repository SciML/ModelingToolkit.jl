using ModelingToolkitStandardLibrary.Electrical, ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkitStandardLibrary.Electrical: _and, _or, _not, _xor
using ModelingToolkitStandardLibrary.Electrical: U, X, F0, F1, Z, W, L, H, DC, Uninitialized
using ModelingToolkitStandardLibrary.Electrical: AndTable, OrTable, NotTable, XorTable
using ModelingToolkitStandardLibrary.Electrical: get_logic_level
using OrdinaryDiffEq: ReturnCode.Success

# using ModelingToolkitStandardLibrary.Electrical: Set, Reset

@testset "Logic, logic-vectors and helpers" begin
    # Logic and helper functions
    @test length(instances(Logic)) == 9
    @test convert.(Logic, [1, 0]) |> typeof == Vector{Logic}
    @test get_logic_level(Z) == 5

    io = IOBuffer()
    show(io, MIME("text/plain"), Uninitialized)
    @test String(take!(io)) == "U"

    # Logic zeros and ones
    @test zero(Logic) == zero(U) == F0
    @test one(Logic) == one(U) == F1
    @test ones(Logic, 2, 2) == [F1 F1
                                F1 F1]

    # Logic vectors
    u_logic = StdULogicVector([U, W, X, 1])
    @test typeof(u_logic.logic) == Vector{Logic}
    @test get_logic_level(u_logic) == [1, 6, 2, 4]

    logic = StdLogicVector([U, W, X, 1])
    @test typeof(logic.logic) == Vector{Logic}
    @test get_logic_level(logic) == [1, 6, 2, 4]

    # Predefined logic vectors
    @test std_ulogic.logic == [U, X, F0, F1, Z, W, L, H, DC]
    @test UX01.logic == [U, X, F0, F1]
    @test UX01Z.logic == [U, X, F0, F1, Z]
    @test X01.logic == [X, F0, F1]
    @test X01Z.logic == [X, F0, F1, Z]

    # Logic vector helpers
    test_logic_matrix = StdULogicVector([U F0
                                         F1 X])
    test_logic_vector = StdLogicVector([U, F0, F1, X])

    size(test_logic_matrix) == (2, 2)
    axes(test_logic_matrix) == (Base.OneTo(2), Base.OneTo(2))

    getindex(test_logic_matrix, 1, 1) == U
    getindex(test_logic_vector, 1) == U

    setindex!(test_logic_matrix, Z, 1, 1)
    @test test_logic_matrix[1, 1] == Z
    setindex!(test_logic_vector, Z, 1)
    @test test_logic_vector[1] == Z

    # Logic helper functions
    @test get_logic_level.([U, X, F0, F1, Z, W, L, H, DC]) == 1:9
    @test convert.(Logic, [1, 0, U]) == [F1, F0, U]
    @test_throws "3 isn't a valid `Logic` value" convert(Logic, 3)
end

@testset "Logic Tables" begin
    # LogicTable vec-or-mat and helpers
    test_not_logic_table = LogicTable([U, X, F1, F0, X, X, F1, F0, X])
    @test test_not_logic_table[1] == U
    @test test_not_logic_table[F1] == F0

    test_not_logic_table[1] = X
    @test test_not_logic_table[1] == X

    @test_throws ArgumentError LogicTable([U; U])
end

@testset "Gate tables and logic gate helpers" begin
    # logic tables and logic gate helpers
    @test size(AndTable) == size(OrTable) == size(XorTable) == (9, 9)
    @test size(NotTable) == (9,)

    # tests (Number, Number), (Logic, Logic), (Logic, Number) inputs
    @test _and(1, 1, U, W, 1, 0) == F0
    @test _or(0, 1, U, 1) == F1
    @test _xor(0, 1, U, U, 1, 1) == U
    # tests (Number, Logic) input
    @test _and(1, F1) == F1
    @test _or(0, F0) == F0
    @test _xor(1, F0) == F1
    # tests Number and Logic (via internal convert)
    @test _not(1) == F0

    A = [U, W, Z, F1, F0]
    B = [U, W, X, F0, DC]
    _xor.(A, B) == _and.(_or.(A, _not.(B)), _or.(_not.(A), B))
end

#=

@named set1 = Set()
@named reset1 = Reset()
@named set2 = Set()
@named reset2 = Reset()
@named out = DigitalPin()

@testset "Not gate" begin
    @named set = Set()
    @named reset = Reset()
    @named not = Not()
    sources = [set, reset]
    for source in sources
        not_eqs = [connect(source.d, not.x)
                   connect(out, not.y)]
        @named not_model = System(not_eqs, t, systems = [out, not, source])
        sys = alias_elimination(not_model)
        u0 = [
            not.y.val => 0.0,
        ]
        prob = ODEProblem(sys, u0, (0, 1.5))
        sol = solve(prob, Rosenbrock23())

        @test .!(sol[not.x.val] .> 0.5) == Bool.(sol[not.y.val])
    end
end

@testset "And, Nand gate" begin
    @named and = And()
    @named nand = Nand()
    for one in [set1, reset1], two in [set2, reset2]
        and_eqs = [connect(one.d, and.x1)
                   connect(two.d, and.x2)
                   connect(out, and.y)]
        @named and_model = System(and_eqs, t, systems = [and, one, two, out])
        sys = alias_elimination(and_model)

        u0 = []
        prob = ODEProblem(sys, u0, (0, 1.5))
        sol = solve(prob, Rosenbrock23())
        @test sol.retcode == Success
        @test sol[and.y.val] == _and.(sol[one.d.val], sol[two.d.val])

        nand_eqs = [connect(one.d, nand.x1)
                    connect(two.d, nand.x2)
                    connect(out, nand.y)]
        @named nand_model = System(nand_eqs, t, systems = [nand, one, two, out])
        sys = alias_elimination(nand_model)

        u0 = []
        prob = ODEProblem(sys, u0, (0, 1.5))
        nsol = solve(prob, Rosenbrock23())
        @test nsol.retcode == Success
        @test nsol[nand.y.val] == _not.(sol[and.y.val])
    end
end

@testset "Or gate" begin
    @named or = Or()
    @named nor = Nor()
    for one in [set1, reset1], two in [set2, reset2]
        or_eqs = [connect(one.d, or.x1)
                  connect(two.d, or.x2)
                  connect(out, or.y)]
        @named or_model = System(or_eqs, t, systems = [or, one, two, out])
        sys = alias_elimination(or_model)

        u0 = []
        prob = ODEProblem(sys, u0, (0, 1.5))
        sol = solve(prob, Rosenbrock23())
        @test sol.retcode == Success
        @test sol[or.y.val] == _or.(sol[one.d.val], sol[two.d.val])

        nor_eqs = [connect(one.d, nor.x1)
                   connect(two.d, nor.x2)
                   connect(out, nor.y)]
        @named nor_model = System(nor_eqs, t, systems = [nor, one, two, out])
        sys = alias_elimination(nor_model)

        u0 = []
        prob = ODEProblem(sys, u0, (0, 1.5))
        nsol = solve(prob, Rosenbrock23())
        @test nsol.retcode == Success
        @test nsol[nor.y.val] == _not.(sol[or.y.val])
    end
end

@testset "Xor gate" begin
    @named xor = Xor()
    @named xnor = Xnor()
    for one in [set1, reset1], two in [set2, reset2]
        xor_eqs = [connect(one.d, xor.x1)
                   connect(two.d, xor.x2)
                   connect(out, xor.y)]
        @named xor_model = System(xor_eqs, t, systems = [xor, one, two, out])
        sys = alias_elimination(xor_model)

        u0 = []
        prob = ODEProblem(sys, u0, (0, 1.5))
        sol = solve(prob, Rosenbrock23())
        @test sol.retcode == Success
        @test sol[xor.y.val] ==
              _or.(_and.(sol[one.d.val], _not.(sol[two.d.val])),
                   _and.(sol[two.d.val], _not.(sol[one.d.val])))

        xnor_eqs = [connect(one.d, xnor.x1)
                    connect(two.d, xnor.x2)
                    connect(out, xnor.y)]
        @named xnor_model = System(xnor_eqs, t, systems = [xnor, one, two, out])
        sys = alias_elimination(xnor_model)

        u0 = []
        prob = ODEProblem(sys, u0, (0, 1.5))
        nsol = solve(prob, Rosenbrock23())
        @test nsol.retcode == Success
        @test nsol[xnor.y.val] == _not.(sol[xor.y.val])
    end
end

@testset "Half Adders" begin
    @named ha = HalfAdder()
    @named out1 = DigitalPin()
    @named out2 = DigitalPin()
    for one in [set1, reset1], two in [set2, reset2]
        ha_eqs = [connect(one.d, ha.x1)
                  connect(two.d, ha.x2)
                  connect(out1, ha.y0)
                  connect(out2, ha.y1)]
        @named ha_model = System(ha_eqs, t, systems = [ha, one, two, out1, out2])
        sys = alias_elimination(ha_model)

        u0 = []

        prob = ODEProblem(sys, u0, (0, 1.5))
        sol = solve(prob, Rosenbrock23())

        @test sol.retcode == Success
        @test sol[ha.y0.val] == sol[ha.sum] == _xor.(sol[one.d.val], sol[two.d.val])
    end
end

@testset "Full Adder" begin
    @named set3 = Set()
    @named reset3 = Reset()
    @named out1 = DigitalPin()
    @named out2 = DigitalPin()
    for one in [set1, reset1], two in [set2, reset2], three in [set3, reset3]
        @named fa = FullAdder()
        fa_eqs = [connect(one.d, fa.x1)
                  connect(two.d, fa.x2)
                  connect(three.d, fa.x3)
                  connect(out1, fa.y0)
                  connect(out2, fa.y1)]
        @named fa_model = System(fa_eqs, t,
                                    systems = [fa, one, two, three, out1, out2])
        sys = mtkcompile(fa_model)

        u0 = []

        prob = ODEProblem(sys, u0, (0, 1.5))
        sol = solve(prob, Rosenbrock23())

        @test sol.retcode == Success
        @test sol[fa.y0.val] == sol[fa.sum] ==
              _xor(sol[one.d.val], sol[two.d.val], sol[three.d.val])
    end
end

@testset "Multiplexers" begin
    @named mux2x1 = MUX(N = 2)
    @named mux4x1 = MUX(N = 4)
    @named demux1x2 = DEMUX(N = 2)
    out = map(0:7) do i
        DigitalPin(; name = Symbol(:out, i))
    end
    input = map(0:7) do i
        rand([Set(name = Symbol(:s, i)), Reset(name = Symbol(:r, i))])
    end
    @named select0 = Set()
    @named select1 = Reset()

    @info "Building a 2:1 MUX and 1:2 DEMUX..."
    mux_eqs = [connect(input[1].d, mux2x1.d0)
               connect(input[2].d, mux2x1.d1)
               connect(select0.d, mux2x1.s0)
               connect(mux2x1.y, demux1x2.d)
               connect(select1.d, demux1x2.s0)
               connect(out[1], demux1x2.y0)
               connect(out[2], demux1x2.y1)]
    @named mux_model = System(mux_eqs, t,
                                 systems = [mux2x1, demux1x2,
                                     out[2], out[1], select0, select1,
                                     input[2], input[1]])
    sys = alias_elimination(mux_model)

    u0 = []
    prob = ODEProblem(sys, u0, (0, 1.5))
    sol = solve(prob, Rosenbrock23())

    @test sol.retcode == Success
    @test sol[mux2x1.y.val] == sol[mux2x1.d1.val]
    @test sol[mux2x1.y.val] == sol[demux1x2.y0.val]

    @info "Building a 4:1 MUX..."
    @named seta = Set()
    @named setb = Set()
    @named reseta = Reset()
    @named resetb = Reset()
    select0, select1 = [reseta, seta], [resetb, setb]
    for idx1 in 1:2, idx0 in 1:2
        mux_eqs = [connect(input[1].d, mux4x1.d0)
                   connect(input[2].d, mux4x1.d1)
                   connect(input[3].d, mux4x1.d2)
                   connect(input[4].d, mux4x1.d3)
                   connect(select0[idx0].d, mux4x1.s0)
                   connect(select1[idx1].d, mux4x1.s1)
                   connect(out[1], mux4x1.y)]
        @named mux_model = System(mux_eqs, t,
                                     systems = [mux4x1, out[1], select0[idx0],
                                         select1[idx1],
                                         input[1], input[2], input[3], input[4]])
        sys = alias_elimination(mux_model)

        u0 = []
        prob = ODEProblem(sys, u0, (0, 1.5))
        sol = solve(prob, Rosenbrock23())

        @test sol.retcode == Success
        (idx1 == 1 && idx0 == 1) && @test sol[mux4x1.y.val] == sol[mux4x1.d0.val]
        (idx1 == 1 && idx0 == 2) && @test sol[mux4x1.y.val] == sol[mux4x1.d1.val]
        (idx1 == 2 && idx0 == 1) && @test sol[mux4x1.y.val] == sol[mux4x1.d2.val]
        (idx1 == 2 && idx0 == 2) && @test sol[mux4x1.y.val] == sol[mux4x1.d3.val]
    end
end

@testset "Demultiplexers" begin
    @named demux1x2 = DEMUX(N = 2)
    @named demux1x4 = DEMUX(N = 4)
    @named demux1x8 = DEMUX(N = 8)
    @named input = Set()
    out = map(1:4) do i
        DigitalPin(; name = Symbol(:out, i))
    end

    @named seta = Set()
    @named setb = Set()
    @named reseta = Reset()
    @named resetb = Reset()
    select0, select1 = [reseta, seta], [resetb, setb]
    for idx1 in 1:2, idx0 in 1:2
        @info "Building 1:4 DEMUX..."
        eqs = [connect(input.d, demux1x4.d)
               connect(select0[idx0].d, demux1x4.s0)
               connect(select1[idx1].d, demux1x4.s1)
               connect(out[1], demux1x4.y0)
               connect(out[2], demux1x4.y1)
               connect(out[3], demux1x4.y2)
               connect(out[4], demux1x4.y3)]
        @named demux_model = System(eqs, t,
                                       systems = [demux1x4, select0[idx0], select1[idx1],
                                           input,
                                           out[1], out[2], out[3], out[4]])
        sys = alias_elimination(demux_model)

        u0 = []
        prob = ODEProblem(sys, u0, (0, 1.5))
        sol = solve(prob, Rosenbrock23())

        @test sol.retcode == Success
        (idx1 == 1 && idx0 == 1) && @test sol[demux1x4.d.val] == sol[demux1x4.y0.val]
        (idx1 == 1 && idx0 == 2) && @test sol[demux1x4.d.val] == sol[demux1x4.y1.val]
        (idx1 == 2 && idx0 == 1) && @test sol[demux1x4.d.val] == sol[demux1x4.y2.val]
        (idx1 == 2 && idx0 == 2) && @test sol[demux1x4.d.val] == sol[demux1x4.y3.val]
    end
end

@testset "Encoder and Decoder" begin
    @named enc4x2 = Encoder(N = 4)
    @named dec2x4 = Decoder(n = 2)
    out = map(1:4) do i
        DigitalPin(; name = Symbol(:out, i))
    end

    input = map(0:7) do i
        rand([Set(name = Symbol(:s, i)), Reset(name = Symbol(:r, i))])
    end

    @info "Building a 4:2 Encoder..."
    enc_eqs = [connect(input[1].d, enc4x2.d0)
               connect(input[2].d, enc4x2.d1)
               connect(input[3].d, enc4x2.d2)
               connect(input[4].d, enc4x2.d3)
               connect(out[1], enc4x2.y0)
               connect(out[2], enc4x2.y1)]
    @named enc_model = System(enc_eqs, t,
                                 systems = [enc4x2, out[1], out[2],
                                     input[1], input[2], input[3], input[4]])
    sys = alias_elimination(enc_model)

    u0 = []
    prob = ODEProblem(sys, u0, (0, 1.5))
    sol = solve(prob, Rosenbrock23())

    @test sol.retcode == Success
    @test sol[enc4x2.y0.val] == _or.(sol[input[4].d.val], sol[input[2].d.val])
    @test sol[enc4x2.y1.val] == _or.(sol[input[4].d.val], sol[input[3].d.val])

    @info "Building a 2:4 Decoder..."
    dec_eqs = [connect(input[1].d, dec2x4.d0)
               connect(input[2].d, dec2x4.d1)
               connect(out[1], dec2x4.y0)
               connect(out[2], dec2x4.y1)
               connect(out[3], dec2x4.y2)
               connect(out[4], dec2x4.y3)]
    @named dec_model = System(dec_eqs, t,
                                 systems = [dec2x4, out[1], out[2],
                                     out[3], out[4],
                                     input[1], input[2]])
    sys = alias_elimination(dec_model)

    u0 = []
    prob = ODEProblem(sys, u0, (0, 1.5))
    sol = solve(prob, Rosenbrock23())

    @test sol.retcode == Success
    @test sol[dec2x4.y0.val] ==
          _and.(_not.(sol[input[2].d.val]),
                _not.(sol[input[1].d.val]))
    @test sol[dec2x4.y1.val] == _and.(_not.(sol[input[2].d.val]),
                                      (sol[input[1].d.val]))
    @test sol[dec2x4.y2.val] == _and.((sol[input[2].d.val]),
                                      _not.(sol[input[1].d.val]))
    @test sol[dec2x4.y3.val] == _and.((sol[input[2].d.val]), (sol[input[1].d.val]))
end

@testset "Sources and DigitalPin" begin
    @named out = DigitalPin()
    @named and = And()
    @named setₐ = Set()
    @named resetₐ = Reset()
    @named setᵦ = Set()
    @named resetᵦ = Reset()

    for α in [setₐ, resetₐ], β in [setᵦ, resetᵦ]
        eqs = [connect(α.d, and.x1)
               connect(β.d, and.x2)
               connect(out, and.y)]

        @named pul = System(eqs, t, systems = [α, β, and, out])
        sys = alias_elimination(pul)
        # sys = mtkcompile(pul)
        u0 = []
        prob = ODEProblem(sys, u0, (0, 1.5))
        sol = solve(prob, Rosenbrock23())

        @test sol.retcode == Success
        @test ModelingToolkitStandardLibrary._and.(sol[and.x1.val], sol[and.x2.val]) ==
              sol[and.y.val]
    end

    # Test Pulse
    @named pulse = Pulse()
    @named pulseD = PulseDiff()
    @named not = Not()
    eqs = [connect(pulseD.d, not.x)
           connect(out, not.y)]
    @named pul = System(eqs, t, systems = [pulseD, not, out])
    sys = alias_elimination(pul)
    # sys = mtkcompile(pul)

    u0 = []
    prob = ODEProblem(sys, u0, (0, 1.5))
    # sol = solve(prob, Rosenbrock23())
end

=#
