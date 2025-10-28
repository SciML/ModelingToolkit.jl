using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Electrical: OnePort
using OrdinaryDiffEq
using OrdinaryDiffEq: ReturnCode.Success
using IfElse: ifelse

@testset "Chua Circuit" begin
    @component function NonlinearResistor(; name, Ga, Gb, Ve)
        @named oneport = OnePort()
        @unpack v, i = oneport
        pars = @parameters Ga=Ga Gb=Gb Ve=Ve
        eqs = [
            i ~ ifelse(v < -Ve,
            Gb * (v + Ve) - Ga * Ve,
            ifelse(v > Ve,
                Gb * (v - Ve) + Ga * Ve,
                Ga * v))
        ]
        extend(System(eqs, t, [], pars; name = name), oneport)
    end

    @named L = Inductor(L = 18, i = 0.0)
    @named Ro = Resistor(R = 12.5e-3)
    @named G = Conductor(G = 0.565)
    @named C1 = Capacitor(C = 10, v = 4)
    @named C2 = Capacitor(C = 100, v = 0.0)
    @named Nr = NonlinearResistor(Ga = -0.757576,
        Gb = -0.409091,
        Ve = 1)
    @named Gnd = Ground()

    connections = [connect(L.p, G.p)
                   connect(G.n, Nr.p)
                   connect(Nr.n, Gnd.g)
                   connect(C1.p, G.n)
                   connect(L.n, Ro.p)
                   connect(G.p, C2.p)
                   connect(C1.n, Gnd.g)
                   connect(C2.n, Gnd.g)
                   connect(Ro.n, Gnd.g)]

    @named model = System(connections, t, systems = [L, Ro, G, C1, C2, Nr, Gnd])
    sys = mtkcompile(model)
    prob = ODEProblem(sys, Pair[], (0, 5e4), saveat = 0.01)
    sol = solve(prob, Rodas4())

    @test sol.retcode == Success
end
