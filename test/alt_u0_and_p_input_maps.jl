using ModelingToolkit, Test

### Prepare test equations ###
@parameters t p1 p2 p3
@variables x1(t) x2(t) x3(t)
D = Differential(t)

eqs = [ D(x1) ~ p1 * (x1 - x2),
        D(x2) ~ p2 * (x2 - x3),
        D(x3) ~ p3 * (x3 - x1)]
tspan = (0.0, 100.0)

u0_num = [x1 => 1.0, x2 => 0.5, x3 => 0.0]
u0_sym = [:x1 => 1.0, :x2 => 0.5, :x3 => 0.0]

p_num = [p1 => 0.9, p2 => 1.0, p3 => 1.1]
p_sym = [:p1 => 0.9, :p2 => 1.0, :p3 => 1.1]


### ODEProblem Tests ###
let
    @named osys = ODESystem(eqs)
    u0_sysnum = [osys.x1 => 1.0, osys.x2 => 0.5, osys.x3 => 0.0]
    p_sysnum = [osys.p1 => 0.9, osys.p2 => 1.0, osys.p3 => 1.1]

    prob_num = ODEProblem(osys, u0_num, tspan, p_num)
    prob_sym = ODEProblem(osys, u0_sym, tspan, p_sym)
    prob_sysnum = ODEProblem(osys, u0_sysnum, tspan, p_sysnum)

    @test prob_num.u0 == prob_sym.u0 == prob_sysnum.u0
    @test prob_num.p == prob_sym.p == prob_sysnum.p
end


### SDEProblem Tests ###
let
    noiseeqs = [0.1*x1, 0.1*x2, 0.1*x3]
    @named ssys = SDESystem(eqs,noiseeqs,t,[x1,x2,x3],[p1,p2,p3]; tspan = tspan)
    u0_sysnum = [ssys.x1 => 1.0, ssys.x2 => 0.5, ssys.x3 => 0.0]
    p_sysnum = [ssys.p1 => 0.9, ssys.p2 => 1.0, ssys.p3 => 1.1]

    prob_num = SDEProblem(ssys, u0_num, tspan, p_num)
    prob_sym = SDEProblem(ssys, u0_sym, tspan, p_sym)
    prob_sysnum = SDEProblem(ssys, u0_sysnum, tspan, p_sysnum)

    @test prob_num.u0 == prob_sym.u0 == prob_sysnum.u0
    @test prob_num.p == prob_sym.p == prob_sysnum.p
end


### DAEProblem Tests ###
let
    @named osys = ODESystem(eqs)
    u0_sysnum = [osys.x1 => 1.0, osys.x2 => 0.5, osys.x3 => 0.0]
    p_sysnum = [osys.p1 => 0.9, osys.p2 => 1.0, osys.p3 => 1.1]

    prob_num = DAEProblem(osys, [1.0,2.0,3.0], u0_num, tspan, p_num)
    prob_sym = DAEProblem(osys, [1.0,2.0,3.0], u0_sym, tspan, p_sym)
    prob_sysnum = DAEProblem(osys, [1.0,2.0,3.0], u0_sysnum, tspan, p_sysnum)

    @test prob_num.u0 == prob_sym.u0 == prob_sysnum.u0
    @test prob_num.p == prob_sym.p == prob_sysnum.p
end


### ODEProblem Tests ###
let
    @named osys = ODESystem(eqs)
    u0_sysnum = [osys.x1 => 1.0, osys.x2 => 0.5, osys.x3 => 0.0]
    p_sysnum = [osys.p1 => 0.9, osys.p2 => 1.0, osys.p3 => 1.1]

    # Technically not a real DDE problem ,but sufficient for tests.
    prob_num = DDEProblem(osys, u0_num, tspan, p_num)
    prob_sym = DDEProblem(osys, u0_sym, tspan, p_sym)
    prob_sysnum = DDEProblem(osys, u0_sysnum, tspan, p_sysnum)

    @test prob_num.u0 == prob_sym.u0 == prob_sysnum.u0
    @test prob_num.p == prob_sym.p == prob_sysnum.p
end


### SDDEProblem Tests ###
let
    noiseeqs = [0.1*x1, 0.1*x2, 0.1*x3]
    @named ssys = SDESystem(eqs,noiseeqs,t,[x1,x2,x3],[p1,p2,p3]; tspan = tspan)
    u0_sysnum = [ssys.x1 => 1.0, ssys.x2 => 0.5, ssys.x3 => 0.0]
    p_sysnum = [ssys.p1 => 0.9, ssys.p2 => 1.0, ssys.p3 => 1.1]

    # Technically not a real SDDE problem ,but sufficient for tests.
    prob_num = SDDEProblem(ssys, u0_num, tspan, p_num)
    prob_sym = SDDEProblem(ssys, u0_sym, tspan, p_sym)
    prob_sysnum = SDDEProblem(ssys, u0_sysnum, tspan, p_sysnum)

    @test prob_num.u0 == prob_sym.u0 == prob_sysnum.u0
    @test prob_num.p == prob_sym.p == prob_sysnum.p
end

