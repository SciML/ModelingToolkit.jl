using ModelingToolkit, Symbolics
using ModelingToolkit: Inferred, merge_inferred, merge_domains, propagate_time_domain, get_time_domain, collect_operator_variables
using Symbolics: value
using OrdinaryDiffEq

# Helper function to extract discrete states from solution in a right or left continuous way
function Base.getindex(sol::ODESolution, var::Num, t::AbstractArray, from::Symbol=:right)
    u = sol[var]
    map(t) do t
        ind = from === :right ? findlast(==(t), sol.t) : findfirst(==(t), sol.t)
        ind === nothing && throw(ArgumentError("t = $t is not a valid time point in the solution."))
        u[ind]
    end
end
function Base.getindex(sol::ODESolution, var::Num, t::Number, from::Symbol=:right)
    ind = from === :right ? findlast(==(t), sol.t) : findfirst(==(t), sol.t)
    ind === nothing && throw(ArgumentError("t = $t is not a valid time point in the solution."))
    sol[var][ind]
end



@variables t x(t)

@testset "merge_inferred" begin
    @info "Testing merge_inferred"

    n = nothing
    i = Inferred()
    d = Clock(t, 1)
    d2 = Clock(t, 2)

    @test merge_inferred(n,n) == n

    @test merge_inferred(n,i) == i
    @test merge_inferred(i,n) == i

    @test merge_inferred(n,d) == d
    @test merge_inferred(d,n) == d

    @test merge_inferred(i,i) == i
    @test merge_inferred(d,d) == d

    @test_throws ClockInferenceException merge_inferred(d,d2)

end


@testset "merge_domains" begin
    @info "Testing merge_domains"

    n = nothing
    i = Inferred()
    id = InferredDiscrete()
    c = Continuous()
    d = Clock(t, 1)
    d2 = Clock(t, 2)

    @test merge_domains(n,n) == n

    @test merge_domains(n,i) == i
    @test merge_domains(i,n) == i

    @test merge_domains(n,d) == d
    @test merge_domains(d,n) == d

    @test merge_domains(c,i,0) == Continuous()
    @test merge_domains(i,c,0) == Continuous()

    @test merge_domains(i,i,0) == i
    @test merge_domains(d,d,0) == d

    @test merge_domains(n,id,0) == id
    @test merge_domains(id,n,0) == id

    @test merge_domains(i,id,0) == id
    @test merge_domains(id,i,0) == id

    @test_throws ClockInferenceException merge_domains(d,c,0)
    @test_throws ClockInferenceException merge_domains(c,d,0)

    @test_throws ClockInferenceException merge_domains(d,d2,0)

end

## Propagate time domain
@testset "propagate_time_domain" begin
    @info "Testing propagate_time_domain"

    @variables t x(t) v[1:2](t)
    @parameters p

    i = Inferred()
    d = Clock(t, 1)
    d2 = Clock(t, 2)


    # explicit time domain
    @variables xd(t) [timedomain=d]
    @variables xc(t) [timedomain=Continuous()]
    @test propagate_time_domain(xd) == (d, Dict(xd => d))
    @test propagate_time_domain(xc) == (Continuous(), Dict(xc => Continuous())) # TODO: se till att få bort -1 etc.

    # test that the domain is propagated correctly for a bunch of different expressions
    exprs = [x, 2x, x^2, x+t, v, 2v, p]
    eqs = [x ~ 1, x ~ p, x ~ v[1]]

    for ex in [exprs; eqs]
        # @show ex

        # when x has no domain set, the domain from above should be returned
        @test propagate_time_domain(ex)[1] == nothing # we can't possibly know the domain of x
        @test propagate_time_domain(ex, Continuous())[1] == Continuous() 
        # @test_broken get_time_domain(x) == Continuous() # Tbis is currently not supported since the metadata is immutable. Might not be needed

        @test propagate_time_domain(ex, i)[1] == i 
        @test propagate_time_domain(ex, d)[1] == d
    end

    ## with operators
    s = Shift(t, 1)
    xd = Sample(t, d)(x)
    xd2 = Sample(t, d2)(x)

    for x in exprs
        # @show x
        @test propagate_time_domain(Differential(t)(x))[1] == Continuous()
    end

    @test propagate_time_domain(xd)[1] == d
    @test propagate_time_domain(Hold()(xd))[1] == Continuous()

    @test propagate_time_domain(s(xd))[1] == d # the domain after a shift is inferred from the shifted variable

    @test_throws ClockInferenceException propagate_time_domain(xd ~ xd2)


    ## resample to change sampling rate
    xd = Sample(t, d)(x)
    xd2 = Sample(t, d2)(xd) # resample xd

    @test propagate_time_domain(xd2)[1] == d2

    ## Inferred clock in sample
    @variables yd(t)
    xd = Sample(t, d)(x)
    eq = yd ~ Sample(t)(x) + xd
    id, vm = propagate_time_domain(eq)
    @test id == d
    @test vm[yd] == d
    @test vm[x] isa ModelingToolkit.UnknownDomain

    # bad hybrid system
    @variables t
    d = Clock(t, 1)
    @variables x(t) xd(t) [timedomain=d]
    eq = Shift(t, 1)(xd) + Shift(t, 3)(xd) ~ Hold(xd) + 1
    @test_throws ClockInferenceException propagate_time_domain(eq)
end




@testset "hybrid system" begin
    @info "Testing hybrid system"
    dt = 0.1
    @variables t x(t) y(t) u(t) yd(t) ud(t) r(t)
    @parameters kp
    k = SampledTime(t, dt)
    D = Differential(t)

    eqs = [
        # controller (time discrete part `dt=0.1`)
        yd ~ Sample(t, dt)(y)
        ud ~ kp * (r - yd)

        # plant (time continuous part)
        u ~ Hold(ud)
        D(x) ~ -x + u
        y ~ x
    ]


    d = Clock(t, dt)
    @test propagate_time_domain(eqs[1]) == (d, Dict(yd => d, y => Inferred()))
    @test propagate_time_domain(eqs[3]) == (Continuous(), Dict(u => Continuous(), ud => InferredDiscrete())) 
    @test propagate_time_domain(eqs[4]) == (Continuous(), Dict(u => Continuous(), x => Continuous())) 


    eqmap, varmap = ModelingToolkit.equation_and_variable_time_domains(eqs)

    @test varmap[yd] == d
    @test varmap[ud] == d
    @test varmap[r] == d
    @test varmap[x] == Continuous()
    @test varmap[y] == Continuous()
    @test varmap[u] == Continuous()

    @test eqmap[1] == d
    @test eqmap[2] == d
    @test eqmap[3] == Continuous()
    @test eqmap[4] == Continuous()
    @test eqmap[5] == Continuous()

end


@testset "multi-rate hybrid system" begin
    @info "Testing multi-rate hybrid system"
    dt = 0.1
    dt2 = 0.2
    @variables t x(t) y(t) u(t) r(t) yd1(t) ud1(t) yd2(t) ud2(t)
    @parameters kp
    k = SampledTime(t, dt)
    D = Differential(t)

    eqs = [
        # controller (time discrete part `dt=0.1`)
        yd1 ~ Sample(t, dt)(y)
        ud1 ~ kp * (Sample(t, dt)(r) - yd1)

        yd2 ~ Sample(t, dt2)(y)
        ud2 ~ kp * (Sample(t, dt2)(r) - yd2)

        # plant (time continuous part)
        u ~ Hold(ud1) + Hold(ud2)
        D(x) ~ -x + u
        y ~ x
    ]


    d = Clock(t, dt)
    d2 = Clock(t, dt2)
    @test propagate_time_domain(eqs[1]) == (d, Dict(yd1 => d, y => Inferred()))
    @test propagate_time_domain(eqs[3]) == (d2, Dict(yd2 => d2, y => Inferred()))


    eqmap, varmap = ModelingToolkit.equation_and_variable_time_domains(eqs)

    @test varmap[yd1] == d
    @test varmap[ud1] == d
    @test varmap[yd2] == d2
    @test varmap[ud2] == d2
    @test varmap[r] == Inferred()
    @test varmap[x] == Continuous()
    @test varmap[y] == Continuous()
    @test varmap[u] == Continuous()

    @test eqmap[1] == d
    @test eqmap[2] == d
    @test eqmap[3] == d2
    @test eqmap[4] == d2
    @test eqmap[5] == Continuous()
    @test eqmap[6] == Continuous()
    @test eqmap[7] == Continuous()

end



@testset "preprocess_hybrid_equations" begin
    @info "Testing preprocess_hybrid_equations"
    @variables t
    d = Clock(t, 1)
    @variables u(t) ud(t) [timedomain=d]

    eq = Shift(t, 1)(u) + Shift(t, 3)(u) ~ 0
    @test isequal(ModelingToolkit.normalize_shifts(eq, Dict()), Shift(t, 1)(u) ~ -Shift(t, -1)(u))

    eq = u ~ Hold(ud) + 1
    @test isequal(ModelingToolkit.strip_operator(eq, Hold), u ~ ud + 1)

    eq = u + Hold(ud) + 1 ~ 0
    @test isequal(ModelingToolkit.strip_operator(eq, Hold), u + ud + 1 ~ 0)

    eq = u ~ Sample(t, 1)(ud) + 1
    @test isequal(ModelingToolkit.strip_operator(eq, Sample), u ~ ud + 1)


    # test with varying number of shift terms
    @variables t
    d = Clock(t, 1)
    @variables ud(t) [timedomain=d]
    eqs = [
            Shift(t, 3)(ud) ~ Shift(t, 2)(ud) + Shift(t, 1)(ud) - Shift(t, 0)(ud)
    ]

    new_eqs, ud_delay = ModelingToolkit.preprocess_hybrid_equations(eqs)
    expected_eqs = [
        Difference(t; dt=1, update=true)(ud) ~ Shift(t, 0)(ud) + ud_delay[1] - ud_delay[2]
        Difference(t; dt=1, update=true)(ud_delay[1]) ~ ud
        Difference(t; dt=1, update=true)(ud_delay[2]) ~ ud_delay[1]
    ]
    @test all(isequal.(new_eqs, expected_eqs))


    eqs = [ # intermediate shift term omitted
            Shift(t, 3)(ud) ~ Shift(t, 1)(ud) - Shift(t, 0)(ud)
    ]

    new_eqs, ud_delay = ModelingToolkit.preprocess_hybrid_equations(eqs)
    expected_eqs = [
        Difference(t; dt=1, update=true)(ud) ~ ud_delay[1] - ud_delay[2]
        Difference(t; dt=1, update=true)(ud_delay[1]) ~ ud
        Difference(t; dt=1, update=true)(ud_delay[2]) ~ ud_delay[1]
    ]
    @test all(isequal.(new_eqs, expected_eqs))


    eqs = [
            Shift(t, 3)(ud) ~ - Shift(t, 0)(ud)
    ]

    new_eqs, ud_delay = ModelingToolkit.preprocess_hybrid_equations(eqs)
    expected_eqs = [
        Difference(t; dt=1, update=true)(ud) ~ - ud_delay[2]
        Difference(t; dt=1, update=true)(ud_delay[1]) ~ ud
        Difference(t; dt=1, update=true)(ud_delay[2]) ~ ud_delay[1]
    ]
    @test all(isequal.(new_eqs, expected_eqs))


    eqs = [ # ud(k) not expressed as a shift
        Shift(t, 3)(ud) ~ - ud
    ]

    new_eqs, ud_delay = ModelingToolkit.preprocess_hybrid_equations(eqs)
    expected_eqs = [
        Difference(t; dt=1, update=true)(ud) ~ - ud_delay[2]
        Difference(t; dt=1, update=true)(ud_delay[1]) ~ ud
        Difference(t; dt=1, update=true)(ud_delay[2]) ~ ud_delay[1]
    ]
    @test all(isequal.(new_eqs, expected_eqs))

    # test with more than one variable
    
    @variables t
    d = Clock(t, 1)
    @variables ud(t) [timedomain=d] yd(t) [timedomain=d]
    eqs = [
            Shift(t, 3)(ud) ~ Shift(t, 2)(ud) + Shift(t, 1)(ud) - Shift(t, 0)(ud)
            Shift(t, 3)(yd) ~ Shift(t, 1)(ud)
    ]
    
    new_eqs, ud_delay = ModelingToolkit.preprocess_hybrid_equations(eqs)
    expected_eqs = [
        Difference(t; dt=1, update=true)(ud) ~ Shift(t, 0)(ud) + ud_delay[1] - ud_delay[2]
        Difference(t; dt=1, update=true)(yd) ~ ud_delay[1]
        Difference(t; dt=1, update=true)(ud_delay[1]) ~ ud
        Difference(t; dt=1, update=true)(ud_delay[2]) ~ ud_delay[1]
    ]
    @test all(isequal.(new_eqs, expected_eqs))

end


## Test hybrid system in simulation
@testset "Simulate hybrid system" begin
    @info "Testing Simulate hybrid system"



dt = 0.5
@variables t x(t)=0 y(t)=0 u(t)=0 yd(t)=0 ud(t)=0 
@parameters kp
D = Differential(t)
timevec = 0:0.1:4

## Eliminate algebraic vars manually
# The following test implements a P controller in discrete time with a continuous-time plant. All algebraic variables have been manually eliminated

eqs = [
    ud ~ kp * (0 - Sample(t, dt)(x))
    D(x) ~ -x + Hold(ud)
]

new_eqs, new_vars = ModelingToolkit.preprocess_hybrid_equations(eqs)

sys = ODESystem(eqs, t, name = :sys)
sysr = ModelingToolkit.hybrid_simplify(sys)



prob = ODEProblem(sysr, [x=>1, kp=>1, ud=>-1], (0.0, 4.0))
sol_manual = solve(prob, Rosenbrock23(), tstops=timevec)
plot(sol_manual)
@test issubset(timevec, sol_manual.t) # this will be false if the solver stopped early

dt_timevec = (prob.tspan[1]:dt:prob.tspan[2]) 
@test sol_manual[ud, dt_timevec] ≈ -sol_manual[x, dt_timevec] atol=1e-4 # this tests that there are no additional delays. 

##
# The following test implements a P controller in discrete time with a continuous-time plant. This time, algebraic variables are present

@test ModelingToolkit.is_discrete_algebraic(Shift(t, 0)(ud, true) ~ 1)
@test !ModelingToolkit.is_discrete_algebraic(Shift(t, 1)(ud) ~ 1)


eqs = [
    # controller (time discrete part `dt=0.1`)
    yd ~ Sample(t, dt)(y)
    ud ~ kp * (0 - yd)

    # plant (time continuous part)
    u ~ Hold(ud)
    D(x) ~ -x + u
    y ~ x
]


new_eqs, new_vars = ModelingToolkit.preprocess_hybrid_equations(eqs)

@test isequal(new_eqs[1], Difference(t; dt=dt, update=true)(yd) ~ y)
@test isequal(new_eqs[2], Differential(t)(x) ~ u - x) 
@test isequal(new_eqs[3], ud ~ -kp*yd) 
@test isequal(new_eqs[4], u ~ ud)
@test isequal(new_eqs[5], y ~ x)

sys = ODESystem(eqs, t, name = :sys)
sysr = ModelingToolkit.hybrid_simplify(sys)
sysrr = structural_simplify(sysr)
# sysr2 = alias_elimination(sysr) # built into hybrid_simplify

    prob = ODEProblem(sysrr, [x=>1, kp=>1, yd=>1, ud=>-1.0], (0.0, 4.0))
    sol = solve(prob, Rosenbrock23(), tstops=timevec)
    isinteractive() && plot(sol)

    @test sol(timevec)[ud] ≈ sol_manual(timevec)[ud] # broken due to extra delay introduced by sampling, the manually reduced system only has one Difference operator
    @test issubset(timevec, sol.t) # this will be false if the solver stopped early


end

# test system with discrete states
@testset "System with discrete states" begin
    @info "Testing System with discrete states"

    dt = 0.5
    @variables t
    d = Clock(t, dt)
    @variables u(t)=1 [timedomain=d]
    timevec = (0:dt:4) 

    eqs = [
        Shift(t, 1)(u) ~ 0.5u
    ]
    sys = ODESystem(eqs, t, name = :sys)
    sysr = ModelingToolkit.hybrid_simplify(sys)
    # sysrr = structural_simplify(sysr)
    prob = ODEProblem(sysr, [], (0.0, 4.1))
    sol = solve(prob, Rosenbrock23(), tstops=timevec)
    isinteractive() && plot(sol)
    @test issubset(timevec, sol.t)

    @test sol[u, timevec[2:end]] ≈ 0.5sol[u, timevec[1:end-1]] # broken due to extra delay introduced by sampling, the manually reduced system only has one Difference operator

end

#=
Many consequtive DiscreteUpdate introduces a delay, so DiscreteUpdate can be used for Sample only, all other algebraic equations must remain algebraic without operators
=#