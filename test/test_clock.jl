using ModelingToolkit, Symbolics
using ModelingToolkit: Inferred, merge_inferred, merge_domains, get_time_domain, collect_operator_variables, get_eq_domain, vars
using Symbolics: value
using OrdinaryDiffEq


eqstates(eqs) = filter(istree, collect(vars(eqs, op=Nothing)))

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

# function Base.getindex(sol::ODESolution, var::Num, t::AbstractArray, from::Symbol=:right)
#     u = sol[var]
#     map(t) do t
#         ind = findfirst(==(t), sol.t)
#         ind === nothing && throw(ArgumentError("t = $t is not a valid time point in the solution."))
#         from === :right && (ind += 1)
#         ind = min(ind, length(sol))
#         u[ind]
#     end
# end
# function Base.getindex(sol::ODESolution, var::Num, t::Number, from::Symbol=:right)
#     ind = findfirst(==(t), sol.t)
#     ind === nothing && throw(ArgumentError("t = $t is not a valid time point in the solution."))
#     from === :right && (ind += 1)
#     ind = min(ind, length(sol))
#     sol[var][ind]
# end

##
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
    @test get_eq_domain(xd) == d
    @test get_eq_domain(xc) == Continuous()

    # test that the domain is propagated correctly for a bunch of different expressions
    exprs = [x, 2x, x^2, x+t, v, 2v, p]
    eqs = [x ~ 1, x ~ p, x ~ v[1]]

    for ex in [exprs; eqs]
        @show ex

        # when x has no domain set, the domain from above should be returned
        @test get_eq_domain(ex) == Inferred() # we can't possibly know the domain of x
        @test get_eq_domain(ex, Continuous()) == Continuous() 
        # @test_broken get_time_domain(x) == Continuous() # Tbis is currently not supported since the metadata is immutable. Might not be needed

        @test get_eq_domain(ex, i) == i 
        @test get_eq_domain(ex, d) == d
    end

    ## with operators
    s = Shift(t, 1)
    xd = Sample(d)(x)
    xd2 = Sample(d2)(x)

    @test get_eq_domain(xd) == d
    @test get_eq_domain(Hold()(xd)) == Continuous()

    @test get_eq_domain(s(xd)) == d # the domain after a shift is inferred from the shifted variable

    @test_throws ClockInferenceException get_eq_domain(xd ~ xd2)


    ## resample to change sampling rate
    xd = Sample(d)(x)
    xd2 = Sample(d2)(xd) # resample xd

    @test get_eq_domain(xd2) == d2

    ## Inferred clock in sample
    @variables yd(t)
    xd = Sample(d)(x)
    eq = yd ~ Sample(x) + xd
    id = get_eq_domain(eq)
    @test id == d

    # bad hybrid system
    @variables t
    d = Clock(t, 1)
    @variables x(t) xd(t) [timedomain=d]
    eq = Shift(t, 1)(xd) + Shift(t, 3)(xd) ~ Hold(xd) + 1
    @test_throws ClockInferenceException get_eq_domain(eq)
end




@testset "hybrid system" begin
    @info "Testing hybrid system"
    dt = 0.1
    @variables t x(t) y(t) u(t) yd(t) ud(t) r(t)
    @parameters kp
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
    @test get_eq_domain(eqs[1]) == d
    @test get_eq_domain(eqs[3]) == Continuous()
    @test get_eq_domain(eqs[4]) == Continuous()


    @unpack eqmap, varmap = ModelingToolkit.clock_inference(eqs)

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
    @test get_eq_domain(eqs[1]) == d
    @test get_eq_domain(eqs[3]) == d2


    @unpack eqmap, varmap = ModelingToolkit.clock_inference(eqs)

    @test varmap[yd1] == d
    @test varmap[ud1] == d
    @test varmap[yd2] == d2
    @test varmap[ud2] == d2
    @test varmap[r] == Continuous()
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

    k = ShiftIndex(d)
    eq = u(k+1) + u(k+3) ~ 0
    @test isequal(ModelingToolkit.normalize_shifts(eq, Dict()), Shift(t, 1)(u) ~ -Shift(t, -1)(u))


    # test with varying number of shift terms
    @variables t
    d = Clock(t, 1)
    @variables ud(t) [timedomain=d]
    eqs = [
            Shift(t, 3)(ud) ~ Shift(t, 2)(ud) + Shift(t, 1)(ud) - Shift(t, 0)(ud)
    ]

    dvs = eqstates(eqs)
    part = ModelingToolkit.preprocess_hybrid_equations(eqs, dvs)
    new_eqs = part.eqs
    ud_delay = part.vars[length(dvs)+1:end]

    expected_eqs = [
        Difference(t; dt=1, update=true)(ud) ~ Shift(t, 0)(ud) + ud_delay[1] - ud_delay[2]
        Difference(t; dt=1, update=true)(ud_delay[1]) ~ ud
        Difference(t; dt=1, update=true)(ud_delay[2]) ~ ud_delay[1]
    ]
    @test all(isequal.(new_eqs, expected_eqs))


    eqs = [ # intermediate shift term omitted
            Shift(t, 3)(ud) ~ Shift(t, 1)(ud) - Shift(t, 0)(ud)
    ]

    dvs = eqstates(eqs)
    part = ModelingToolkit.preprocess_hybrid_equations(eqs, dvs)
    new_eqs = part.eqs
    ud_delay = part.vars[length(dvs)+1:end]
    expected_eqs = [
        Difference(t; dt=1, update=true)(ud) ~ ud_delay[1] - ud_delay[2]
        Difference(t; dt=1, update=true)(ud_delay[1]) ~ ud
        Difference(t; dt=1, update=true)(ud_delay[2]) ~ ud_delay[1]
    ]
    @test all(isequal.(new_eqs, expected_eqs))


    eqs = [
            Shift(t, 3)(ud) ~ - Shift(t, 0)(ud)
    ]

    dvs = eqstates(eqs)
    part = ModelingToolkit.preprocess_hybrid_equations(eqs, dvs)
    new_eqs = part.eqs
    ud_delay = part.vars[length(dvs)+1:end]
    expected_eqs = [
        Difference(t; dt=1, update=true)(ud) ~ - ud_delay[2]
        Difference(t; dt=1, update=true)(ud_delay[1]) ~ ud
        Difference(t; dt=1, update=true)(ud_delay[2]) ~ ud_delay[1]
    ]
    @test all(isequal.(new_eqs, expected_eqs))

    eqs = [
        ud(k+3) ~ - ud(k)
    ]

    dvs = eqstates(eqs)
    part = ModelingToolkit.preprocess_hybrid_equations(eqs, dvs)
    new_eqs = part.eqs
    ud_delay = part.vars[length(dvs)+1:end]
    expected_eqs = [
        Difference(t; dt=1, update=true)(ud) ~ - ud_delay[2]
        Difference(t; dt=1, update=true)(ud_delay[1]) ~ ud
        Difference(t; dt=1, update=true)(ud_delay[2]) ~ ud_delay[1]
    ]
    @test all(isequal.(new_eqs, expected_eqs))


    eqs = [ # ud(k) not expressed as a shift
        Shift(t, 3)(ud) ~ - ud
    ]

    dvs = eqstates(eqs)
    part = ModelingToolkit.preprocess_hybrid_equations(eqs, dvs)
    new_eqs = part.eqs
    ud_delay = part.vars[length(dvs)+1:end]
    expected_eqs = [
        Difference(t; dt=1, update=true)(ud) ~ - ud_delay[2]
        Difference(t; dt=1, update=true)(ud_delay[1]) ~ ud
        Difference(t; dt=1, update=true)(ud_delay[2]) ~ ud_delay[1]
    ]
    @test all(isequal.(new_eqs, expected_eqs))

    eqs = [ # ud(k) not expressed as a shift
        ud(k+3) ~ - ud
    ]

    dvs = eqstates(eqs)
    part = ModelingToolkit.preprocess_hybrid_equations(eqs, dvs)
    new_eqs = part.eqs
    ud_delay = part.vars[length(dvs)+1:end]
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
    
    dvs = eqstates(eqs)
    part = ModelingToolkit.preprocess_hybrid_equations(eqs, dvs)
    new_eqs = part.eqs
    ud_delay = part.vars[length(dvs)+1:end]
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

    dvs = eqstates(eqs)
    new_eqs = ModelingToolkit.preprocess_hybrid_equations(eqs, dvs).eqs

    @test isequal(new_eqs[1], Differential(t)(x) ~ u - x) 
    @test isequal(new_eqs[2], Difference(t; dt=dt, update=true)(yd) ~ y)
    @test isequal(new_eqs[3], u ~ ud)
    @test isequal(new_eqs[4], y ~ x)
    @test isequal(new_eqs[5], ud ~ -kp*yd) 

    sys = ODESystem(eqs, t, name = :sys)
    sysr = ModelingToolkit.hybrid_simplify(sys)
    sysrr = structural_simplify(sysr)

    prob = ODEProblem(sysrr, [x=>1, kp=>1, yd=>1, ud=>-1.0], (0.0, 4.0))
    sol = solve(prob, Rosenbrock23(), tstops=timevec)
    isinteractive() && plot(sol) |> display

    @test sol(timevec)[ud] ≈ sol_manual(timevec)[ud] 
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
    isinteractive() && plot(sol) |> display
    @test issubset(timevec, sol.t)

    @test sol[u, timevec[2:end]] ≈ 0.5sol[u, timevec[1:end-1]] 

end

#=
Many consequtive DiscreteUpdate introduces a delay, so DiscreteUpdate can be used for Sample only, all other algebraic equations must remain algebraic without operators
=#


## Compose systems together

using ModelingToolkit.BipartiteGraphs
using ModelingToolkit.Graphs


dt = 0.5
@variables t
d = Clock(t, dt)
k = ShiftIndex(d)
timevec = 0:0.1:4

function plant(; name)
    @variables  x(t)=1 u(t)=0 [input=true] y(t)=0 [output=true]
    D = Differential(t)
    eqs = [
        D(x) ~ -x + u
        y ~ x
    ]
    ODESystem(eqs, t; name=name)
end

function filt(; name)
    @variables x(t)=0 u(t)=0 [input=true] y(t)=0 [output=true]
    a = 1/exp(dt)
    eqs = [
        x(k+1) ~ a*x + (1-a)*u(k)
        y ~ x
    ]
    ODESystem(eqs, t, name=name)
end

function controller(kp; name)
    @variables y(t)=0 r(t)=0 ud(t)=0 yd(t)=0
    @parameters kp=kp
    eqs = [
        yd ~ Sample(y)
        ud ~ kp * (r - yd)
    ]
    ODESystem(eqs, t; name=name)
end



@named f = filt()
@named c = controller(1)
@named p = plant()

connections = [
    f.u ~ Sample(d)(-1)#(t >= 1)  # step input
    f.y ~ c.r # filtered reference to controller reference
    Hold(c.ud) ~ p.u # controller output to plant input
    p.y ~ c.y # feedback
]

@named cl = ODESystem(connections, t, systems=[f,c,p])

eqs = equations(cl)
cres = ModelingToolkit.clock_inference(eqs)

@test cres.varmap[f.x ] == Clock(t, 0.5)
@test cres.varmap[p.x ] == Continuous()
@test cres.varmap[p.y ] == Continuous()
@test cres.varmap[c.ud] == Clock(t, 0.5)
@test cres.varmap[c.yd] == Clock(t, 0.5)
@test cres.varmap[c.y ] == Continuous()
@test cres.varmap[f.y ] == Clock(t, 0.5)
@test cres.varmap[f.u ] == Clock(t, 0.5)
@test cres.varmap[p.u ] == Continuous()
@test cres.varmap[c.r ] == Clock(t, 0.5)



eqs = equations(cl)
dvs = collect(vars(eqs))
# part = ModelingToolkit.preprocess_hybrid_equations(eqs, dvs)


sysr = ModelingToolkit.hybrid_simplify(cl, param=true)
cont, disc = ModelingToolkit.get_clocked_partitions(sysr)
@test cont isa ODESystem
@test length(equations(structural_simplify(cont))) == 2 # the diff.eq. and the connection to the discrete input
@test length(disc) == 1




# eqmap, varmap = ModelingToolkit.clock_inference(equations(cl))
# ModelingToolkit.substitute_algebraic_eqs(equations(cl), varmap)

sysrr = structural_simplify(sysr)
prob = ODAEProblem(sysrr, [], (0.0, 4.1))
# prob = ODEProblem(sysr, [], (0.0, 4.1), check_length=false)

##

sol = solve(prob, Rosenbrock23(), tstops=timevec)
@test issubset(timevec, sol.t)
isinteractive() && plot(sol) |> display
@test_broken all(sol[f.u, timevec[2:end]] .== -1) 
# @test sol[u, timevec[2:end]] ≈ 0.5sol[u, timevec[1:end-1]]

#=
QUESTIONS:
Does equation order matter in a discrete partition? No, in modelica
> The order between the equations in a when-equation does not matter

Modelica does impose causality on the equations in a when condition:
> The equations within the when-equation must have one of the following forms:
• v = expr;
• (out1, out2, out3, ...) = function_call_name(in1, in2, ...);
8.3.5.3 https://specification.modelica.org/v3.4/Ch8.html#when-equations

What's the causality here? ud should obviously come after yd, but is that always the case?
yd ~ sample(y)
ud ~ (r - yd)
If we consider ud an algebraic equation, it's always valid, and we do not have to care more about it. If it's a discrete algebraic equation evluated together with everything else in the clock partition, it depends on yd

alias_elimination moves equations to observed, but also messes with the remaining discrete algebraic equations so that 
ud ~ (r - yd)
becomes
0 ~ (r - yd) - ud


It's better if the hybrid_simplify substitutes the rhs of all discrete-algebraic equations into difference eqs where the lhs appears
=#

#=
the discrete states can be replaced by parameters which are updated by the discrete affect function? We still need to handle storage of the discrete states in the solution. push!(sol, stuff)?

1. replace discrete vars with params. Done
2. ODAEProblem needs to write the integrator function using only discrete states
3. difference affect must take care of discrete dynamics
=#