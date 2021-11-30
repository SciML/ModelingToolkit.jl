using ModelingToolkit, Symbolics
using ModelingToolkit: Inferred, merge_inferred, merge_domains, propagate_time_domain, get_time_domain
using Symbolics: value
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
        @show x
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
    @test propagate_time_domain(eqs[3]) == (Continuous(), Dict(u => Continuous(), ud => Inferred())) 
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