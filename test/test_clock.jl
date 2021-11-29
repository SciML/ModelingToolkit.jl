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

    @test_throws ClockInferenceException merge_domains(c,i,0)
    @test_throws ClockInferenceException merge_domains(i,c,0)

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
    @test propagate_time_domain(xd) == d
    @test propagate_time_domain(xc) == Continuous()

    # test that the domain is propagated correctly for a bunch of different expressions
    exprs = [x, 2x, x^2, x+t, v, 2v, p]
    eqs = [x ~ 1, x ~ p, x ~ v[1]]

    for x in [exprs; eqs]
        # @show x

        # when x has no domain set, the domain from above should be returned
        @test propagate_time_domain(x) == nothing # we can't possibly know the domain of x
        @test propagate_time_domain(x, Continuous()) == Continuous() 
        # @test_broken get_time_domain(x) == Continuous() # Tbis is currently not supported since the metadata is immutable. Might not be needed

        @test propagate_time_domain(x, i) == i 
        @test propagate_time_domain(x, d) == d
    end

    ## with operators
    s = Shift(t, 1)
    xd = Sample(t, d)(x)
    xd2 = Sample(t, d2)(x)

    for x in exprs
        @show x
        @test propagate_time_domain(Differential(t)(x)) == Continuous()
    end

    @test propagate_time_domain(xd) == d
    @test propagate_time_domain(Hold()(xd)) == Continuous()

    @test propagate_time_domain(s(xd)) == d # the domain after a shift is inferred from the shifted variable

    @test_throws ClockInferenceException propagate_time_domain(xd ~ xd2)


    ## resample to change sampling rate
    xd = Sample(t, d)(x)
    xd2 = Sample(t, d2)(xd) # resample xd

    @test propagate_time_domain(xd2) == d2
end
