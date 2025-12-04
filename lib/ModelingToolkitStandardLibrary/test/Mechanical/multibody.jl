using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Mechanical.MultiBody2D
using ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition
using OrdinaryDiffEq
# using Setfield
using Test

@named link1 = Link(; m = 1, l = 10, I = 84, g = -9.807)
@named link2 = Link(; m = 1, l = 10, I = 84, g = -9.807, x1_0 = 10)
@named cart = Mass(; m = 1, s = 0)
# @named force = SineForce(;amp=3e3, freq=15)
@named fixed = Fixed()
# @named m1 = Mass(;m=0.5)
# @named m2 = Mass(;m=0.5)

eqs = [connect(link1.TX1, cart.flange) #, force.flange)
       connect(link1.TY1, fixed.flange)
       connect(link1.TX2, link2.TX1)
       connect(link1.TY2, link2.TY1)]

@named model = System(eqs, t, [], []; systems = [link1, link2, cart, fixed])

sys = mtkcompile(model)
@test length(unknowns(sys)) == 6

# The below code does work...
#=
unset_vars = setdiff(unknowns(sys), keys(ModelingToolkit.defaults(sys)))
prob = ODEProblem(sys, unset_vars .=> 0.0, (0.0, 20), []; jac = true)
sol = solve(prob, Rodas5P())

@test sol[cart.s][end] ≈ 4.767 atol=1e-3

plot(sol, idxs = [cart.s])
=#

#=
using CairoMakie
f = Figure()
a = Axis(f[1,1],xlabel="time [s]", ylabel="cart x pos. [m]")
lines!(a, sol.t, sol[cart.s])
f

function plot_link(sol, sys, tmax)
    tm = Observable(0.0)
    idx = Dict(reverse.(enumerate(unknowns(sys))))

    fig = Figure()
    a = Axis(fig[1,1], aspect=DataAspect(), )
    hidedecorations!(a)
    s = @lift(sol($tm, idxs=[link1.x1, link1.x2, link2.x1, link2.x2, link1.y1, link1.y2, link2.y1, link2.y2]))

    m1x1 = @lift($s[1])
    m1x2 = @lift($s[2])
    m2x1 = @lift($s[3])
    m2x2 = @lift($s[4])

    m1y1 = @lift($s[5])
    m1y2 = @lift($s[6])
    m2y1 = @lift($s[7])
    m2y2 = @lift($s[8])

    sz1 = 0.5
    # lines!(a, [-sz1, sz1, sz1, -sz1, -sz1], @lift([$m1x1, $m1x1, $m1x1+sz1*2, $m1x1+sz1*2, $m1x1]))

    lines!(a, @lift([$m1x1, $m1x2]), @lift([$m1y1, $m1y2]), linewidth=10, color=:blue)
    lines!(a, @lift([$m2x1, $m2x2]), @lift([$m2y1, $m2y2]), linewidth=10, color=:red)

    CairoMakie.ylims!(a, -40, 20)
    CairoMakie.xlims!(a, -20, 40)

    # a = Axis(fig[1, 1], xlabel="time [s]", ylabel="position [m]")
    # lines!(a, sol.t, sol[2,:])
    # lines!(a, sol.t, sol[4,:])

    # scatter!(a, tm, m1x)
    # scatter!(a, tm, m2x)
    # ylims!(a, -60, 30)

    framerate = 30
    timestamps = range(0, tmax, step=1/framerate)

    record(fig, "links.mp4", timestamps;
            framerate = framerate) do t
        tm[] = t
    end

    #=
    CairoMakie.Makie.Record(fig, timestamps; framerate=framerate) do t
        tm[] = t
    end
    =#

    nothing
end

plot_link(sol, sys, 20)

=#
