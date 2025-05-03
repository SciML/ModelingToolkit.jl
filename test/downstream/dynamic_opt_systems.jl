function build_lotkavolterra(; with_constraint = false)
    @parameters α=1.5 β=1.0 γ=3.0 δ=1.0
    @variables x(..) y(..)
    t = M.t_nounits
    D = M.D_nounits

    eqs = [D(x(t)) ~ α * x(t) - β * x(t) * y(t),
        D(y(t)) ~ -γ * y(t) + δ * x(t) * y(t)]

    tspan = (0.0, 1.0)
    parammap = [α => 1.5, β => 1.0, γ => 3.0, δ => 1.0]

    if with_constraint
        constr = [x(0.6) ~ 3.5, x(0.3) ~ 7.0]
        guess = [x(t) => 4.0, y(t) => 2.0]
        u0map = Pair[]
    else
        constr = nothing
        guess = Pair[]
        u0map = [x(t) => 4.0, y(t) => 2.0]
    end

    @mtkbuild sys = ODESystem(eqs, t; constraints = constr)
    sys, u0map, tspan, parammap, guess
end

function rocket_fft()
    
end

function rocket()
    
end

function cartpole()
    
end
