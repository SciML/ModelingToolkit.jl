using ModelingToolkit, SciMLBase, Serialization, Unitful, DiffEqNoiseProcess

@parameters t
@variables x(t)
D = Differential(t)

sys = ODESystem([D(x) ~ -0.5*x], defaults=Dict(x=>1.0))
for prob in [
    eval(ModelingToolkit.ODEProblem{false}(sys, nothing, nothing, SciMLBase.NullParameters())),
    eval(ModelingToolkit.ODEProblemExpr{false}(sys, nothing, nothing, SciMLBase.NullParameters()))
]
    _fn = tempname()

    open(_fn, "w") do f
        serialize(f, prob)
    end

    _cmd = "using ModelingToolkit, Serialization; deserialize(\"$_fn\")"

    run(`$(Base.julia_cmd()) -e $(_cmd)`)
end

include("../examples/rc_model.jl")
io = IOBuffer()
write(io, rc_model)
sys = include_string(@__MODULE__, String(take!(io)))
@test sys == flatten(rc_model)

# test metadata is preserved in `toexpr`
W = WienerProcess(0, 0, 0)
@parameters begin 
    # (t=0), [unit=u"s"]
    t, [unit = u"s"]
    (σ = 28.), [description = "sigma"] 
    (ρ = 10)
    (β = 8 / 3)
end

@variables begin 
    (x(t) = 0), 
            [unit = u"m/s"
            description = "rate of convection"
            noise = W]
    (y(t) = 0), [unit = u"m/s"; connect = Flow]
    (z(t) = 0), [unit = u"m/s"]
end

D = Differential(t)

eqs = [D(x) ~ σ * (y - x),
       D(y) ~ x * (ρ - z) - y,
       D(z) ~ x * y - β * z]

sys = ODESystem(eqs)
ps = parameters(sys)
sigma_metadata = ps[1].metadata
expr = toexpr(sys)
sys2 = eval(expr)
ps2 = parameters(sys2)
sigma_metadata2 = ps2[1].metadata
@test sigma_metadata2 == sigma_metadata # fails if :default is not in push_metadata!

sts2 = states(sys2)
for (i, st) in enumerate(states(sys))
    @test st.metadata == sts2[i].metadata 
end

ps2 = parameters(sys2)
for (i, p) in enumerate(parameters(sys))
    @test p.metadata == ps2[i].metadata 
end

iv = independent_variable(sys)
iv2 = independent_variable(sys2)
@test iv.metadata == iv2.metadata 

# [tests nothing] 
# do we want to recurse into ps, sts, etc to ensure they have equivalent metadata?
io = IOBuffer()
write(io, sys)
sys2 = include_string(@__MODULE__, String(take!(io)))
@test sys2 == flatten(sys) 
