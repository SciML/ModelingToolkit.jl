using ModelingToolkit, SciMLBase, Serialization

@parameters t
@variables x(t)
D = Differential(t)

@named sys = ODESystem([D(x) ~ -0.5*x], defaults=Dict(x=>1.0))
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
str = String(take!(io))
sys = include_string(@__MODULE__, str)
@test_broken sys == flatten(rc_model) # this actually kind of works, but the variables would have different identities.
