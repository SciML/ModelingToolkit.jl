using ModelingToolkit, SciMLBase, Serialization

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
