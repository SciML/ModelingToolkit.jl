using ModelingToolkit, OrdinaryDiffEq, Test

@parameters t σ ρ β
@variables x(t) y(t) z(t) F(t) u(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x) + F,
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

aliases = [u ~ x + y - z]
lorenz1 = ODESystem(eqs,inputs=[F],outputs=aliases,name=:lorenz1)
lorenz2 = ODESystem(eqs,inputs=[F],outputs=aliases,name=:lorenz2)

connections = [lorenz1.F ~ lorenz2.u,
               lorenz2.F ~ lorenz1.u]
connected = ODESystem(Equation[],t,[],[],outputs=connections,systems=[lorenz1,lorenz2])
sys = connected

@variables lorenz1₊F lorenz2₊F
@test inputs(connected) == Variable[lorenz1₊F, lorenz2₊F]
@show equations(connected)
@show outputs(connected)

function collapse_inputs!(sys::ModelingToolkit.AbstractSystem)
    outputs = []
    for x in inputs(sys)
        idxs = findall(y->convert(Variable,y.lhs).name == x.name,fulleqs)
        idxs === nothing && error("Not all inputs are connected")
        outputeq = x(sys.iv()) ~ reduce(+,getproperty.(getindex.((fulleqs,),idxs),:rhs))
        push!(outputs,outputeq)
    end
end
