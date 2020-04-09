using ModelingToolkit

@parameters t a b
@variables x(t) y(t) z(t)
rxs = [Reaction(a,[x],[y]),
       Reaction(b,[y],[x])]
rs = ReactionSystem(rxs,t,[x,y],[a,b])
sys = convert(ODESystem,rs)
