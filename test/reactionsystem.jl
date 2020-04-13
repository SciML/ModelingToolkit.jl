using ModelingToolkit

@parameters t a b c
@variables x(t) y(t) z(t)
rxs = [Reaction(a,[x],[y]),
       Reaction(b,[y],[x]),
       Reaction(c,[x,y],[x])]
rs = ReactionSystem(rxs,t,[x,y],[a,b])
odesys = convert(ODESystem,rs)
sdesys = convert(SDESystem,rs)
