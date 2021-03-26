using ModelingToolkit, AbstractTrees, Test

include("rc_model.jl")

io = IOBuffer()
print_tree(io, rc_model)
ser = String(take!(io))
str = 
"""rc
├─ resistor
│  ├─ p
│  └─ n
├─ capacitor
│  ├─ p
│  └─ n
├─ source
│  ├─ p
│  └─ n
└─ ground
   └─ g
"""
@test strip(ser) == strip(str)
