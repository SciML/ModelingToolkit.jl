
using ModelingToolkit, AbstractTrees, Test

include("components.jl")

open("rc_tree.txt", "w") do io
    print_tree(io, rc_model)
end;

str = 
"""rc_model
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
@test read("rc_tree.txt", String) == str

rm("rc_tree.txt")
