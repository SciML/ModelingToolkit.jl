using ModelingToolkit, AbstractTrees, Test

include("common/rc_model.jl")

io = IOBuffer()
print_tree(io, rc_model)
ser = String(take!(io))
str = """rc_model
      ├─ resistor
      │  ├─ p
      │  └─ n
      ├─ capacitor
      │  ├─ p
      │  └─ n
      ├─ shape
      │  └─ output
      ├─ source
      │  ├─ p
      │  ├─ n
      │  └─ V
      └─ ground
         └─ g
      """
@test strip(ser) == strip(str)
