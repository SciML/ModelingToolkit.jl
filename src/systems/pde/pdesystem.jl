export PDESystem

struct PDESystem <: ModelingToolkit.AbstractSystem
  eq
  bcs
  domain
  indvars
  depvars
end
