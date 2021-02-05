struct PDESystem <: ModelingToolkit.AbstractSystem
  eq
  bcs
  domain
  indvars
  depvars
end

Base.getproperty(sys::PDESystem, sym::Symbol) = getfield(sys, sym)
