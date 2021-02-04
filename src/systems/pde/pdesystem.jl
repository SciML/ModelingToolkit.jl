struct PDESystem <: ModelingToolkit.AbstractSystem
  eq
  bcs
  domain
  indvars
  depvars
  
  @add_kwonly function PDESystem(eq, bcs, domain, indvars, depvars)
      new(eq, bcs, domain, indvars, depvars)
  end
end
