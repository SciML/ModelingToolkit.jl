struct PDESystem
  eqs
  bcs
  domain
  indvars
  depvars

  @add_kwonly function PDESystem(eqs, bcs, domain, indvars, depvars) 
      new(eqs, bcs, domain, indvars, depvars)
  end
end



