struct PDESystem
  eq
  bcs
  domain
  indvars
  depvars
  
  @add_kwonly function PDESystem(eq, bcs, domain, indvars, depvars)
      new(eq, bcs, domain, indvars, depvars)
  end
end


struct PDESystem{iip}
  eq
  bcs
  domain
  indvars
  depvars

  @add_kwonly function PDESystem{iip}(eq, bcs, domain, indvars, depvars) where iip
      new{iip}(eq, bcs, domain, indvars, depvars)
  end
end

PDESystem(args...) = PDESystem{true}(args...)

SciMLBase.isinplace(prob::PDESystem{iip}) where iip = iip


