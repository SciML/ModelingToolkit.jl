struct PDESystem <: ModelingToolkit.AbstractSystem
  eqs
  bcs
  domain
  indvars
  depvars
  ps
  default_p
  @add_kwonly function PDESystem(eqs, bcs, domain, indvars, depvars, ps = SciMLBase.NullParameters(), default_p = nothing)
      new(eqs, bcs, domain, indvars, depvars, ps, default_p)
  end
end

Base.getproperty(x::PDESystem, sym::Symbol) = getfield(x, sym)

Base.summary(prob::PDESystem) = string(nameof(typeof(prob)))
function Base.show(io::IO, sys::PDESystem)
  println(io,summary(sys))
  println(io,"Equations: ", sys.eqs)
  println(io,"Boundary Conditions: ", sys.bcs)
  println(io,"Domain: ", sys.domain)
  println(io,"Dependent Variables: ", sys.depvars)
  println(io,"Independent Variables: ", sys.indvars)
  println(io,"Parameters: ", sys.ps)
  println(io,"Default Parameter Values", sys.default_p)
  return nothing
end
