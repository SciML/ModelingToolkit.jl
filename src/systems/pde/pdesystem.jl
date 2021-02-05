struct PDESystem <: ModelingToolkit.AbstractSystem
  eqs
  bcs
  domain
  indvars
  depvars
  
  @add_kwonly function PDESystem(eqs, bcs, domain, indvars, depvars)
      new(eqs, bcs, domain, indvars, depvars)
  end
end

Base.getproperty(x::PDESystem, sym::Symbol) = getfield(x, sym)

Base.summary(prob::PDESystem) = string(nameof(typeof(prob)))
function Base.show(io::IO, sys::PDESystem)
  println(io,summary(sys))
  println(io,"eqs: ", sys.eqs)
  println(io,"bcs: ", sys.bcs)
  println(io,"domain: ", sys.domain)
  println(io,"depvars: ", sys.depvars)
  println(io,"indvars: ", sys.indvars)
  return nothing
end

