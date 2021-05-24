using DomainSets
struct VarDomainPairing
  variables
  domain::Domain
end
Base.:∈(variable::ModelingToolkit.Num,domain::Domain) = VarDomainPairing(value(variable),domain)
Base.:∈(variable::ModelingToolkit.Num,domain::Interval) = VarDomainPairing(value(variable),domain)
Base.:∈(variables::NTuple{N,ModelingToolkit.Num},domain::Domain) where N = VarDomainPairing(value.(variables),domain)

@deprecate IntervalDomain(a,b) Interval(a,b)
@deprecate CircleDomain() Ball()