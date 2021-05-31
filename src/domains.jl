import DomainSets: Interval, Ball, infimum, supremum

@deprecate IntervalDomain(a,b) Interval(a,b)
@deprecate CircleDomain() Ball()

# type piracy on Interval for downstream compatibility to be reverted once upgrade is complete
function Base.getproperty(domain::Interval, sym::Symbol)
    if sym === :lower
        return infimum(domain) # or domain.left also defined by IntervalSets.jl
    elseif sym === :upper
        return supremum(domain) # or domain.right
    else
        return getfield(domain, sym)
    end
end
