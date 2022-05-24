import DomainSets: Interval, Ball, infimum, supremum

@deprecate IntervalDomain(a, b) Interval(a, b)
@deprecate CircleDomain() Ball()

# type piracy on Interval for downstream compatibility to be reverted once upgrade is complete
function Base.getproperty(domain::Interval, sym::Symbol)
    if sym === :lower
        @warn "domain.lower is deprecated, use infimum(domain) instead"
        return infimum(domain)
    elseif sym === :upper
        @warn "domain.upper is deprecated, use supremum(domain) instead"
        return supremum(domain)
    else
        return getfield(domain, sym)
    end
end
