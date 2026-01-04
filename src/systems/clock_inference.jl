function MTKTearing.input_timedomain(::Sample, _::MTKTearing.IOTimeDomainArgsT = nothing)
    return MTKTearing.InputTimeDomainElT[ContinuousClock()]
end
function MTKTearing.output_timedomain(s::Sample, _::MTKTearing.IOTimeDomainArgsT = nothing)
    return s.clock
end

function MTKTearing.input_timedomain(::Hold, args::MTKTearing.IOTimeDomainArgsT = nothing)
    if args !== nothing
        arg = args[1]
        if MTKTearing.has_time_domain(arg)
            return MTKTearing.InputTimeDomainElT[MTKTearing.get_time_domain(arg)]
        end
    end
    return MTKTearing.InputTimeDomainElT[MTKTearing.InferredDiscrete()] # the Hold accepts any discrete
end

function MTKTearing.output_timedomain(::Hold, _::MTKTearing.IOTimeDomainArgsT = nothing)
    return ContinuousClock()
end
