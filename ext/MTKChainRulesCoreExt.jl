module MTKChainRulesCoreExt

import ModelingToolkit as MTK
import ChainRulesCore

function ChainRulesCore.rrule(::Type{MTK.MTKParameters}, tunables, args...)
    function mtp_pullback(dt)
        (NoTangent(), dt.tunable[1:length(tunables)], ntuple(_ -> NoTangent(), length(args))...)
    end
    MTK.MTKParameters(tunables, args...), mtp_pullback
end

end
