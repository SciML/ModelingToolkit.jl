module MTKTrackerExt

import ModelingToolkitBase as MTKBase
import Tracker

function MTKBase.promote_type_with_nothing(::Type{Tracker.TrackedReal{T}}, x::Tracker.TrackedArray{T}) where {T}
    return Tracker.TrackedReal{T}
end

function MTKBase.promote_with_nothing(::Type{Tracker.TrackedReal{T}}, x::Tracker.TrackedArray{T}) where {T}
    return x 
end

function MTKBase.__iip_u0_ad_wrapper(u0::Tracker.TrackedVector)
    convert(Vector{eltype(u0)}, u0)
end

end
