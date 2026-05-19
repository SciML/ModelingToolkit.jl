module MTKTrackerExt

import ModelingToolkitBase as MTKBase
import Tracker

function MTKBase.promote_type_with_nothing(::Type{Tracker.TrackedReal{T}}, x::Tracker.TrackedArray{T}) where {T}
    return Tracker.TrackedReal{T}
end

function MTKBase.promote_with_nothing(::Type{Tracker.TrackedReal{T}}, x::Tracker.TrackedArray{T}) where {T}
    return x 
end

end
