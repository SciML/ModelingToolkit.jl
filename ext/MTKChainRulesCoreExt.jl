module MTKChainRulesCoreExt

import ModelingToolkit as MTK
import ChainRulesCore
import ChainRulesCore: Tangent, ZeroTangent, NoTangent, zero_tangent, unthunk

function ChainRulesCore.rrule(::Type{MTK.MTKParameters}, tunables, args...)
    function mtp_pullback(dt)
        dt = unthunk(dt)
        dtunables = dt isa AbstractArray ? dt : dt.tunable
        (NoTangent(), dtunables[1:length(tunables)],
            ntuple(_ -> NoTangent(), length(args))...)
    end
    MTK.MTKParameters(tunables, args...), mtp_pullback
end

function subset_idxs(idxs, portion, template)
    ntuple(Val(length(template))) do subi
        [Base.tail(idx.idx) for idx in idxs if idx.portion == portion && idx.idx[1] == subi]
    end
end

selected_tangents(::NoTangent, _) = ()
selected_tangents(::ZeroTangent, _) = ZeroTangent()
function selected_tangents(
        tangents::AbstractArray{T}, idxs::Vector{Tuple{Int}}) where {T <: Number}
    selected_tangents(tangents, map(only, idxs))
end
function selected_tangents(tangents::AbstractArray{T}, idxs...) where {T <: Number}
    newtangents = copy(tangents)
    view(newtangents, idxs...) .= zero(T)
    newtangents
end
function selected_tangents(
        tangents::AbstractVector{T}, idxs) where {S <: Number, T <: AbstractArray{S}}
    newtangents = copy(tangents)
    for i in idxs
        j, k... = i
        if k == ()
            newtangents[j] = zero(newtangents[j])
        else
            newtangents[j] = selected_tangents(newtangents[j], k...)
        end
    end
    newtangents
end
function selected_tangents(tangents::AbstractVector{T}, idxs) where {T <: AbstractArray}
    newtangents = similar(tangents, Union{T, NoTangent})
    copyto!(newtangents, tangents)
    for i in idxs
        j, k... = i
        if k == ()
            newtangents[j] = NoTangent()
        else
            newtangents[j] = selected_tangents(newtangents[j], k...)
        end
    end
    newtangents
end
function selected_tangents(
        tangents::Union{Tangent{<:Tuple}, Tangent{T, <:Tuple}}, idxs) where {T}
    ntuple(Val(length(tangents))) do i
        selected_tangents(tangents[i], idxs[i])
    end
end

function ChainRulesCore.rrule(
        ::typeof(MTK.remake_buffer), indp, oldbuf::MTK.MTKParameters, idxs, vals)
    if idxs isa AbstractSet
        idxs = collect(idxs)
    end
    idxs = map(idxs) do i
        i isa MTK.ParameterIndex ? i : MTK.parameter_index(indp, i)
    end
    newbuf = MTK.remake_buffer(indp, oldbuf, idxs, vals)
    tunable_idxs = reduce(
        vcat, (idx.idx for idx in idxs if idx.portion isa MTK.SciMLStructures.Tunable))
    disc_idxs = subset_idxs(idxs, MTK.SciMLStructures.Discrete(), oldbuf.discrete)
    const_idxs = subset_idxs(idxs, MTK.SciMLStructures.Constants(), oldbuf.constant)
    nn_idxs = subset_idxs(idxs, MTK.NONNUMERIC_PORTION, oldbuf.nonnumeric)

    pullback = let idxs = idxs
        function remake_buffer_pullback(buf′)
            buf′ = unthunk(buf′)
            f′ = NoTangent()
            indp′ = NoTangent()

            tunable = selected_tangents(buf′.tunable, tunable_idxs)
            discrete = selected_tangents(buf′.discrete, disc_idxs)
            constant = selected_tangents(buf′.constant, const_idxs)
            nonnumeric = selected_tangents(buf′.nonnumeric, nn_idxs)
            oldbuf′ = Tangent{typeof(oldbuf)}(; tunable, discrete, constant, nonnumeric)
            idxs′ = NoTangent()
            vals′ = map(i -> MTK._ducktyped_parameter_values(buf′, i), idxs)
            return f′, indp′, oldbuf′, idxs′, vals′
        end
    end
    newbuf, pullback
end

end
