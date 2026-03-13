@deprecate structural_simplify(sys; kwargs...) mtkcompile(sys; kwargs...)
@deprecate structural_simplify(sys, io; kwargs...) mtkcompile(
    sys; inputs = io[1], outputs = io[2], kwargs...
)

macro mtkbuild(exprs...)
    return quote
        Base.depwarn("`@mtkbuild` is deprecated. Use `@mtkcompile` instead.", :mtkbuild)
        $(
            Expr(
                :macrocall, var"@mtkcompile",
                LineNumberNode(@__LINE__, @__FILE__), exprs...
            )
        )
    end |> esc
end

const ODESystem = IntermediateDeprecationSystem

function IntermediateDeprecationSystem(args...; kwargs...)
    Base.depwarn(
        "`ODESystem(args...; kwargs...)` is deprecated. Use `System(args...; kwargs...) instead`.",
        :ODESystem
    )

    return System(args...; kwargs...)
end

for T in [:NonlinearSystem, :DiscreteSystem, :ImplicitDiscreteSystem]
    @eval @deprecate $T(args...; kwargs...) System(args...; kwargs...)
end

for T in [
        :ODEProblem, :DDEProblem, :SDEProblem, :SDDEProblem, :DAEProblem,
        :BVProblem, :DiscreteProblem, :ImplicitDiscreteProblem,
    ]
    @eval @fallback_iip_specialize function SciMLBase.$T{iip, spec}(sys::System, u0, tspan, p; kw...) where {iip, spec}
        @warn """
        `$($T)(sys, u0, tspan, p)` is deprecated. Use `$($T)(sys, op, tspan)` instead and provide
        both unknown and parameter values in the operating point `op`.
        """
        if u0 === nothing
            u0 = Dict()
        elseif u0 isa AbstractDict
            u0 = u0
        elseif u0 isa AbstractArray{<:Pair}
            u0 = Dict(u0)
        elseif u0 isa AbstractArray
            u0 = isempty(u0) ? Dict() : Dict(unknowns(sys) .=> u0)
        end
        if p === nothing || p isa SciMLBase.NullParameters
            p = Dict()
        elseif p isa AbstractDict
            p = p
        elseif p isa AbstractArray{<:Pair}
            p = Dict(p)
        elseif p isa AbstractArray
            p = isempty(p) ? Dict() : Dict(parameters(sys) .=> p)
        end
        return SciMLBase.$T{iip, spec}(sys, merge(u0, p), tspan; kw...)
    end
end

for T in [
        :NonlinearProblem, :NonlinearLeastSquaresProblem,
        :SCCNonlinearProblem, :OptimizationProblem, :SteadyStateProblem,
    ]
    @eval @fallback_iip_specialize function SciMLBase.$T{iip, spec}(sys::System, u0, p; kw...) where {iip, spec}
        @warn """
        `$($T)(sys, u0, p)` is deprecated. Use `$($T)(sys, op)` instead and provide
        both unknown and parameter values in the operating point `op`.
        """
        if u0 === nothing
            u0 = Dict()
        elseif u0 isa AbstractDict
            u0 = u0
        elseif u0 isa AbstractArray{<:Pair}
            u0 = Dict(u0)
        elseif u0 isa AbstractArray
            u0 = isempty(u0) ? Dict() : Dict(unknowns(sys) .=> u0)
        end
        if p === nothing || p isa SciMLBase.NullParameters
            p = Dict()
        elseif p isa AbstractDict
            p = p
        elseif p isa AbstractArray{<:Pair}
            p = Dict(p)
        elseif p isa AbstractArray
            p = isempty(p) ? Dict() : Dict(parameters(sys) .=> p)
        end
        return SciMLBase.$T{iip, spec}(sys, merge(u0, p); kw...)
    end
end

macro brownian(xs...)
    return quote
        Base.depwarn(
            "`@brownian` is deprecated. Use `@brownians` instead", :brownian_macro
        )
        $(@__MODULE__).@brownians $(xs...)
    end |> esc
end
