@deprecate structural_simplify(sys; kwargs...) mtkcompile(sys; kwargs...)
@deprecate structural_simplify(sys, io; kwargs...) mtkcompile(
    sys; inputs = io[1], outputs = io[2], kwargs...)

macro mtkbuild(exprs...)
    return quote
        Base.depwarn("`@mtkbuild` is deprecated. Use `@mtkcompile` instead.", :mtkbuild)
        $(Expr(:macrocall, var"@mtkcompile",
            LineNumberNode(@__LINE__, @__FILE__), exprs...))
    end |> esc
end

const ODESystem = IntermediateDeprecationSystem

function IntermediateDeprecationSystem(args...; kwargs...)
    Base.depwarn(
        "`ODESystem(args...; kwargs...)` is deprecated. Use `System(args...; kwargs...) instead`.",
        :ODESystem)

    return System(args...; kwargs...)
end

for T in [:NonlinearSystem, :DiscreteSystem, :ImplicitDiscreteSystem]
    @eval @deprecate $T(args...; kwargs...) System(args...; kwargs...)
end

for T in [:ODEProblem, :DDEProblem, :SDEProblem, :SDDEProblem, :DAEProblem,
    :BVProblem, :DiscreteProblem, :ImplicitDiscreteProblem]
    for (pType, pCanonical) in [
            (AbstractDict, :p),
            (AbstractArray{<:Pair}, :(Dict(p))),
            (AbstractArray, :(isempty(p) ? Dict() : Dict(parameters(sys) .=> p)))
        ],
        (uType, uCanonical) in [
            (Nothing, :(Dict())),
            (AbstractDict, :u0),
            (AbstractArray{<:Pair}, :(Dict(u0))),
            (AbstractArray, :(isempty(u0) ? Dict() : Dict(unknowns(sys) .=> u0)))
        ]

        @eval function SciMLBase.$T(sys::System, u0::$uType, tspan, p::$pType; kw...)
            ctor = string($T)
            uCan = string($(QuoteNode(uCanonical)))
            pCan = string($(QuoteNode(pCanonical)))
            @warn """
            `$ctor(sys, u0, tspan, p; kw...)` is deprecated. Use
            `$ctor(sys, merge($uCan, $pCan), tspan)` instead.
            """
            SciMLBase.$T(sys, merge($uCanonical, $pCanonical), tspan; kw...)
        end
        @eval function SciMLBase.$T{iip}(
                sys::System, u0::$uType, tspan, p::$pType; kw...) where {iip}
            ctor = string($T{iip})
            uCan = string($(QuoteNode(uCanonical)))
            pCan = string($(QuoteNode(pCanonical)))
            @warn """
            `$ctor(sys, u0, tspan, p; kw...)` is deprecated. Use
            `$ctor(sys, merge($uCan, $pCan), tspan)` instead.
            """
            return SciMLBase.$T{iip}(sys, merge($uCanonical, $pCanonical), tspan; kw...)
        end
        @eval function SciMLBase.$T{iip, spec}(
                sys::System, u0::$uType, tspan, p::$pType; kw...) where {iip, spec}
            ctor = string($T{iip, spec})
            uCan = string($(QuoteNode(uCanonical)))
            pCan = string($(QuoteNode(pCanonical)))
            @warn """
            `$ctor(sys, u0, tspan, p; kw...)` is deprecated. Use
            `$ctor(sys, merge($uCan, $pCan), tspan)` instead.
            """
            return $T{iip, spec}(sys, merge($uCanonical, $pCanonical), tspan; kw...)
        end
    end

    for pType in [SciMLBase.NullParameters, Nothing], uType in [Any, Nothing]

        @eval function SciMLBase.$T(sys::System, u0::$uType, tspan, p::$pType; kw...)
            ctor = string($T)
            pT = string($(QuoteNode(pType)))
            @warn """
            `$ctor(sys, u0, tspan, p::$pT; kw...)` is deprecated. Use
            `$ctor(sys, u0, tspan)` instead.
            """
            $T(sys, u0, tspan; kw...)
        end
        @eval function SciMLBase.$T{iip}(
                sys::System, u0::$uType, tspan, p::$pType; kw...) where {iip}
            ctor = string($T{iip})
            pT = string($(QuoteNode(pType)))
            @warn """
            `$ctor(sys, u0, tspan, p::$pT; kw...)` is deprecated. Use
            `$ctor(sys, u0, tspan)` instead.
            """
            return $T{iip}(sys, u0, tspan; kw...)
        end
        @eval function SciMLBase.$T{iip, spec}(
                sys::System, u0::$uType, tspan, p::$pType; kw...) where {iip, spec}
            ctor = string($T{iip, spec})
            pT = string($(QuoteNode(pType)))
            @warn """
            `$ctor(sys, u0, tspan, p::$pT; kw...)` is deprecated. Use
            `$ctor(sys, u0, tspan)` instead.
            """
            return $T{iip, spec}(sys, u0, tspan; kw...)
        end
    end
end

for T in [:NonlinearProblem, :NonlinearLeastSquaresProblem,
    :SCCNonlinearProblem, :OptimizationProblem, :SteadyStateProblem]
    for (pType, pCanonical) in [
            (AbstractDict, :p),
            (AbstractArray{<:Pair}, :(Dict(p))),
            (AbstractArray, :(isempty(p) ? Dict() : Dict(parameters(sys) .=> p)))
        ],
        (uType, uCanonical) in [
            (Nothing, :(Dict())),
            (AbstractDict, :u0),
            (AbstractArray{<:Pair}, :(Dict(u0))),
            (AbstractArray, :(isempty(u0) ? Dict() : Dict(unknowns(sys) .=> u0)))
        ]

        @eval function SciMLBase.$T(sys::System, u0::$uType, p::$pType; kw...)
            ctor = string($T)
            uCan = string($(QuoteNode(uCanonical)))
            pCan = string($(QuoteNode(pCanonical)))
            @warn """
            `$ctor(sys, u0, p; kw...)` is deprecated. Use `$ctor(sys, merge($uCan, $pCan))`
            instead.
            """
            $T(sys, merge($uCanonical, $pCanonical); kw...)
        end
        @eval function SciMLBase.$T{iip}(
                sys::System, u0::$uType, p::$pType; kw...) where {iip}
            ctor = string($T{iip})
            uCan = string($(QuoteNode(uCanonical)))
            pCan = string($(QuoteNode(pCanonical)))
            @warn """
            `$ctor(sys, u0, p; kw...)` is deprecated. Use `$ctor(sys, merge($uCan, $pCan))`
            instead.
            """
            return $T{iip}(sys, merge($uCanonical, $pCanonical); kw...)
        end
        @eval function SciMLBase.$T{iip, spec}(
                sys::System, u0::$uType, p::$pType; kw...) where {iip, spec}
            ctor = string($T{iip, spec})
            uCan = string($(QuoteNode(uCanonical)))
            pCan = string($(QuoteNode(pCanonical)))
            @warn """
            `$ctor(sys, u0, p; kw...)` is deprecated. Use `$ctor(sys, merge($uCan, $pCan))`
            instead.
            """
            return $T{iip, spec}(sys, merge($uCanonical, $pCanonical); kw...)
        end
    end
    for pType in [SciMLBase.NullParameters, Nothing], uType in [Any, Nothing]

        @eval function SciMLBase.$T(sys::System, u0::$uType, p::$pType; kw...)
            ctor = string($T)
            pT = string($(QuoteNode(pType)))
            @warn """
            `$ctor(sys, u0, p::$pT; kw...)` is deprecated. Use `$ctor(sys, u0)` instead
            """
            $T(sys, u0; kw...)
        end
        @eval function SciMLBase.$T{iip}(
                sys::System, u0::$uType, p::$pType; kw...) where {iip}
            ctor = string($T{iip})
            pT = string($(QuoteNode(pType)))
            @warn """
            `$ctor(sys, u0, p::$pT; kw...)` is deprecated. Use `$ctor(sys, u0)` instead
            """
            return $T{iip}(sys, u0; kw...)
        end
        @eval function SciMLBase.$T{iip, spec}(
                sys::System, u0::$uType, p::$pType; kw...) where {iip, spec}
            ctor = string($T{iip, spec})
            pT = string($(QuoteNode(pType)))
            @warn """
            `$ctor(sys, u0, p::$pT; kw...)` is deprecated. Use `$ctor(sys, u0)` instead
            """
            return $T{iip, spec}(sys, u0; kw...)
        end
    end
end

macro brownian(xs...)
    return quote
        Base.depwarn(
            "`@brownian` is deprecated. Use `@brownians` instead", :brownian_macro)
        $(@__MODULE__).@brownians $(xs...)
    end |> esc
end
