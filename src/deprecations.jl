@deprecate structural_simplify(sys; kwargs...) mtkcompile(sys; kwargs...)
@deprecate structural_simplify(sys, io; kwargs...) mtkcompile(
    sys; inputs = io[1], outputs = io[2], kwargs...)

macro mtkbuild(exprs...)
    return quote
        Base.depwarn("`@mtkbuild` is deprecated. Use `@mtkcompile` instead.", :mtkbuild)
        @mtkcompile $(exprs...)
    end |> esc
end
