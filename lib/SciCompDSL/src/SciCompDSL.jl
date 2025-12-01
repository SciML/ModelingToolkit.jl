module SciCompDSL

using OrderedCollections
using Symbolics
using Symbolics: getname, wrap
using SymbolicUtils: unwrap
import ModelingToolkitBase as MTKBase
using MLStyle
using URIs
using PrecompileTools
using DocStringExtensions
using SymbolicIndexingInterface
import ModelingToolkitBase: observed
using Setfield

let allnames = names(MTKBase; all = true),
    banned_names = Set{Symbol}([:eval, :include, :Variable])

    using_expr = Expr(:using, Expr(:(:), Expr(:., :ModelingToolkitBase)))
    inner_using_expr = using_expr.args[1]

    for name in allnames
        name in banned_names && continue
        startswith(string(name), '#') && continue
        push!(inner_using_expr.args, Expr(:., name))
    end
    @eval SciCompDSL $using_expr
end

using ModelingToolkitBase: COMMON_SENTINEL, COMMON_NOTHING, COMMON_MISSING,
                           COMMON_TRUE, COMMON_FALSE, COMMON_INF

@recompile_invalidations begin
    include("model_parsing.jl")
end

export @mtkmodel

@compile_workload begin
    @mtkmodel __testmod__ begin
        @constants begin
            c = 1.0
        end
        @structural_parameters begin
            structp = false
        end
        if structp
            @variables begin
                x(t) = 0.0, [description = "foo", guess = 1.0]
            end
        else
            @variables begin
                x(t) = 0.0, [description = "foo w/o structp", guess = 1.0]
            end
        end
        @parameters begin
            a = 1.0, [description = "bar"]
            if structp
                b = 2 * a, [description = "if"]
            else
                c
            end
        end
        @equations begin
            x ~ a + b
        end
    end
end

end # module SciCompDSL
