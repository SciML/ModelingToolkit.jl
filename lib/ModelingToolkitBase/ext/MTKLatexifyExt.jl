module MTKLatexifyExt

using ModelingToolkitBase
import SymbolicIndexingInterface: getname
using Latexify

@latexrecipe function f(c::ModelingToolkitBase.Connection)
    index --> :subscript
    fn = eltype(c.systems) <: ModelingToolkitBase.AbstractSystem ? nameof : getname
    return Expr(:call, :connect, map(fn, c.systems)...)
end

function Base.show(io::IO, ::MIME"text/latex", ap::ModelingToolkitBase.Connection)
    return print(io, latexify(ap))
end

@latexrecipe function f(sys::ModelingToolkitBase.AbstractSystem)
    return latexify(equations(sys))
end

function Base.show(io::IO, ::MIME"text/latex", x::ModelingToolkitBase.AbstractSystem)
    return print(io, "\$\$ " * latexify(x) * " \$\$")
end

@latexrecipe function f(ap::ModelingToolkitBase.AnalysisPoint)
    index --> :subscript
    snakecase --> true
    ap.input === nothing && return 0
    outs = Expr(:vect)
    append!(outs.args, ModelingToolkitBase.ap_var.(ap.outputs))
    return Expr(:call, :AnalysisPoint, ModelingToolkitBase.ap_var(ap.input), ap.name, outs)
end

function Base.show(io::IO, ::MIME"text/latex", ap::ModelingToolkitBase.AnalysisPoint)
    return print(io, latexify(ap))
end

end
