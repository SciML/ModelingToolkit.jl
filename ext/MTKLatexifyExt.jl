module MTKLatexifyExt

using ModelingToolkit
import SymbolicIndexingInterface: getname
using Latexify

@latexrecipe function f(c::ModelingToolkit.Connection)
    index --> :subscript
    fn = eltype(c.systems) <: ModelingToolkit.AbstractSystem ? nameof : getname
    return Expr(:call, :connect, map(fn, c.systems)...)
end

function Base.show(io::IO, ::MIME"text/latex", ap::ModelingToolkit.Connection)
    print(io, latexify(ap))
end

@latexrecipe function f(sys::ModelingToolkit.AbstractSystem)
    return latexify(equations(sys))
end

function Base.show(io::IO, ::MIME"text/latex", x::ModelingToolkit.AbstractSystem)
    print(io, "\$\$ " * latexify(x) * " \$\$")
end

@latexrecipe function f(ap::ModelingToolkit.AnalysisPoint)
    index --> :subscript
    snakecase --> true
    ap.input === nothing && return 0
    outs = Expr(:vect)
    append!(outs.args, ModelingToolkit.ap_var.(ap.outputs))
    return Expr(:call, :AnalysisPoint, ModelingToolkit.ap_var(ap.input), ap.name, outs)
end

function Base.show(io::IO, ::MIME"text/latex", ap::ModelingToolkit.AnalysisPoint)
    print(io, latexify(ap))
end

end
