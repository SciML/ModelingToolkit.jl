import SymbolicUtils
import SymbolicUtils: FnType

# ModelingToolkit -> SymbolicUtils
SymbolicUtils.istree(x::Operation) = true
function SymbolicUtils.operation(x::Operation)
    if x.op isa Variable
        T = FnType{NTuple{length(x.args), Any}, vartype(x.op)}
        SymbolicUtils.Variable{T}(x.op.name)
    else
        x.op
    end
end

# This is required to infer the right type for
# Operation(Variable{Parameter{Number}}(:foo), [])
# While keeping the metadata that the variable is a parameter.
SymbolicUtils.promote_symtype(f::SymbolicUtils.Variable{FnType{X,Parameter{Y}}},
                              xs...) where {X, Y} = Y

SymbolicUtils.arguments(x::Operation) = x.args

# SymbolicUtils wants raw numbers
SymbolicUtils.to_symbolic(x::Constant) = x.value
SymbolicUtils.to_symbolic(x::Variable{T}) where {T} = SymbolicUtils.Variable{T}(x.name)

# Optional types of vars
# Once converted to SymbolicUtils Variable, a Parameter needs to hide its metadata
_vartype(x::Variable{<:Parameter{T}}) where {T} = T
_vartype(x::Variable{T}) where {T} = T
SymbolicUtils.symtype(x::Variable) = _vartype(x) # needed for a()
SymbolicUtils.symtype(x::SymbolicUtils.Variable{<:Parameter{T}}) where {T} = T

# returning Any causes SymbolicUtils to infer the type using `promote_symtype`
# But we are OK with Number here for now I guess
SymbolicUtils.symtype(x::Expression) = Number


# SymbolicUtils -> ModelingToolkit

function simplify_constants(expr)
    SymbolicUtils.simplify(expr) |> to_mtk
end

to_mtk(x) = x
to_mtk(v::SymbolicUtils.Variable{T}) where {T} = Variable{T}(nameof(v))
to_mtk(v::SymbolicUtils.Variable{FnType{X,Y}}) where {X,Y} = Variable{Y}(nameof(v))
function to_mtk(expr::SymbolicUtils.Term)
    Operation(to_mtk(SymbolicUtils.operation(expr)),
              map(to_mtk, SymbolicUtils.arguments(expr)))
end
