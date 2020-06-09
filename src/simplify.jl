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
# Operation(Variable{Parameter{Real}}(:foo), [])
# While keeping the metadata that the variable is a parameter.
SymbolicUtils.promote_symtype(f::SymbolicUtils.Sym{FnType{X,Parameter{Y}}},
                              xs...) where {X, Y} = Y

SymbolicUtils.arguments(x::Operation) = x.args

# SymbolicUtils wants raw numbers
SymbolicUtils.to_symbolic(x::Constant) = x.value
SymbolicUtils.to_symbolic(x::Variable{T}) where {T} = SymbolicUtils.Sym{T}(x.name)

# Optional types of vars
# Once converted to SymbolicUtils Variable, a Parameter needs to hide its metadata
_vartype(x::Variable{<:Parameter{T}}) where {T} = T
_vartype(x::Variable{T}) where {T} = T
SymbolicUtils.symtype(x::Variable) = _vartype(x) # needed for a()
SymbolicUtils.symtype(x::SymbolicUtils.Sym{<:Parameter{T}}) where {T} = T

# returning Any causes SymbolicUtils to infer the type using `promote_symtype`
# But we are OK with Real here for now I guess
SymbolicUtils.symtype(x::Expression) = Real


# SymbolicUtils -> ModelingToolkit

simplify(expr::Expression) = SymbolicUtils.simplify(expr) |> to_mtk
simplify(x::Equation) = simplify(x.lhs) ~ simplify(x.rhs)
simplify(expr) = expr |> to_mtk

@deprecate simplify_constants(ex) simplify(ex)

to_mtk(x) = x
to_mtk(x::Real) = Constant(x)
to_mtk(v::SymbolicUtils.Sym{T}) where {T} = Variable{T}(nameof(v))
to_mtk(v::SymbolicUtils.Sym{FnType{X,Y}}) where {X,Y} = Variable{Y}(nameof(v))
function to_mtk(expr::SymbolicUtils.Term)
    Operation(to_mtk(SymbolicUtils.operation(expr)),
              map(to_mtk, SymbolicUtils.arguments(expr)))
end
