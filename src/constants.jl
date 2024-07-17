import SymbolicUtils: symtype, term, hasmetadata, issym
struct MTKConstantCtx end

isconstant(x::Num) = isconstant(unwrap(x))
"""
Test whether `x` is a constant-type Sym.
"""
function isconstant(x)
    x = unwrap(x)
    x isa Symbolic && getmetadata(x, MTKConstantCtx, false)
end

"""
    toconstant(s)

Maps the parameter to a constant. The parameter must have a default.
"""
function toconstant(s)
    hasmetadata(s, Symbolics.VariableDefaultValue) ||
        throw(ArgumentError("Constant `$(s)` must be assigned a default value."))
    setmetadata(s, MTKConstantCtx, true)
end

toconstant(s::Num) = wrap(toconstant(value(s)))

"""
$(SIGNATURES)

Define one or more constants.

See also [`@independent_variables`](@ref), [`@parameters`](@ref) and [`@variables`](@ref).
"""
macro constants(xs...)
    Symbolics._parse_vars(:constants,
        Real,
        xs,
        toconstant) |> esc
end

"""
Substitute all `@constants` in the given expression
"""
function subs_constants(eqs)
    consts = collect_constants(eqs)
    if !isempty(consts)
        csubs = Dict(c => getdefault(c) for c in consts)
        eqs = substitute(eqs, csubs)
    end
    return eqs
end
