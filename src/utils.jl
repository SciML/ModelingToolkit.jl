function isvariable(x)
    x isa SymbolicT || return false
end
function collect_vars!(unknowns::OrderedSet{SymbolicT}, parameters::OrderedSet{SymbolicT}, eq::Union{Equation, Inequality}, iv::Union{SymbolicT, Nothing};
        depth = 0, op = Symbolics.Operator)
end
function flatten_equations(eqs::Vector{Equation})
    return eqs
end
const JumpType = Union{VariableRateJump, ConstantRateJump, MassActionJump}
