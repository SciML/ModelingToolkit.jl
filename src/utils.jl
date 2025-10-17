function check_parameters(ps, iv)
    for p in ps
    end
end
function is_delay_var(iv::SymbolicT, var::SymbolicT)
    Moshi.Match.@match var begin
    end
end
function check_variables(dvs, iv)
    for dv in dvs
    end
end
function check_lhs(eq::Equation, op, dvs::Set)
    for eq in eqs
    end
end
@noinline function throw_bad_namespacing(systems, idxs)
    if op isa target_op
        x = if target_op <: Differential
        end
    end
end
function setdefault(v, val)
end
function collect_defaults!(defs::SymmapT, vars::Vector{SymbolicT})
    for v in vars
        if def !== nothing
        end
        Moshi.Match.@match v begin
            BSImpl.Term(; f, args) && if f === getindex end => begin
            end
        end
        Moshi.Match.@match v begin
            BSImpl.Term(; f, args) && if f === getindex end => begin
            end
        end
    end
end
function collect_var_to_name!(vars::Dict{Symbol, SymbolicT}, xs::Vector{SymbolicT})
    for x in xs
        x = Moshi.Match.@match x begin
        end
    end
end
@noinline function throw_invalid_operator(opvar, eq, op::Type)
    foreach(expr -> _check_operator_variables(eq, op, expr),
        SymbolicUtils.arguments(expr))
    for eq in eqs
        if length(tmp) == 1
        end
    end
end
function isvariable(x)
    x isa SymbolicT || return false
end
function collect_applied_operators(x, ::Type{op}) where {op}
    if !is_variable_floatingpoint(arg)
    end
end
struct ContinuousOperatorDiscreteArgumentError <: Exception
end
function Base.showerror(io::IO, err::ContinuousOperatorDiscreteArgumentError)
end
struct OperatorIndepvarMismatchError <: Exception
end
function Base.showerror(io::IO, err::OperatorIndepvarMismatchError)
    Moshi.Match.@match ex begin
    end
end
struct OperatorIsAtomic{O} end
function collect_vars!(unknowns::OrderedSet{SymbolicT}, parameters::OrderedSet{SymbolicT}, expr::SymbolicT, iv::Union{SymbolicT, Nothing}; depth = 0, op = Symbolics.Operator)
    Moshi.Match.@match expr begin
    end
end
function collect_vars!(unknowns::OrderedSet{SymbolicT}, parameters::OrderedSet{SymbolicT}, eq::Union{Equation, Inequality}, iv::Union{SymbolicT, Nothing};
        depth = 0, op = Symbolics.Operator)
end
function collect_var!(unknowns::OrderedSet{SymbolicT}, parameters::OrderedSet{SymbolicT}, var::SymbolicT, iv::Union{SymbolicT, Nothing}; depth = 0)
    if Symbolics.iswrapped(var)
        error("""
        """)
    end
    if iscalledparameter(var)
    end
    if scope isa LocalScope
    end
end
function empty_substitutions(sys)
end
function similar_variable(var::BasicSymbolic, name = :anon; use_gensym = true)
    if use_gensym
    end
end
function flatten_equations(eqs::Vector{Equation})
    _eqs = Equation[]
    for eq in eqs
        if !SU.is_array_shape(SU.shape(eq.lhs))
            push!(_eqs, eq)
            continue
        end
        for (l, r) in zip(lhs, rhs)
        end
    end
    return _eqs
end
const JumpType = Union{VariableRateJump, ConstantRateJump, MassActionJump}
