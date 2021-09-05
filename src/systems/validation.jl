struct ValidationError <: Exception
    message::String
end

function unitfactor(u,tu) #Just need to make the error type consistent
    try
        return Unitful.convfact(u,tu)
    catch err
        throw(ValidationError("Unable to convert [$tu] to [$u]"))
    end
end

"Throw exception on invalid unit types, otherwise return argument."
function screen_unit(result)
    result isa Symbolic && return result
    result isa Unitful.Unitlike || throw(ValidationError("Unit must be a subtype of Unitful.Unitlike, not $(typeof(result))."))
    result isa Unitful.ScalarUnits || throw(ValidationError("Non-scalar units such as $result are not supported. Use a scalar unit instead."))
    result
end

"""Test unit equivalence.

Example of implemented behavior:
```julia
using ModelingToolkit, Unitful
MT = ModelingToolkit
@parameters γ P [unit = u"MW"] E [unit = u"kJ"] τ [unit = u"ms"]
@test MT.equivalent(u"MW" ,u"kJ/ms") # Understands prefixes
@test !MT.equivalent(u"m", u"cm") # Units must be same magnitude
@test MT.equivalent(MT.get_unit(P^γ), MT.get_unit((E/τ)^γ)) # Handles symbolic exponents
```
"""
equivalent(x,y) = isequal(1*x,1*y)
unitless = Unitful.unit(1)

#For dispatching get_unit
Literal = Union{Sym,Symbolics.ArrayOp,Symbolics.Arr,Symbolics.CallWithMetadata}
Conditional = Union{typeof(ifelse),typeof(IfElse.ifelse)}
Comparison = Union{typeof.([==, !=, ≠, <, <=, ≤, >, >=, ≥])...}

set_unitless(x::Vector) = [_has_unit(y) ? y : SymbolicUtils.setmetadata(y,VariableUnit,unitless) for y in x]

constructunit(x::Num) = constructunit(value(x))
function constructunit(x::Unitful.Quantity)
    return Constant(x.val,Dict(VariableUnit=>Unitful.unit(x)))
end

struct Constant{T,M} <: SymbolicUtils.Symbolic{T}
    val::T
    metadata::M
end

Constant(x) = Constant(x, Dict(VariableUnit => Unitful.unit(x)))
Base.:*(x::Num, y::Unitful.Quantity) = value(x) * y
Base.:*(x::Unitful.Quantity, y::Num) = x * value(y)
Base.show(io::IO, v::Constant) = Base.show(io, v.val)

Unitless = Union{typeof.([exp, log, sinh, asinh, asin,
                                  cosh, acosh, acos,
                                  tanh, atanh, atan,
                                  coth, acoth, acot,
                                  sech, asech, asec,
                                  csch, acsch, acsc])...}
isunitless(f::Unitless) = true

function unitcoerce(u::Unitful.Unitlike, x) 
    st = constructunit(x)
    tu = _get_unit(st)
    output = unitfactor(u, tu) * st
    return SymbolicUtils.setmetadata(output,VariableUnit,u)
end

function constructunit(x) #This is where it all starts
    if _has_unit(x)
        return x
    elseif !SymbolicUtils.istree(x) || operation(x) isa Sym # If a bare symbol doesn't have units, it's unitless
        return SymbolicUtils.setmetadata(x, VariableUnit, unitless)
    else
        op = operation(x)
        if op isa Term
            gp = getmetadata(x, Symbolics.GetindexParent, nothing) # Like x[1](t)
            tu = screen_unit(getmetadata(gp, VariableUnit, unitless))
            return SymbolicUtils.setmetadata(x, VariableUnit, tu)
        end
        args = arguments(x)
        constructunit(op, args)
    end
end

function constructunit(op, args) # Fallback
    if isunitless(op)
        try
            args = unitcoerce.(unitless, args)
            return SymbolicUtils.setmetadata(op(args...), VariableUnit, unitless)
        catch err
            if err isa Unitful.DimensionError
                argunits = get_unit.(args)
                throw(ValidationError("Unable to coerce $args to dimensionless from $argunits for function $op."))
            else
                rethrow(err)
            end
        end
    else
        throw(ValidationError("Unknown function $op supplied with $args with units $argunits"))
    end
end

function constructunit(op::typeof(getindex), subterms) #for symbolic array access
    arr = subterms[1]
    arrunit = _get_unit(arr) #It had better be there!
    output = op(subterms...)
    return SymbolicUtils.setmetadata(output, VariableUnit, arrunit)
end

function uniformize(subterms)
    newterms = Vector{Any}(undef, size(subterms))
    firstunit = nothing
    for (idx, st) in enumerate(subterms)
        if !isequal(st, 0)
            st = constructunit(st)
            tu = _get_unit(st)
            if firstunit === nothing
                firstunit = tu
            end
            newterms[idx] = unitfactor(firstunit, tu) * st
        else
            newterms[idx] = 0
        end
    end
    return newterms
end

function constructunit(op::typeof(+), subterms)
    newterms = uniformize(subterms)
    output = +(newterms...)
    return SymbolicUtils.setmetadata(output, VariableUnit, _get_unit(newterms[1]))
end

Literal = Union{Sym,Symbolics.ArrayOp, Symbolics.Arr, Symbolics.CallWithMetadata}
Conditional = Union{typeof(ifelse), typeof(IfElse.ifelse)}
Comparison = Union{typeof.([==, !=, ≠, <, <=, ≤, >, >=, ≥])...}
Trig = Union{typeof.([sin, cos, tan, sec, csc, cot])...}

function constructunit(op::Trig, subterms)
    arg = constructunit(only(subterms))
    argunit = _get_unit(arg)
    if equivalent(argunit, u"°")
        arg = pi/180 * arg
    end
    return SymbolicUtils.setmetadata(op(arg), VariableUnit, unitless)
end

function constructunit(op::Conditional, subterms)
    newterms = Vector{Any}(undef, 3)
    firstunit = nothing
    newterms[1] = constructunit(subterms[1])
    newterms[2:3] = uniformize(subterms[2:3])
    output = op(newterms...)
    return SymbolicUtils.setmetadata(output, VariableUnit, _get_unit(newterms[2]))
end

function constructunit(op::Union{Differential,Difference}, subterms)
    numerator = constructunit(only(subterms))
    nu = _get_unit(numerator)
    denominator  = op isa Differential ? constructunit(op.x) : constructunit(op.t) #TODO: make consistent!
    du = _get_unit(denominator)
    output = op isa Differential ? Differential(denominator)(numerator) : Difference(denominator)(numerator)
    return SymbolicUtils.setmetadata(output, VariableUnit, nu/du)
end

function constructunit(op::typeof(^), subterms)
    base, exponent = subterms
    base = constructunit(base)
    bu = _get_unit(base)
    exponent = constructunit(exponent)
    exponent = unitfactor(unitless, _get_unit(exponent)) * exponent
    output = base^exponent
    output_unit = bu == unitless ? unitless : (exponent isa Real ? bu^exponent : (1*bu)^exponent)
    return SymbolicUtils.setmetadata(output, VariableUnit, output_unit)
end

function constructunit(op::Comparison, subterms)
    newterms = uniformize(subterms)
    output = op(newterms...)
    return SymbolicUtils.setmetadata(output, VariableUnit, unitless)
end

function constructunit(op::typeof(*), subterms)
    newterms = Vector{Any}(undef, size(subterms))
    pu = unitless
    for (idx, st) in enumerate(subterms)
        st = constructunit(st)
        pu *= _get_unit(st)
        newterms[idx] = st
    end
    output = op(newterms...)
    return SymbolicUtils.setmetadata(output, VariableUnit, pu)
end

#_has_unit(x::Equation) = getmetadata(x,VariableUnit) Doesn't work yet, equations don't have metadata.
_has_unit(x::Real) = true
_has_unit(x::Num) = _has_unit(value(x))
_has_unit(x::Symbolic) = hasmetadata(x,VariableUnit)

_get_unit(x::Real) = unitless
_get_unit(x::Num) = _get_unit(value(x))
_get_unit(x::Symbolic) = screen_unit(getmetadata(x,VariableUnit,unitless))

get_unit(x::Num) = get_unit(value(x))
get_unit(x) = _has_unit(x) ?  _get_unit(x) :  _get_unit(constructunit(x))

function functionize(pt)
    syms = Symbolics.get_variables(pt)
    eval(build_function(constructunit(pt), syms, expression = Val{false}))
end

function constructunit(eq::ModelingToolkit.Equation)
    newterms = uniformize([eq.lhs,eq.rhs])
    return ~(newterms...)
    #return SymbolicUtils.setmetadata(output,VariableUnit,firstunit) #Fix this once Symbolics.jl Equations accept units
end

function validate(jump::Union{ModelingToolkit.VariableRateJump, ModelingToolkit.ConstantRateJump}, t::Symbolic; info::String = "")
    newinfo = replace(info,"eq."=>"jump")
    _validate([jump.rate, 1/t], ["rate", "1/t"], info = newinfo) && # Assuming the rate is per time units
    validate(jump.affect!,info = newinfo)
end

function validate(jump::ModelingToolkit.MassActionJump, t::Symbolic; info::String = "")
    left_symbols = [x[1] for x in jump.reactant_stoch] #vector of pairs of symbol,int -> vector symbols
    net_symbols = [x[1] for x in jump.net_stoch]
    all_symbols = vcat(left_symbols,net_symbols)
    allgood = _validate(all_symbols, string.(all_symbols), info = info)
    n = sum(x->x[2],jump.reactant_stoch,init = 0)
    base_unitful = all_symbols[1] #all same, get first
    allgood && _validate([jump.scaled_rates, 1/(t*base_unitful^n)], ["scaled_rates", "1/(t*reactants^$n))"], info = info)
end

function validate(jumps::ArrayPartition{<:Union{Any, Vector{<:JumpType}}}, t::Symbolic)
    labels = ["in Mass Action Jumps,", "in Constant Rate Jumps,", "in Variable Rate Jumps,"]
    all([validate(jumps.x[idx], t, info = labels[idx]) for idx in 1:3])
end

validate(eq::ModelingToolkit.Equation; info::String = "") = _validate([eq.lhs, eq.rhs], ["left", "right"], info = info)
validate(eq::ModelingToolkit.Equation, term::Union{Symbolic,Unitful.Quantity,Num}; info::String = "") = _validate([eq.lhs, eq.rhs, term], ["left","right","noise"], info = info)
validate(eq::ModelingToolkit.Equation, terms::Vector; info::String = "") = _validate(vcat([eq.lhs, eq.rhs], terms), vcat(["left", "right"], "noise  #".*string.(1:length(terms))), info = info)

"Returns true iff units of equations are valid."
validate(eqs::Vector, noise::Vector; info::String = "") = all([validate(eqs[idx], noise[idx], info = info*" in eq. #$idx") for idx in 1:length(eqs)])
validate(eqs::Vector, noise::Matrix; info::String = "") = all([validate(eqs[idx], noise[idx, :], info = info*" in eq. #$idx") for idx in 1:length(eqs)])
validate(eqs::Vector, term::Symbolic; info::String = "") = all([validate(eqs[idx], term, info = info*" in eq. #$idx") for idx in 1:length(eqs)])
validate(term::Symbolics.SymbolicUtils.Symbolic) = safe_get_unit(term,"") !== nothing

"Throws error if units of equations are invalid."
function rewrite_units(eqs::Vector{Equation})
    output = similar(eqs)
    allgood = true
    for (idx, eq) in enumerate(eqs)
        try
            output[idx] = constructunit(eq)
        catch err
            allgood = false
            err isa ValidationError ? @warn("in eq [$idx], "*err.message) : rethrow(err)
        end
    end
    allgood || throw(ValidationError("Some equations had invalid units. See warnings for details."))
    return output
end

all_dimensionless(states) = all(map(x-> _get_unit(x) in (unitless,nothing),states))
