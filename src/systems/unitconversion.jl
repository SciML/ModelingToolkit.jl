"Wrapper for Unitful.convfact that returns a Constant & throws ValidationError instead of DimensionError."
function unitfactor(u, t)
    try
        cf = Unitful.convfact(u, t)
        if cf == 1
            1
        else
            Constant(cf*u/t)
        end
    catch err
        throw(ValidationError("Unable to convert [$t] to [$u]"))
    end
end

"Turn an expression into a Julia function w/ correct units behavior." # mostly for testing
function functionize(pt)
    syms = Symbolics.get_variables(pt)
    eval(build_function(constructunit(pt), syms, expression = Val{false}))
end

"Represent a constant as a Symbolic (esp. for lifting units to metadata level)."
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

"Convert symbolic expression `x` to have units `u` if possible."
function unitcoerce(u::Unitful.Unitlike, x::Symbolic)
    st =  _has_unit(x) ? x : constructunit(x)
    tu = _get_unit(st)
    output = unitfactor(u, tu) * st
    return SymbolicUtils.setmetadata(output,VariableUnit,u)
end

#Should run this at the end of @variables and @parameters
set_unitless(x::Vector) = [_has_unit(y) ? y : SymbolicUtils.setmetadata(y,VariableUnit,unitless) for y in x]

constructunit(x::Num) = constructunit(value(x))
function constructunit(x::Unitful.Quantity)
    return Constant(x.val,Dict(VariableUnit=>Unitful.unit(x)))
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

function constructunit(eq::ModelingToolkit.Equation)
    newterms = uniformize([eq.lhs,eq.rhs])
    return ~(newterms...)
    #return SymbolicUtils.setmetadata(output,VariableUnit,firstunit) #Fix this once Symbolics.jl Equations accept units
end

"Rewrite a set of equations by inserting appropriate unit conversion factors."
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
