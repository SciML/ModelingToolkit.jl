Base.:*(x::Union{Num,Symbolic},y::Unitful.AbstractQuantity) = x * y
Base.:/(x::Union{Num,Symbolic},y::Unitful.AbstractQuantity) = x / y

struct ValidationError <: Exception
    message::String
end

function screen_unit(result)
    result isa Unitful.Unitlike || throw(ValidationError("Unit must be a subtype of Unitful.Unitlike, not $(typeof(result))."))
    result isa Unitful.ScalarUnits || throw(ValidationError("Non-scalar units such as $result are not supported. Use a scalar unit instead."))
    result == u"°" && throw(ValidationError("Degrees are not supported. Use radians instead."))
    result
end
equivalent(x,y) = isequal(1*x,1*y)
unitless = Unitful.unit(1)

#For dispatching get_unit
Literal = Union{Sym,Symbolics.ArrayOp,Symbolics.Arr,Symbolics.CallWithMetadata}
Conditional = Union{typeof(ifelse),typeof(IfElse.ifelse)}
Comparison = Union{typeof(Base.:>), typeof(Base.:<), typeof(==)}

"Find the unit of a symbolic item."
get_unit(x::Real) = unitless
get_unit(x::Unitful.Quantity) = screen_unit(Unitful.unit(x))
get_unit(x::AbstractArray) = map(get_unit,x)
get_unit(x::Num) = get_unit(value(x))
get_unit(x::Literal) = screen_unit(getmetadata(x,VariableUnit, unitless))
get_unit(op::Differential, args) = get_unit(args[1]) / get_unit(op.x)
get_unit(op::Difference, args) =   get_unit(args[1]) / get_unit(op.t) #why are these not identical?!?
get_unit(op::typeof(getindex),args) = get_unit(args[1]) 
function get_unit(op,args) #Fallback
    result = op(1 .* get_unit.(args)...)
    try 
        unit(result)
    catch 
        throw(ValidationError("Unable to get unit for operation $op with arguments $args."))
    end
end

function get_unit(op::Integral,args)
    unit = 1
    if op.x isa Vector
        for u in op.x
            unit *= get_unit(u)
        end
    else
        unit *= get_unit(op.x)
    end
    return get_unit(args[1]) * unit
end

function get_unit(x::Pow)
    pargs = arguments(x)
    base,expon = get_unit.(pargs)
    @assert expon isa Unitful.DimensionlessUnits
    if base == unitless
        unitless
    else
        pargs[2] isa Number ? base^pargs[2] : (1*base)^pargs[2]
    end
end

function get_unit(x::Add)
    terms = get_unit.(arguments(x))
    firstunit = terms[1]
    for other in terms[2:end]
        termlist = join(map(repr, terms), ", ")
        equivalent(other, firstunit) || throw(ValidationError(", in sum $x, units [$termlist] do not match."))
    end
    return firstunit
end

function get_unit(op::Conditional, args)
    terms = get_unit.(args)
    terms[1] == unitless || throw(ValidationError(", in $x, [$(terms[1])] is not dimensionless."))
    equivalent(terms[2], terms[3]) || throw(ValidationError(", in $x, units [$(terms[2])] and [$(terms[3])] do not match."))
    return terms[2]
end

function get_unit(op::typeof(Symbolics._mapreduce),args)
    if args[2] == +
        get_unit(args[3])
    else
        throw(ValidationError("Unsupported array operation $op"))
    end
end

function get_unit(op::Comparison, args)
    terms = get_unit.(args)
    equivalent(terms[1], terms[2]) || throw(ValidationError(", in comparison $x, units [$(terms[1])] and [$(terms[2])] do not match."))
    return unitless
end

function get_unit(x::Symbolic) 
    if SymbolicUtils.istree(x)
        op = operation(x)
        if op isa Sym # Not a real function call, just a dependent variable. Unit is on the Sym.
            return screen_unit(getmetadata(x, VariableUnit, unitless)) 
        elseif op isa Term && !(operation(op) isa Term) #
            gp = getmetadata(x,Symbolics.GetindexParent,nothing)
            return screen_unit(getmetadata(gp, VariableUnit, unitless))
        elseif op isa Term
            return screen_unit(getmetadata(x, VariableUnit, unitless))
        end        
        args = arguments(x)
        return get_unit(op, args)
    else #This function should only be reached by Terms, for which `istree` is true
        throw(ArgumentError("Unsupported value $x."))
    end
end

"Get unit of term, returning nothing & showing warning instead of throwing errors."
function safe_get_unit(term, info)
    side = nothing
    try
        side = get_unit(term)
    catch err
        if err isa Unitful.DimensionError
            @warn("$info: $(err.x) and $(err.y) are not dimensionally compatible.")
        elseif err isa ValidationError
            @warn(info*err.message)
        elseif err isa MethodError #Warning: Unable to get unit for operation x[1] with arguments SymbolicUtils.Sym{Real, Base.ImmutableDict{DataType, Any}}[t].
            @warn("$info: no method matching $(err.f) for arguments $(typeof.(err.args)).")
        else
            rethrow()
        end
    end
    side
end

function _validate(terms::Vector, labels::Vector{String}; info::String = "")
    valid = true
    first_unit = nothing
    first_label = nothing
    for (term,label) in zip(terms,labels)
        equnit = safe_get_unit(term, info*label)
        if equnit === nothing
            valid = false
        elseif !isequal(term,0)
            if first_unit === nothing
                first_unit = equnit
                first_label = label
            elseif !equivalent(first_unit, equnit)
                valid = false
                @warn("$info: units [$(first_unit)] for $(first_label) and [$(equnit)] for $(label) do not match.")
            end
        end
    end
    valid
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
validate(eqs::Vector; info::String = "") = all([validate(eqs[idx], info = info*" in eq. #$idx") for idx in 1:length(eqs)])
validate(eqs::Vector, noise::Vector; info::String = "") = all([validate(eqs[idx], noise[idx], info = info*" in eq. #$idx") for idx in 1:length(eqs)])
validate(eqs::Vector, noise::Matrix; info::String = "") = all([validate(eqs[idx], noise[idx, :], info = info*" in eq. #$idx") for idx in 1:length(eqs)])
validate(eqs::Vector, term::Symbolic; info::String = "") = all([validate(eqs[idx], term, info = info*" in eq. #$idx") for idx in 1:length(eqs)])
validate(term::Symbolics.SymbolicUtils.Symbolic) = safe_get_unit(term,"") !== nothing

"Throws error if units of equations are invalid."
check_units(eqs...) = validate(eqs...) || throw(ValidationError("Some equations had invalid units. See warnings for details."))
all_dimensionless(states) = all(map(x->safe_get_unit(x,"") in (unitless,nothing),states))