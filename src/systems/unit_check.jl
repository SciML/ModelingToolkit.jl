#For dispatching get_unit
const Conditional = Union{typeof(ifelse)}
const Comparison = Union{typeof.([==, !=, ≠, <, <=, ≤, >, >=, ≥])...}

struct ValidationError <: Exception
    message::String
end

check_units(::Nothing, _...) = true

function __get_literal_unit(x)
    if x isa Pair
        x = x[1]
    end
    if !(x isa Union{Num, Symbolic})
        return nothing
    end
    v = value(x)
    u = getmetadata(v, VariableUnit, nothing)
    u isa DQ.AbstractQuantity ? screen_unit(u) : u
end
function __get_scalar_unit_type(v)
    u = __get_literal_unit(v)
    if u isa DQ.AbstractQuantity
        return Val(:DynamicQuantities)
    elseif u isa Unitful.Unitlike
        return Val(:Unitful)
    end
    return nothing
end
function __get_unit_type(vs′...)
    for vs in vs′
        if vs isa AbstractVector
            for v in vs
                u = __get_scalar_unit_type(v)
                u === nothing || return u
            end
        else
            v = vs
            u = __get_scalar_unit_type(v)
            u === nothing || return u
        end
    end
    return nothing
end

function screen_unit(result)
    if result isa DQ.AbstractQuantity
        d = DQ.dimension(result)
        if d isa DQ.Dimensions
            return result
        elseif d isa DQ.SymbolicDimensions
            return DQ.uexpand(oneunit(result))
        else
            throw(ValidationError("$result doesn't have a recognized unit"))
        end
    else
        throw(ValidationError("$result doesn't have any unit."))
    end
end

const unitless = DQ.Quantity(1.0)
get_literal_unit(x) = screen_unit(something(__get_literal_unit(x), unitless))

"""
Find the unit of a symbolic item.
"""
get_unit(x::Real) = unitless
get_unit(x::DQ.AbstractQuantity) = screen_unit(x)
get_unit(x::AbstractArray) = map(get_unit, x)
get_unit(x::Num) = get_unit(unwrap(x))
get_unit(x::Symbolics.Arr) = get_unit(unwrap(x))
get_unit(op::Differential, args) = get_unit(args[1]) / get_unit(op.x)
get_unit(op::Difference, args) = get_unit(args[1]) / get_unit(op.t)
get_unit(op::typeof(getindex), args) = get_unit(args[1])
get_unit(x::SciMLBase.NullParameters) = unitless
get_unit(op::typeof(instream), args) = get_unit(args[1])

function get_unit(op, args) # Fallback
    result = oneunit(op(get_unit.(args)...))
    try
        get_unit(result)
    catch
        throw(ValidationError("Unable to get unit for operation $op with arguments $args."))
    end
end

function get_unit(::Union{typeof(+), typeof(-)}, args)
    u = get_unit(args[1])
    if all(i -> get_unit(args[i]) == u, 2:length(args))
        return u
    end
end

function get_unit(op::Integral, args)
    unit = 1
    if op.domain.variables isa Vector
        for u in op.domain.variables
            unit *= get_unit(u)
        end
    else
        unit *= get_unit(op.domain.variables)
    end
    return get_unit(args[1]) * unit
end

equivalent(x, y) = isequal(x, y)
function get_unit(op::Conditional, args)
    terms = get_unit.(args)
    terms[1] == unitless ||
        throw(ValidationError(", in $op, [$(terms[1])] is not dimensionless."))
    equivalent(terms[2], terms[3]) ||
        throw(ValidationError(", in $op, units [$(terms[2])] and [$(terms[3])] do not match."))
    return terms[2]
end

function get_unit(op::typeof(Symbolics._mapreduce), args)
    if args[2] == +
        get_unit(args[3])
    else
        throw(ValidationError("Unsupported array operation $op"))
    end
end

function get_unit(op::Comparison, args)
    terms = get_unit.(args)
    equivalent(terms[1], terms[2]) ||
        throw(ValidationError(", in comparison $op, units [$(terms[1])] and [$(terms[2])] do not match."))
    return unitless
end

function get_unit(x::Symbolic)
    if (u = __get_literal_unit(x)) !== nothing
        screen_unit(u)
    elseif issym(x)
        get_literal_unit(x)
    elseif isadd(x)
        terms = get_unit.(arguments(x))
        firstunit = terms[1]
        for other in terms[2:end]
            termlist = join(map(repr, terms), ", ")
            equivalent(other, firstunit) ||
                throw(ValidationError(", in sum $x, units [$termlist] do not match."))
        end
        return firstunit
    elseif ispow(x)
        pargs = arguments(x)
        base, expon = get_unit.(pargs)
        @assert oneunit(expon) == unitless
        if base == unitless
            unitless
        else
            pargs[2] isa Number ? base^pargs[2] : (1 * base)^pargs[2]
        end
    elseif iscall(x)
        op = operation(x)
        if issym(op) || (iscall(op) && iscall(operation(op))) # Dependent variables, not function calls
            return screen_unit(getmetadata(x, VariableUnit, unitless)) # Like x(t) or x[i]
        elseif iscall(op) && !iscall(operation(op))
            gp = getmetadata(x, Symbolics.GetindexParent, nothing) # Like x[1](t)
            return screen_unit(getmetadata(gp, VariableUnit, unitless))
        end  # Actual function calls:
        args = arguments(x)
        return get_unit(op, args)
    else # This function should only be reached by Terms, for which `iscall` is true
        throw(ArgumentError("Unsupported value $x."))
    end
end

"""
Get unit of term, returning nothing & showing warning instead of throwing errors.
"""
function safe_get_unit(term, info)
    side = nothing
    try
        side = get_unit(term)
    catch err
        if err isa DQ.DimensionError
            @warn("$info: $(err.x) and $(err.y) are not dimensionally compatible.")
        elseif err isa ValidationError
            @warn(info*err.message)
        elseif err isa MethodError
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
    for (term, label) in zip(terms, labels)
        equnit = safe_get_unit(term, info * label)
        if equnit === nothing
            valid = false
        elseif !isequal(term, 0)
            if first_unit === nothing
                first_unit = equnit
                first_label = label
            elseif !equivalent(first_unit, equnit)
                valid = false
                str = "$info: units [$(first_unit)] for $(first_label) and [$(equnit)] for $(label) do not match."
                if oneunit(first_unit) == oneunit(equnit)
                    str *= " If there are non-SI units in the system, please use symbolic units like `us\"ms\"`"
                end
                @warn(str)
            end
        end
    end
    valid
end

function _validate(conn::Connection; info::String = "")
    valid = true
    syss = get_systems(conn)
    sys = first(syss)
    st = unknowns(sys)
    for i in 2:length(syss)
        s = syss[i]
        sst = unknowns(s)
        if length(st) != length(sst)
            valid = false
            @warn("$info: connected systems $(nameof(sys)) and $(nameof(s)) have $(length(st)) and $(length(sst)) unknowns, cannot connect.")
            continue
        end
        for (i, x) in enumerate(st)
            j = findfirst(isequal(x), sst)
            if j == nothing
                valid = false
                @warn("$info: connected systems $(nameof(sys)) and $(nameof(s)) do not have the same unknowns.")
            else
                aunit = safe_get_unit(x, info * string(nameof(sys)) * "#$i")
                bunit = safe_get_unit(sst[j], info * string(nameof(s)) * "#$j")
                if !equivalent(aunit, bunit)
                    valid = false
                    str = "$info: connected system unknowns $x ($aunit) and $(sst[j]) ($bunit) have mismatched units."
                    if oneunit(aunit) == oneunit(bunit)
                        str *= " If there are non-SI units in the system, please use symbolic units like `us\"ms\"`"
                    end
                    @warn(str)
                end
            end
        end
    end
    valid
end

function validate(jump::Union{VariableRateJump,
            ConstantRateJump}, t::Symbolic;
        info::String = "")
    newinfo = replace(info, "eq." => "jump")
    _validate([jump.rate, 1 / t], ["rate", "1/t"], info = newinfo) && # Assuming the rate is per time units
        validate(jump.affect!, info = newinfo)
end

function validate(jump::MassActionJump, t::Symbolic; info::String = "")
    left_symbols = [x[1] for x in jump.reactant_stoch] #vector of pairs of symbol,int -> vector symbols
    net_symbols = [x[1] for x in jump.net_stoch]
    all_symbols = vcat(left_symbols, net_symbols)
    allgood = _validate(all_symbols, string.(all_symbols); info)
    n = sum(x -> x[2], jump.reactant_stoch, init = 0)
    base_unitful = all_symbols[1] #all same, get first
    allgood && _validate([jump.scaled_rates, 1 / (t * base_unitful^n)],
        ["scaled_rates", "1/(t*reactants^$n))"]; info)
end

function validate(jumps::ArrayPartition{<:Union{Any, Vector{<:JumpType}}}, t::Symbolic)
    labels = ["in Mass Action Jumps,", "in Constant Rate Jumps,", "in Variable Rate Jumps,"]
    all([validate(jumps.x[idx], t, info = labels[idx]) for idx in 1:3])
end

function validate(eq::Equation; info::String = "")
    if typeof(eq.lhs) == Connection
        _validate(eq.rhs; info)
    else
        _validate([eq.lhs, eq.rhs], ["left", "right"]; info)
    end
end
function validate(eq::Equation,
        term::Union{Symbolic, DQ.AbstractQuantity, Num}; info::String = "")
    _validate([eq.lhs, eq.rhs, term], ["left", "right", "noise"]; info)
end
function validate(eq::Equation, terms::Vector; info::String = "")
    _validate(vcat([eq.lhs, eq.rhs], terms),
        vcat(["left", "right"], "noise  #" .* string.(1:length(terms))); info)
end

"""
Returns true iff units of equations are valid.
"""
function validate(eqs::Vector; info::String = "")
    all([validate(eqs[idx], info = info * " in eq. #$idx") for idx in 1:length(eqs)])
end
function validate(eqs::Vector, noise::Vector; info::String = "")
    all([validate(eqs[idx], noise[idx], info = info * " in eq. #$idx")
         for idx in 1:length(eqs)])
end
function validate(eqs::Vector, noise::Matrix; info::String = "")
    all([validate(eqs[idx], noise[idx, :], info = info * " in eq. #$idx")
         for idx in 1:length(eqs)])
end
function validate(eqs::Vector, term::Symbolic; info::String = "")
    all([validate(eqs[idx], term, info = info * " in eq. #$idx") for idx in 1:length(eqs)])
end
validate(term::Symbolics.SymbolicUtils.Symbolic) = safe_get_unit(term, "") !== nothing

"""
Throws error if units of equations are invalid.
"""
function check_units(::Val{:DynamicQuantities}, eqs...)
    validate(eqs...) ||
        throw(ValidationError("Some equations had invalid units. See warnings for details."))
end
