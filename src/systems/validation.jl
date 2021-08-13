Base.:*(x::Union{Num,Symbolic},y::Unitful.AbstractQuantity) = x * y

"Find the units of a symbolic item."
get_units(x::Real) = 1
get_units(x::Unitful.Quantity) = 1 * Unitful.unit(x)
get_units(x::Num) = get_units(value(x))
function get_units(x::Symbolic)
    if x isa Sym || operation(x) isa Sym || (operation(x) isa Term && operation(x).f == getindex) || x isa Symbolics.ArrayOp
        if x.metadata !== nothing
            symunits = get(x.metadata, VariableUnit, 1)
        else
            symunits = 1
        end
        return oneunit(1 * symunits)
    elseif operation(x) isa Differential || operation(x) isa Difference
        return get_units(arguments(x)[1]) / get_units(arguments(arguments(x)[1])[1])
    elseif x isa Pow
        pargs = arguments(x)
        base,expon = get_units.(pargs)
        uconvert(NoUnits, expon) # This acts as an assertion
        return base == 1 ? 1 : operation(x)(base, pargs[2])
    elseif x isa Add # Cannot simply add the units b/c they may differ in magnitude (eg, kg vs g)
        terms = get_units.(arguments(x))
        firstunit = 1*unit(terms[1])
        for other in terms[2:end]
            isequal(1 * unit(other), firstunit) || throw(Unitful.DimensionError(x, terms))
        end
        return 1 * firstunit
    elseif operation(x) == Symbolics._mapreduce 
        if x.arguments[2] == +
            get_units(x.arguments[3])
        else
            throw(ArgumentError("Unknown array operation $x"))
        end
    else
        return oneunit(operation(x)(get_units.(arguments(x))...))
    end
end

"Get units of term, returning nothing & showing warning instead of throwing errors."
function safe_get_units(term, info)
    side = nothing
    try
        side = get_units(term)
    catch err
        if err isa Unitful.DimensionError 
            @warn("$info: $(err.x) and $(err.y) are not dimensionally compatible.")
        elseif err isa MethodError
            @warn("$info: no method matching $(err.f) for arguments $(typeof.(err.args)).")
        else
            rethrow()
        end
    end
    side
end

function _validate(terms::Vector, labels::Vector{String}; info::String = "")
    equnits = safe_get_units.(terms, info*", ".*labels)
    allthere = all(map(x -> x!==nothing, equnits))
    allmatching = true
    first_unit = nothing
    if allthere
        for idx in 1:length(equnits)
            if !isequal(terms[idx],0)
                if first_unit === nothing
                    first_unit = equnits[idx]
                elseif !isequal(first_unit, equnits[idx])
                    allmatching = false
                    @warn("$info: units $(equnits[1]) for $(labels[1]) and $(equnits[idx]) for $(labels[idx]) do not match.")
           
                end
            end
        end
    end
    allthere && allmatching
end

function validate(jump::Union{ModelingToolkit.VariableRateJump, ModelingToolkit.ConstantRateJump}, t::Symbolic; info::String = "")
    _validate([jump.rate, 1/t], ["rate", "1/t"], info = info) && # Assuming the rate is per time units
    validate(jump.affect!,info = info) 
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
    labels = ["in Mass Action Jumps, ", "in Constant Rate Jumps, ", "in Variable Rate Jumps, "]
    all([validate(jumps.x[idx], t, info = labels[idx]) for idx in 1:3])
end

validate(eq::ModelingToolkit.Reaction; info::String = "") = _validate([oderatelaw(eq)], ["",], info = info)
validate(eq::ModelingToolkit.Equation; info::String = "") = _validate([eq.lhs, eq.rhs], ["left", "right"], info = info)
validate(eq::ModelingToolkit.Equation, term::Union{Symbolic,Unitful.Quantity,Num}; info::String = "") = _validate([eq.lhs, eq.rhs, term], ["left","right","noise"], info = info)
validate(eq::ModelingToolkit.Equation, terms::Vector; info::String = "") = _validate(vcat([eq.lhs, eq.rhs], terms), vcat(["left", "right"], "noise #".*string.(1:length(terms))), info = info)

"Returns true iff units of equations are valid."
validate(eqs::Vector; info::String = "") = all([validate(eqs[idx], info = info*"in eq. #$idx") for idx in 1:length(eqs)])
validate(eqs::Vector, noise::Vector; info::String = "") = all([validate(eqs[idx], noise[idx], info = info*"in eq. #$idx") for idx in 1:length(eqs)])
validate(eqs::Vector, noise::Matrix; info::String = "") = all([validate(eqs[idx], noise[idx, :], info = info*"in eq. #$idx") for idx in 1:length(eqs)])
validate(eqs::Vector, term::Symbolic; info::String = "") = all([validate(eqs[idx], term, info = info*"in eq. #$idx") for idx in 1:length(eqs)])

"Throws error if units of equations are invalid."
check_units(eqs...) = validate(eqs...) || throw(ArgumentError("Some equations had invalid units. See warnings for details."))
