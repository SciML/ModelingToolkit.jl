Base.:*(x::Union{Num,Symbolic},y::Unitful.AbstractQuantity) = x * y

"Find the units of a symbolic item."
get_units(x) = 1
get_units(x::Unitful.Quantity) = 1 * Unitful.unit(x)
get_units(x::Num) = get_units(value(x))
function get_units(x::Symbolic)
    vx = value(x)
    if vx isa Sym || operation(vx) isa Sym || (operation(vx) isa Term && operation(vx).f == getindex) || vx isa Symbolics.ArrayOp
        if x.metadata !== nothing
            symunits = haskey(x.metadata, VariableUnit) ? x.metadata[VariableUnit] : 1
        else
            symunits = 1
        end
        return oneunit(1 * symunits)
    elseif operation(vx) isa Differential || operation(vx) isa Difference
        return get_units(arguments(vx)[1]) / get_units(arguments(arguments(vx)[1])[1])
    elseif vx isa Pow
        pargs = arguments(vx)
        base,expon = get_units.(pargs)
        uconvert(NoUnits, expon) # This acts as an assertion
        return base == 1 ? 1 : operation(vx)(base, pargs[2])
    elseif vx isa Add # Cannot simply add the units b/c they may differ in magnitude (eg, kg vs g)
        terms = get_units.(arguments(vx))
        firstunit = unit(terms[1])
        @assert all(map(x -> ustrip(firstunit, x) == 1, terms[2:end]))
        return 1 * firstunit
    elseif operation(vx) == Symbolics._mapreduce 
        if vx.arguments[2] == +
            get_units(vx.arguments[3])
        else
            throw(ArgumentError("Unknown array operation $vx"))
        end
    else
        return oneunit(operation(vx)(get_units.(arguments(vx))...))
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
        elseif err isa MethodError #TODO: filter for only instances where the arguments are unitful
            @warn("$info: no method matching $(err.f) for arguments $(err.args).")
        else
            rethrow()
        end
    end
    side
end

function _validate(terms::Vector, labels::Vector; info::String = "")
    equnits = safe_get_units.(terms, info*", ".*labels)
    allthere = all(map(x -> x!==nothing, equnits))
    allmatching = true
    if allthere
        for idx in 2:length(equnits)
            if !isequal(equnits[1], equnits[idx])
                allmatching = false
                @warn("$info: units $(equnits[1]) for $(labels[1]) and $(equnits[idx]) for $(labels[idx]) do not match.")
            end
        end
    end
    allthere && allmatching
end

function validate(eq::ModelingToolkit.Equation; info::String = "")
    labels = ["left-hand side", "right-hand side"]
    terms = [eq.lhs, eq.rhs]
    _validate(terms, labels, info = info)
end

function validate(eq::ModelingToolkit.Equation, noiseterm; info::String = "")
    labels = ["left-hand side", "right-hand side", "noise term"]
    terms = [eq.lhs, eq.rhs, noiseterm]
    _validate(terms, labels, info = info)
end

function validate(eq::ModelingToolkit.Equation, noisevec::Vector; info::String = "")
    labels = vcat(["left-hand side", "right-hand side"], "noise term #".* string.(1:length(noisevec)))
    terms = vcat([eq.lhs, eq.rhs], noisevec)
    _validate(terms, labels, info = info)
end

function validate(eqs::Vector{ModelingToolkit.Equation})
    all([validate(eqs[idx], info = "In eq. #$idx") for idx in 1:length(eqs)])
end

function validate(eqs::Vector{ModelingToolkit.Equation}, noise::Vector)
    all([validate(eqs[idx], noise[idx], info = "In eq. #$idx") for idx in 1:length(eqs)])
end

function validate(eqs::Vector{ModelingToolkit.Equation}, noise::Matrix)
    all([validate(eqs[idx], noise[idx, :], info = "In eq. #$idx") for idx in 1:length(eqs)])
end

"Returns true iff units of equations are valid."
validate(eqs::Vector) = validate(convert.(ModelingToolkit.Equation, eqs))

"Throws error if units of equations are invalid."
check_units(eqs...) = validate(eqs...) || throw(ArgumentError("Some equations had invalid units. See warnings for details."))
