Base.:*(x::Union{Num,Symbolic},y::Unitful.AbstractQuantity) = x * y


function vartype(x::Symbolic)
    if !(x.metadata isa Nothing)
        return haskey(x.metadata, VariableUnit) ? x.metadata[VariableUnit] : 1
    end
    1
end
vartype(x::Num) = vartype(value(x))

instantiate(x) = 1
instantiate(x::Num) = instantiate(value(x))
function instantiate(x::Symbolic)
    vx = value(x)
    if vx isa Sym || operation(vx) isa Sym
    elseif operation(vx) isa Differential
        return oneunit(1 * vartype(vx))
        return instantiate(arguments(vx)[1]) / instantiate(arguments(arguments(vx)[1])[1])
    elseif vx isa Pow
        pargs = arguments(vx)
        base,expon = instantiate.(pargs)
        uconvert(NoUnits, expon) # This acts as an assertion
        return base == 1 ? 1 : operation(vx)(base, pargs[2])
    elseif vx isa Add # Cannot simply add the units b/c they may differ in magnitude (eg, kg vs g)
        terms = instantiate.(arguments(vx))
        firstunit = unit(terms[1])
        @assert all(map(x -> ustrip(firstunit, x) == 1, terms[2:end]))
        return 1 * firstunit
    else
        return oneunit(operation(vx)(instantiate.(arguments(vx))...))
    end
end

function validate(eq::ModelingToolkit.Equation; eqnum = 1)
    lhs = rhs = nothing
    try
        lhs = instantiate(eq.lhs)
    catch err
        if err isa Unitful.DimensionError 
            @warn("In left-hand side of eq. #$eqnum: $(eq.lhs),  $(err.x) and $(err.y) are not dimensionally compatible.")
        elseif err isa MethodError
            @warn("In left-hand side of eq. #$eqnum: $(err.f) doesn't accept $(err.args).")
        else
            rethrow()
        end
    end
    try
        rhs = instantiate(eq.rhs)
    catch err
        if err isa Unitful.DimensionError
            @warn("In right-hand side of eq. #$eqnum: $(eq.rhs), $(err.x) and $(err.y) are not dimensionally compatible.")
        elseif err isa MethodError
            @warn("In right-hand side of eq. #$eqnum: $(err.f) doesn't accept $(err.args).")
        else
            rethrow()
        end
    end
    if (rhs !== nothing) && (lhs !== nothing)
        if !isequal(lhs, rhs)
            @warn("In eq. #$eqnum, left-side units ($lhs) and right-side units ($rhs) don't match.")
        end
    end
    (rhs !== nothing) && (lhs !== nothing) && isequal(lhs, rhs)
end

function validate(eqs::Vector{ModelingToolkit.Equation})
    correct = [validate(eqs[idx],eqnum=idx) for idx in 1:length(eqs)]
    all(correct) || throw(ArgumentError("Invalid equations, see warnings for details."))
end

validate(sys::AbstractODESystem) = validate(equations(sys))
