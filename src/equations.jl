export Equation


mutable struct Equation
    lhs::Expression
    rhs::Expression
end
Base.broadcastable(eq::Equation) = Ref(eq)

Base.:~(lhs::Expression, rhs::Expression) = Equation(lhs, rhs)
Base.:~(lhs::Expression, rhs::Number    ) = Equation(lhs, rhs)
Base.:~(lhs::Number    , rhs::Expression) = Equation(lhs, rhs)


function extract_elements(eqs, targetmap, default = nothing)
    elems = Dict{Symbol, Vector{Variable}}()
    names = Dict{Symbol, Set{Symbol}}()
    if default == nothing
        targets =  unique(collect(values(targetmap)))
    else
        targets = [unique(collect(values(targetmap))), default]
    end
    for target in targets
        elems[target] = Vector{Variable}()
        names[target] = Set{Symbol}()
    end
    for eq in eqs
        extract_elements!(eq, elems, names, targetmap, default)
    end
    Tuple(elems[target] for target in targets)
end
# Walk the tree recursively and push variables into the right set
function extract_elements!(op, elems, names, targetmap, default)
    args = isa(op, Equation) ? Expression[op.lhs, op.rhs] : op.args

    for arg in args
        if arg isa Operation
            extract_elements!(arg, elems, names, targetmap, default)
        elseif arg isa Variable
            if default == nothing
                target = haskey(targetmap, arg.subtype) ? targetmap[arg.subtype] : continue
            else
                target = haskey(targetmap, arg.subtype) ? targetmap[arg.subtype] : default
            end
            if !in(arg.name, names[target])
                push!(names[target], arg.name)
                push!(elems[target], arg)
            end
        end
    end
end
