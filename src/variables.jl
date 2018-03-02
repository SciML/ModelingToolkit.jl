# <: Real to make tracing easier. Maybe a bad idea?
struct Variable <: Real
    name::Symbol
    subtype::Symbol
    value
    value_type::DataType
    diff::Union{AbstractOperator,Void}
end

Variable(name,subtype,value = nothing,value_type = typeof(value)) =
                                 Variable(name,subtype,value,value_type,nothing)
Parameter(name,args...) = Variable(name,:Parameter,args...)
Constant(value) = Variable(:None,:Constant,value,typeof(value))
DependentVariable(name,args...) = Variable(name,:DependentVariable,args...)
IndependentVariable(name,args...) = Variable(name,:IndependentVariable,args...)

export Variable,Parameter,Constant,DependentVariable,IndependentVariable

# Variables use isequal for equality since == is an Operation
function Base.isequal(x::Variable,y::Variable)
    x.name == y.name && x.subtype == y.subtype && x.value == y.value &&
    x.value_type == y.value_type && x.diff == y.diff
end

function Base.Expr(x::Variable)
    if x.diff == nothing
        return :($(x.name))
    else
        return :($(Symbol("$(x.name)_$(x.diff.x.name)")))
    end
end

function Base.show(io::IO, A::Variable)
    if A.subtype == :Constant
        print(io,"Constant($(A.value))")
    else
        str = "$(A.subtype)($(A.name))"
        if A.value != nothing
            str *= ", value = " * string(A.value)
        end

        if A.diff != nothing
            str *= ", diff = " * string(A.diff)
        end

        print(io,str)
    end
end


function extract_elements(ops, eltypes)
    elems = Dict{Symbol, Vector{Variable}}()
    names = Dict{Symbol, Set{Symbol}}()
    for el in eltypes
        elems[el] = Vector{Variable}()
        names[el] = Set{Symbol}()
    end
    for op in ops
        extract_elements!(op, elems, names)
    end
    Tuple(elems[el] for el in eltypes)
end
# Walk the tree recursively and push variables into the right set
function extract_elements!(op::Operation, elems, names)
    for arg in op.args
        if arg isa Operation
            extract_elements!(arg, elems, names)
        elseif arg isa Variable && haskey(elems, arg.subtype) && !in(arg.name, names[arg.subtype])
            arg.subtype == :DependentVariable && arg.diff != nothing && continue
            push!(names[arg.subtype], arg.name)
            push!(elems[arg.subtype], arg)
        end
    end
end
