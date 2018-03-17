# <: Real to make tracing easier. Maybe a bad idea?
struct Variable <: Expression
    name::Symbol
    value
    value_type::DataType
    subtype::Symbol
    diff::Union{AbstractOperator,Void}
    dependents::Vector{Variable}
    description::String
    flow::Bool
    domain
    size
    context
end

Variable(name,
         value = nothing,
         value_type = typeof(value);
         subtype::Symbol=:Variable,
         dependents::Vector{Variable} = Variable[],
         flow::Bool = false,
         description::String = "",
         domain = nothing,
         size = nothing,
         context = nothing) =
         Variable(name,value,value_type,subtype,nothing,
                  dependents,description,flow,domain,size,context)
Variable(name,args...;kwargs...) = Variable(name,args...;subtype=:Variable,kwargs...)
Parameter(name,args...;kwargs...) = Variable(name,args...;subtype=:Parameter,kwargs...)
Constant(value::Number) = Variable(Symbol(value),value,typeof(value);subtype=:Constant)
Constant(name,args...;kwargs...) = Variable(name,args...;subtype=:Constant,kwargs...)
DependentVariable(name,args...;kwargs...) = Variable(name,args...;subtype=:DependentVariable,kwargs...)
StateVariable(name,args...;kwargs...) = Variable(name,args...;subtype=:StateVariable,kwargs...)
ControlVariable(name,args...;kwargs...) = Variable(name,args...;subtype=:ControlVariable,kwargs...)
IndependentVariable(name,args...;kwargs...) = Variable(name,args...;subtype=:IndependentVariable,kwargs...)
JumpVariable(name,args...;kwargs...) = Variable(name,args...;subtype=:JumpVariable,kwargs...)
NoiseVariable(name,args...;kwargs...) = Variable(name,args...;subtype=:NoiseVariable,kwargs...)

export Variable,Parameter,Constant,DependentVariable,IndependentVariable,JumpVariable,NoiseVariable,
       @Var, @DVar, @IVar, @Param, @Const

# Variables use isequal for equality since == is an Operation
function Base.:(==)(x::Variable,y::Variable)
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

extract_idv(eq) = eq.args[1].diff.x

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
function extract_elements!(op::AbstractOperation, elems, names)
    for arg in op.args
        if arg isa Operation
            extract_elements!(arg, elems, names)
        elseif arg isa Variable && haskey(elems, arg.subtype) && !in(arg.name, names[arg.subtype])
            push!(names[arg.subtype], arg.name)
            push!(elems[arg.subtype], arg)
        end
    end
end

# Build variables more easily
function _parse_vars(macroname, fun, x)
    ex = Expr(:block)
    lhss = Symbol[]
    # if parsing things in the form of
    # begin
    #     x
    #     y
    #     z = exp(2)
    # end
    x = flatten_expr!(x)
    for _var in x
        iscall = typeof(_var) <: Expr && _var.head == :call
        @show iscall
        if iscall
            dependents = Variable[:($(esc(d))) for d in _var.args[2:end]]
            var = _var.args[1]
        else
            dependents = Variable[]
            var = _var
        end
        issym    = var isa Symbol
        isassign = issym ? false : var.head == :(=)
        @assert issym || isassign "@$macroname expects a tuple of expressions!\nE.g. `@$macroname x y z=1`"
        if issym
            lhs = var
            push!(lhss, lhs)
            expr = :( $lhs = $fun( Symbol($(String(lhs))) ,
                      dependents = $dependents))
        end
        if isassign
            lhs = var.args[1]
            rhs = var.args[2]
            push!(lhss, lhs)
            expr = :( $lhs = $fun( Symbol($(String(lhs))) , $rhs))
        end
        push!(ex.args, expr)
    end
    push!(ex.args, Expr(:tuple, lhss...))
    ex
end

for funs in ((:DVar, :DependentVariable), (:IVar, :IndependentVariable),
             (:Var, :Variable),
             (:Param, :Parameter))
    @eval begin
        macro ($(funs[1]))(x...)
            esc(_parse_vars(String($funs[1]), $funs[2], x))
        end
    end
end

function _const_assign(x)
    ex = Expr(:block)
    lhss = Symbol[]
    x = flatten_expr!(x)
    for eq in x
        @assert eq isa Expr && eq.head == :(=) "@Const expects a tuple of assignments!\nE.g. `@Const D=t W=g`"
        lhs = eq.args[1]
        push!(lhss, lhs)
        rhs = eq.args[2]
        expr = :($lhs = Constant($rhs))
        push!(ex.args,  expr)
    end
    push!(ex.args, Expr(:tuple, lhss...))
    ex
end

macro Const(x...)
    esc(_const_assign(x))
end
