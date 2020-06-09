const IndexMap = Dict{Char,Char}(
            '0' => '₀',
            '1' => '₁',
            '2' => '₂',
            '3' => '₃',
            '4' => '₄',
            '5' => '₅',
            '6' => '₆',
            '7' => '₇',
            '8' => '₈',
            '9' => '₉')
function map_subscripts(indices)
    str = string(indices)
    join(IndexMap[c] for c in str)
end

"""
$(TYPEDEF)

A named variable which represents a numerical value. The variable is uniquely
identified by its `name`, and all variables with the same `name` are treated
as equal.

# Fields
$(FIELDS)

For example, the following code defines an independent variable `t`, a parameter
`α`, a function parameter `σ`, a variable `x`, which depends on `t`, a variable
`y` with no dependents, a variable `z`, which depends on `t`, `α`, and `x(t)`
and parameters `β₁` and `β₂`.


```julia
t = Variable(:t)()  # independent variables are treated as known
α = Variable(:α)()  # parameters are known
σ = Variable(:σ)    # left uncalled, since it is used as a function
w = Variable(:w)   # unknown, left uncalled
x = Variable(:x)(t)  # unknown, depends on `t`
y = Variable(:y)()   # unknown, no dependents
z = Variable(:z)(t, α, x)  # unknown, multiple arguments
β₁ = Variable(:β, 1)() # with index 1
β₂ = Variable(:β, 2)() # with index 2

expr = β₁ * x + y^α + σ(3) * (z - t) - β₂ * w(t - 1)
```
"""
struct Variable{T} <: Function
    """The variable's unique name."""
    name::Symbol
    Variable(name) = new{Real}(name)
    Variable{T}(name) where T = new{T}(name)
    function Variable{T}(name, indices...) where T
        var_name = Symbol("$(name)$(join(map_subscripts.(indices), "ˏ"))")
        Variable{T}(var_name)
    end
end

function Variable(name, indices...)
    var_name = Symbol("$(name)$(join(map_subscripts.(indices), "ˏ"))")
    Variable(var_name)
end

vartype(::Variable{T}) where T = T
(x::Variable)(args...) = Operation(x, collect(Expression, args))
rename(x::Variable{T},name) where T = Variable{T}(name)

Base.isequal(x::Variable, y::Variable) = x.name == y.name
Base.print(io::IO, x::Variable) = show(io, x)
Base.show(io::IO, x::Variable) = print(io, x.name)
function Base.show(io::IO, ::MIME"text/plain", x::Variable)
    print(io, x.name)
end


"""
$(TYPEDEF)

An expression which wraps a constant numerical value.
"""
struct Constant <: Expression
    """The constant value"""
    value
end

"""
$(TYPEDSIGNATURES)

Get the value of a [`ModelingToolkit.Constant`](@ref).
"""
Base.get(c::Constant) = c.value

Base.iszero(c::Constant) = iszero(c.value)

Base.isone(ex::Expression)  = isa(ex, Constant) && isone(ex.value)

# Variables use isequal for equality since == is an Operation
Base.isequal(c::Constant, n::Number) = c.value == n
Base.isequal(n::Number, c::Constant) = c.value == n
Base.isequal(a::Constant, b::Constant) = a.value == b.value

Base.convert(::Type{Expr}, c::Constant) = c.value


# Build variables more easily
function _parse_vars(macroname, type, x)
    ex = Expr(:block)
    var_names = Symbol[]
    # if parsing things in the form of
    # begin
    #     x
    #     y
    #     z
    # end
    x = x isa Tuple && first(x) isa Expr && first(x).head == :tuple ? first(x).args : x # tuple handling
    x = flatten_expr!(x)
    for _var in x
        iscall = isa(_var, Expr) && _var.head == :call
        isarray = isa(_var, Expr) && _var.head == :ref
        issym  = _var isa Symbol
        @assert iscall || isarray || issym "@$macroname expects a tuple of expressions or an expression of a tuple (`@$macroname x y z(t) v[1:3] w[1:2,1:4]` or `@$macroname x, y, z(t) v[1:3] w[1:2,1:4]`)"

        if iscall
            var_name, expr = _construct_vars(_var.args[1], type, _var.args[2:end])
        else
            var_name, expr = _construct_vars(_var, type, nothing)
        end
        push!(var_names, var_name)
        push!(ex.args, expr)
    end
    push!(ex.args, build_expr(:tuple, var_names))
    return ex
end

function _construct_vars(_var, type, call_args)
    issym  = _var isa Symbol
    isarray = isa(_var, Expr) && _var.head == :ref
    if isarray
        var_name = _var.args[1]
        indices = _var.args[2:end]
        expr = _construct_array_vars(var_name, type, call_args, indices...)
    else
        # Implicit 0-args call
        var_name = _var
        expr = _construct_var(var_name, type, call_args)
    end
    var_name, :($var_name = $expr)
end

function _construct_var(var_name, type, call_args)
    if call_args === nothing
        :(Variable{$type}($(Meta.quot(var_name)))())
    elseif !isempty(call_args) && call_args[end] == :..
        :(Variable{$type}($(Meta.quot(var_name))))
    else
        :(Variable{$type}($(Meta.quot(var_name)))($(call_args...)))
    end
end

function _construct_var(var_name, type, call_args, ind)
    if call_args === nothing
        :(Variable{$type}($(Meta.quot(var_name)), $ind...)())
    elseif !isempty(call_args) && call_args[end] == :..
        :(Variable{$type}($(Meta.quot(var_name)), $ind...))
    else
        :(Variable{$type}($(Meta.quot(var_name)), $ind...)($(call_args...)))
    end
end


function _construct_array_vars(var_name, type, call_args, indices...)
    :(map(Iterators.product($(indices...))) do ind
        $(_construct_var(var_name, type, call_args, :ind))
    end)
end


"""
$(SIGNATURES)

Define one or more unknown variables.

```julia
@parameters t α σ(..) β[1:2]
@variables w(..) x(t) y() z(t, α, x)

expr = β₁* x + y^α + σ(3) * (z - t) - β₂ * w(t - 1)
```

Note that `@parameters` and `@variables` implicitly add `()` to values that
are not given a call. The former specifies the values as known, while the
latter specifies it as unknown. `(..)` signifies that the value should be
left uncalled.

Sometimes it is convenient to define arrays of variables to model things like `x₁,…,x₃`.
The `@variables` and `@parameters` macros support this with the following syntax:

```julia
@variables x[1:3];
x

3-element Array{Operation,1}:
 x₁()
 x₂()
 x₃()

# support for arbitrary ranges and tensors
@variables y[2:3,1:5:6];
y

2×2 Array{Operation,2}:
    y₂̒₁() y₂̒₆()
    y₃̒₁() y₃̒₆()

# also works for dependent variables
@parameters t; @variables z[1:3](t);
z

3-element Array{Operation,1}:
 z₁(t())
 z₂(t())
 z₃(t())
```
"""
macro variables(xs...)
    esc(_parse_vars(:variables, Real, xs))
end

"""
$(TYPEDSIGNATURES)

Renames the variable `x` to have `name`.
"""
function rename(x::Variable,name::Symbol)
    Variable{vartype(x)}(name)
end

TreeViews.hastreeview(x::Variable) = true
function TreeViews.treelabel(io::IO,x::Variable,
                             mime::MIME"text/plain" = MIME"text/plain"())
  show(io,mime,Text(x.name))
end

"""
varmap_to_vars(varmap,varlist)

Takes a list of pairs of variables=>values and an ordered list of variables and
creates the array of values in the correct order
"""
function varmap_to_vars(varmap,varlist)
    out = similar(varmap,typeof(last(first(varmap))))
    for i in 1:length(varmap)
        ivar = convert(Variable,varmap[i][1])
        j = findfirst(x->ivar.name == convert(Variable,x).name,varlist)
        out[j] = varmap[i][2]
    end

    # Make output match varmap in type and shape
    # Does things like MArray->SArray
    ArrayInterface.restructure(varmap,out)
end

varmap_to_vars(varmap::DiffEqBase.NullParameters,varlist) = varmap
varmap_to_vars(varmap::Nothing,varlist) = varmap
