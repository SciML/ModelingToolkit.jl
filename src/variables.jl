using SymbolicUtils: FnType

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

function subscripted(T::Type, name::Symbol, indices...)
    var_name = Symbol("$(name)$(join(map_subscripts.(indices), "ˏ"))")
    Sym{T}(var_name)
end

function subscripted(name::Symbol, indices...)
    subscripted(Real, name, indices...)
end

rename(x::Sym{T},name) where T = Sym{T}(name)
rename(x::Term, name) where T = x.op isa Sym ? rename(x.op, name)(x.args...) : error("can't rename $x to $name")

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
    rhs = build_expr(:tuple, var_names)
    push!(ex.args, :(($(var_names...),) = $rhs))
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
    expr = if call_args === nothing
        :($Num($Sym{$type}($(Meta.quot(var_name)))))
    elseif !isempty(call_args) && call_args[end] == :..
        :($Num($Sym{$FnType{Tuple, $type}}($(Meta.quot(var_name))))) # XXX: using Num as output
    else
        :($Num($Sym{$FnType{NTuple{$(length(call_args)), Any}, $type}}($(Meta.quot(var_name)))($(map(x->:($value($x)), call_args)...))))
    end
end

function _construct_var(var_name, type, call_args, ind)
    # TODO: just use Sym here
    if call_args === nothing
        :($Num($subscripted($type, $(Meta.quot(var_name)), $ind...)))
    elseif !isempty(call_args) && call_args[end] == :..
        :($Num($subscripted($FnType{Tuple{Any}, $type}, $(Meta.quot(var_name)), $ind...))) # XXX: using Num as output
    else
        :($Num($subscripted($FnType{NTuple{$(length(call_args)), Any}, $type}, $(Meta.quot(var_name)), $ind...)($(map(x->:($value($x)), call_args)...))))
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

TreeViews.hastreeview(x::Sym) = true
function TreeViews.treelabel(io::IO,x::Sym,
                             mime::MIME"text/plain" = MIME"text/plain"())
  show(io,mime,Text(x.name))
end

"""
varmap_to_vars(varmap,varlist)

Takes a list of pairs of variables=>values and an ordered list of variables and
creates the array of values in the correct order
"""
function varmap_to_vars(varmap::AbstractArray{<:Pair},varlist)
    out = similar(varmap,typeof(last(first(varmap))))
    for (ivar, ival) in varmap
        j = findfirst(isequal(ivar),varlist)
        if isnothing(j)
            throw(ArgumentError("Value $(ivar) provided in map not found in $(varlist)"))
        elseif j > length(varmap)
            throw(ArgumentError("Missing value in $(varmap), need $(varlist)"))
        end
        out[j] = ival
    end

    # Make output match varmap in type and shape
    # Does things like MArray->SArray
    ArrayInterface.restructure(varmap,out)
end
varmap_to_vars(varmap::AbstractArray,varlist) = varmap
varmap_to_vars(varmap::DiffEqBase.NullParameters,varlist) = varmap
varmap_to_vars(varmap::Nothing,varlist) = varmap
