using MacroTools


Base.convert(::Type{Expression}, ex::Expr) = convert(Term, ex)
Base.convert(::Type{Expression}, x::Expression) = x
Base.convert(::Type{Expression}, sym::Symbol, args) = convert(Term, Expr(:call, eval(sym), args...))
Base.convert(::Type{Expression}, x::Variable) = x
Base.convert(::Type{Expression}, x::Number) = convert(Term, x)


function expr_arr_to_block(exprs)
    block = :(begin end)
    foreach(expr -> push!(block.args, expr), exprs)
    block
end

# used in parsing
isblock(x) = length(x) == 1 && x[1] isa Expr && x[1].head == :block
function flatten_expr!(x)
    isb = isblock(x)
    if isb
        x = MacroTools.striplines(x[1])
        filter!(z->z isa Symbol || z.head != :line, x.args)
        x = (x.args...,)
    end
    x
end

toexpr(ex) = MacroTools.postwalk(x -> isa(x, Expression) ? convert(Expr, x) : x, ex)

function partition(f, xs)
    idxs = map(f, xs)
    not_idxs = eachindex(xs) .∉ (idxs,)
    return (xs[idxs], xs[not_idxs])
end

function unpack(t::Term)
    @assert root(t) === :call
    _children = children(t)
    fn, args = root(_children[1]), _children[2:end]
    return (fn, args)
end

pack(fn, args) = convert(Term, :($fn($(args...))))

is_constant(t::Term) = !is_branch(t) && !isa(root(t), Variable)
is_constant(::Any) = false

has_dependent(t::Variable) = Base.Fix2(has_dependent, t)
has_dependent(x::Variable, t::Variable) =
    t ∈ x.dependents || any(has_dependent(t), x.dependents)


extract_idv(eq) = eq.args[1].diff.x

function extract_elements(ops, targetmap, default = nothing)
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
    for op in ops
        extract_elements!(op, elems, names, targetmap, default)
    end
    Tuple(elems[target] for target in targets)
end

# Walk the tree recursively and push variables into the right set
function extract_elements!(op::Term, elems, names, targetmap, default)  # FIXME
    for arg in op.args
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
