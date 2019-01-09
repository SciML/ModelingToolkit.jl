using MacroTools
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

toexpr(ex) = MacroTools.postwalk(x->x isa Union{Expression,Operation} ? Expr(x) : x, ex)

is_constant(x::Variable) = x.subtype === :Constant
is_constant(::Any) = false

is_operation(::Operation) = true
is_operation(::Any) = false

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
function extract_elements!(op::Operation, elems, names, targetmap, default)
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
