"""
$(TYPEDEF)

ModelingToolkit component or connector with metadata

# Fields
$(FIELDS)
"""
struct Model{F, S}
    """The constructor that returns ODESystem."""
    f::F
    """
    The dictionary with metadata like keyword arguments (:kwargs), base
    system this Model extends (:extend), sub-components of the Model (:components),
    variables (:variables), parameters (:parameters), structural parameters
    (:structural_parameters) and equations (:equations).
    """
    structure::S
    """
        This flag is `true` when the Model is a connector and is `false` when it is
        a component
        """
    isconnector::Bool
end
(m::Model)(args...; kw...) = m.f(args...; kw...)

for f in (:connector, :mtkmodel)
    isconnector = f == :connector ? true : false
    @eval begin
        macro $f(name::Symbol, body)
            esc($(:_model_macro)(__module__, name, body, $isconnector))
        end
    end
end

function _model_macro(mod, name, expr, isconnector)
    exprs = Expr(:block)
    dict = Dict{Symbol, Any}()
    dict[:kwargs] = Dict{Symbol, Any}()
    comps = Symbol[]
    ext = Ref{Any}(nothing)
    eqs = Expr[]
    icon = Ref{Union{String, URI}}()
    ps, sps, vs, = [], [], []
    kwargs = Set()

    push!(exprs.args, :(variables = []))
    push!(exprs.args, :(parameters = []))
    push!(exprs.args, :(systems = ODESystem[]))
    push!(exprs.args, :(equations = Equation[]))

    Base.remove_linenums!(expr)
    for arg in expr.args
        if arg.head == :macrocall
            parse_model!(exprs.args, comps, ext, eqs, icon, vs, ps,
                sps, dict, mod, arg, kwargs)
        elseif arg.head == :block
            push!(exprs.args, arg)
        elseif arg.head == :if
            MLStyle.@match arg begin
                Expr(:if, condition, x) => begin
                    parse_conditional_model_statements(comps, dict, eqs, exprs, kwargs,
                        mod, ps, vs, parse_top_level_branch(condition, x.args)...)
                end
                Expr(:if, condition, x, y) => begin
                    parse_conditional_model_statements(comps, dict, eqs, exprs, kwargs,
                        mod, ps, vs, parse_top_level_branch(condition, x.args, y)...)
                end
                _ => error("Got an invalid argument: $arg")
            end
        elseif isconnector
            # Connectors can have variables listed without `@variables` prefix or
            # begin block.
            parse_variable_arg!(exprs.args, vs, dict, mod, arg, :variables, kwargs)
        else
            error("$arg is not valid syntax. Expected a macro call.")
        end
    end

    iv = get(dict, :independent_variable, nothing)
    if iv === nothing
        iv = dict[:independent_variable] = variable(:t)
    end

    push!(exprs.args, :(push!(equations, $(eqs...))))
    push!(exprs.args, :(push!(parameters, $(ps...))))
    push!(exprs.args, :(push!(systems, $(comps...))))
    push!(exprs.args, :(push!(variables, $(vs...))))

    gui_metadata = isassigned(icon) > 0 ? GUIMetadata(GlobalRef(mod, name), icon[]) :
                   GUIMetadata(GlobalRef(mod, name))

    sys = :($ODESystem($Equation[equations...], $iv, variables, parameters;
        name, systems, gui_metadata = $gui_metadata))

    if ext[] === nothing
        push!(exprs.args, :(var"#___sys___" = $sys))
    else
        push!(exprs.args, :(var"#___sys___" = $extend($sys, $(ext[]))))
    end

    isconnector && push!(exprs.args,
        :($Setfield.@set!(var"#___sys___".connector_type=$connector_type(var"#___sys___"))))

    f = :($(Symbol(:__, name, :__))(; name, $(kwargs...)) = $exprs)
    :($name = $Model($f, $dict, $isconnector))
end

function parse_variable_def!(dict, mod, arg, varclass, kwargs;
        def = nothing, indices::Union{Vector{UnitRange{Int}}, Nothing} = nothing)
    metatypes = [(:connection_type, VariableConnectType),
        (:description, VariableDescription),
        (:unit, VariableUnit),
        (:bounds, VariableBounds),
        (:noise, VariableNoiseType),
        (:input, VariableInput),
        (:output, VariableOutput),
        (:irreducible, VariableIrreducible),
        (:state_priority, VariableStatePriority),
        (:misc, VariableMisc),
        (:disturbance, VariableDisturbance),
        (:tunable, VariableTunable),
        (:dist, VariableDistribution),
        (:binary, VariableBinary),
        (:integer, VariableInteger)]

    arg isa LineNumberNode && return
    MLStyle.@match arg begin
        a::Symbol => begin
            push!(kwargs, Expr(:kw, a, nothing))
            var = generate_var!(dict, a, varclass; indices)
            dict[:kwargs][getname(var)] = def
            (var, def)
        end
        Expr(:call, a, b) => begin
            push!(kwargs, Expr(:kw, a, nothing))
            var = generate_var!(dict, a, b, varclass; indices)
            dict[:kwargs][getname(var)] = def
            (var, def)
        end
        Expr(:(=), a, b) => begin
            Base.remove_linenums!(b)
            def, meta = parse_default(mod, b)
            var, def = parse_variable_def!(dict, mod, a, varclass, kwargs; def)
            dict[varclass][getname(var)][:default] = def
            if meta !== nothing
                for (type, key) in metatypes
                    if (mt = get(meta, key, nothing)) !== nothing
                        key == VariableConnectType && (mt = nameof(mt))
                        if dict[varclass] isa Vector
                            dict[varclass][1][getname(var)][type] = mt
                        else
                            dict[varclass][getname(var)][type] = mt
                        end
                    end
                end
                var = set_var_metadata(var, meta)
            end
            (var, def)
        end
        Expr(:tuple, a, b) => begin
            var, def = parse_variable_def!(dict, mod, a, varclass, kwargs)
            meta = parse_metadata(mod, b)
            if meta !== nothing
                for (type, key) in metatypes
                    if (mt = get(meta, key, nothing)) !== nothing
                        key == VariableConnectType && (mt = nameof(mt))
                        # @info dict 164
                        if dict[varclass] isa Vector
                            dict[varclass][1][getname(var)][type] = mt
                        else
                            dict[varclass][getname(var)][type] = mt
                        end
                    end
                end
                var = set_var_metadata(var, meta)
            end
            (var, def)
        end
        Expr(:ref, a, b...) => begin
            indices = map(i -> UnitRange(i.args[2], i.args[end]), b)
            parse_variable_def!(dict, mod, a, varclass, kwargs;
                def, indices)
        end
        _ => error("$arg cannot be parsed")
    end
end

function generate_var(a, varclass;
        indices::Union{Vector{UnitRange{Int}}, Nothing} = nothing)
    var = indices === nothing ? Symbolics.variable(a) : first(@variables $a[indices...])
    if varclass == :parameters
        var = toparam(var)
    end
    var
end

function generate_var!(dict, a, varclass;
        indices::Union{Vector{UnitRange{Int}}, Nothing} = nothing)
    vd = get!(dict, varclass) do
        Dict{Symbol, Dict{Symbol, Any}}()
    end
    vd isa Vector && (vd = first(vd))
    vd[a] = Dict{Symbol, Any}()
    indices !== nothing && (vd[a][:size] = Tuple(lastindex.(indices)))
    generate_var(a, varclass; indices)
end

function generate_var!(dict, a, b, varclass;
        indices::Union{Vector{UnitRange{Int}}, Nothing} = nothing)
    iv = generate_var(b, :variables)
    prev_iv = get!(dict, :independent_variable) do
        iv
    end
    @assert isequal(iv, prev_iv) "Multiple independent variables are used in the model"
    vd = get!(dict, varclass) do
        Dict{Symbol, Dict{Symbol, Any}}()
    end
    vd isa Vector && (vd = first(vd))
    vd[a] = Dict{Symbol, Any}()
    var = if indices === nothing
        Symbolics.variable(a, T = SymbolicUtils.FnType{Tuple{Real}, Real})(iv)
    else
        vd[a][:size] = Tuple(lastindex.(indices))
        first(@variables $a(iv)[indices...])
    end
    if varclass == :parameters
        var = toparam(var)
    end
    var
end

function parse_default(mod, a)
    a = Base.remove_linenums!(deepcopy(a))
    MLStyle.@match a begin
        Expr(:block, x) => parse_default(mod, x)
        Expr(:tuple, x, y) => begin
            def, _ = parse_default(mod, x)
            meta = parse_metadata(mod, y)
            (def, meta)
        end
        ::Symbol || ::Number => (a, nothing)
        Expr(:call, a...) => begin
            def = parse_default.(Ref(mod), a)
            expr = Expr(:call)
            for (d, _) in def
                push!(expr.args, d)
            end
            (expr, nothing)
        end
        Expr(:if, condition, x, y) => (a, nothing)
        _ => error("Cannot parse default $a $(typeof(a))")
    end
end

function parse_metadata(mod, a)
    MLStyle.@match a begin
        Expr(:vect, eles...) => Dict(parse_metadata(mod, e) for e in eles)
        Expr(:(=), a, b) => Symbolics.option_to_metadata_type(Val(a)) => get_var(mod, b)
        _ => error("Cannot parse metadata $a")
    end
end

function set_var_metadata(a, ms)
    for (m, v) in ms
        a = wrap(set_scalar_metadata(unwrap(a), m, v))
    end
    a
end

function get_var(mod::Module, b)
    if b isa Symbol
        getproperty(mod, b)
    elseif b isa Expr
        Core.eval(mod, b)
    else
        b
    end
end

function parse_model!(exprs, comps, ext, eqs, icon, vs, ps, sps,
        dict, mod, arg, kwargs)
    mname = arg.args[1]
    body = arg.args[end]
    if mname == Symbol("@components")
        parse_components!(exprs, comps, dict, body, kwargs)
    elseif mname == Symbol("@extend")
        parse_extend!(exprs, ext, dict, mod, body, kwargs)
    elseif mname == Symbol("@variables")
        parse_variables!(exprs, vs, dict, mod, body, :variables, kwargs)
    elseif mname == Symbol("@parameters")
        parse_variables!(exprs, ps, dict, mod, body, :parameters, kwargs)
    elseif mname == Symbol("@structural_parameters")
        parse_structural_parameters!(exprs, sps, dict, mod, body, kwargs)
    elseif mname == Symbol("@equations")
        parse_equations!(exprs, eqs, dict, body)
    elseif mname == Symbol("@icon")
        isassigned(icon) && error("This model has more than one icon.")
        parse_icon!(body, dict, icon, mod)
    else
        error("$mname is not handled.")
    end
end

function parse_structural_parameters!(exprs, sps, dict, mod, body, kwargs)
    Base.remove_linenums!(body)
    for arg in body.args
        MLStyle.@match arg begin
            Expr(:(=), a, b) => begin
                push!(sps, a)
                push!(kwargs, Expr(:kw, a, b))
                dict[:kwargs][a] = b
            end
            a => begin
                push!(sps, a)
                push!(kwargs, a)
                dict[:kwargs][a] = nothing
            end
        end
    end
end

function extend_args!(a, b, dict, expr, kwargs, varexpr, has_param = false)
    # Whenever `b` is a function call, skip the first arg aka the function name.
    # Whenever it is a kwargs list, include it.
    start = b.head == :call ? 2 : 1
    for i in start:lastindex(b.args)
        arg = b.args[i]
        arg isa LineNumberNode && continue
        MLStyle.@match arg begin
            x::Symbol => begin
                if b.head != :parameters
                    if has_param
                        popat!(b.args, i)
                        push!(b.args[2].args, x)
                    else
                        b.args[i] = Expr(:parameters, x)
                    end
                end
                push!(kwargs, Expr(:kw, x, nothing))
                dict[:kwargs][x] = nothing
            end
            Expr(:kw, x) => begin
                push!(kwargs, Expr(:kw, x, nothing))
                dict[:kwargs][x] = nothing
            end
            Expr(:kw, x, y) => begin
                b.args[i] = Expr(:kw, x, x)
                push!(varexpr.args, :($x = $x === nothing ? $y : $x))
                push!(kwargs, Expr(:kw, x, nothing))
                dict[:kwargs][x] = nothing
            end
            Expr(:parameters, x...) => begin
                has_param = true
                extend_args!(a, arg, dict, expr, kwargs, varexpr, has_param)
            end
            _ => error("Could not parse $arg of component $a")
        end
    end
end

const EMPTY_DICT = Dict()
const EMPTY_VoVoSYMBOL = Vector{Symbol}[]

function Base.names(model::Model)
    vars = keys(get(model.structure, :variables, EMPTY_DICT))
    vars = union(vars, keys(get(model.structure, :parameters, EMPTY_DICT)))
    vars = union(vars,
        map(first, get(model.structure, :components, EMPTY_VoVoSYMBOL)))
    collect(vars)
end

function _parse_extend!(ext, a, b, dict, expr, kwargs, varexpr, vars)
    extend_args!(a, b, dict, expr, kwargs, varexpr)
    ext[] = a
    push!(b.args, Expr(:kw, :name, Meta.quot(a)))
    push!(expr.args, :($a = $b))

    dict[:extend] = [Symbol.(vars.args), a, b.args[1]]

    push!(expr.args, :(@unpack $vars = $a))
end

function parse_extend!(exprs, ext, dict, mod, body, kwargs)
    expr = Expr(:block)
    varexpr = Expr(:block)
    push!(exprs, varexpr)
    push!(exprs, expr)
    body = deepcopy(body)
    MLStyle.@match body begin
        Expr(:(=), a, b) => begin
            if Meta.isexpr(b, :(=))
                vars = a
                if !Meta.isexpr(vars, :tuple)
                    error("`@extend` destructuring only takes an tuple as LHS. Got $body")
                end
                a, b = b.args
                _parse_extend!(ext, a, b, dict, expr, kwargs, varexpr, vars)
            else
                error("When explicitly destructing in `@extend` please use the syntax: `@extend a, b = oneport = OnePort()`.")
            end
        end
        Expr(:call, a′, _...) => begin
            a = Symbol(Symbol("#mtkmodel"), :__anonymous__, a′)
            b = body
            if (model = getproperty(mod, b.args[1])) isa Model
                vars = Expr(:tuple)
                append!(vars.args, names(model))
                _parse_extend!(ext, a, b, dict, expr, kwargs, varexpr, vars)
            else
                error("Cannot infer the exact `Model` that `@extend $(body)` refers." *
                      " Please specify the names that it brings into scope by:" *
                      " `@extend a, b = oneport = OnePort()`.")
            end
        end
        _ => error("`@extend` only takes an assignment expression. Got $body")
    end
    return nothing
end

function parse_variable_arg!(exprs, vs, dict, mod, arg, varclass, kwargs)
    name, ex = parse_variable_arg(dict, mod, arg, varclass, kwargs)
    push!(vs, name)
    push!(exprs, ex)
end

function parse_variable_arg(dict, mod, arg, varclass, kwargs)
    vv, def = parse_variable_def!(dict, mod, arg, varclass, kwargs)
    name = getname(vv)
    return vv isa Num ? name : :($name...),
    :($name = $name === nothing ? $setdefault($vv, $def) : $setdefault($vv, $name))
end

function handle_conditional_vars!(arg, conditional_branch, mod, varclass, kwargs)
    conditional_dict = Dict(:kwargs => Dict(),
        :parameters => Any[Dict{Symbol, Dict{Symbol, Any}}()],
        :variables => Any[Dict{Symbol, Dict{Symbol, Any}}()])
    for _arg in arg.args
        name, ex = parse_variable_arg(conditional_dict, mod, _arg, varclass, kwargs)
        push!(conditional_branch.args, ex)
        push!(conditional_branch.args, :(push!($varclass, $name)))
    end
    conditional_dict
end

function prune_conditional_dict!(conditional_tuple::Tuple)
    prune_conditional_dict!.(collect(conditional_tuple))
end
function prune_conditional_dict!(conditional_dict::Dict)
    for k in [:parameters, :variables]
        length(conditional_dict[k]) == 1 && isempty(first(conditional_dict[k])) &&
            delete!(conditional_dict, k)
    end
    isempty(conditional_dict[:kwargs]) && delete!(conditional_dict, :kwargs)
end
prune_conditional_dict!(_) = return nothing

function get_conditional_dict!(conditional_dict, conditional_y_tuple::Tuple)
    k = get_conditional_dict!.(Ref(conditional_dict), collect(conditional_y_tuple))
    push_something!(conditional_dict,
        k...)
    conditional_dict
end

function get_conditional_dict!(conditional_dict::Dict, conditional_y_tuple::Dict)
    merge!(conditional_dict[:kwargs], conditional_y_tuple[:kwargs])
    for key in [:parameters, :variables]
        merge!(conditional_dict[key][1], conditional_y_tuple[key][1])
    end
    conditional_dict
end

get_conditional_dict!(a, b) = (return nothing)

function push_conditional_dict!(dict, condition, conditional_dict,
        conditional_y_tuple, varclass)
    vd = get!(dict, varclass) do
        Dict{Symbol, Dict{Symbol, Any}}()
    end
    for k in keys(conditional_dict[varclass][1])
        vd[k] = copy(conditional_dict[varclass][1][k])
        vd[k][:condition] = (:if, condition, conditional_dict, conditional_y_tuple)
    end
    conditional_y_dict = Dict(:kwargs => Dict(),
        :parameters => Any[Dict{Symbol, Dict{Symbol, Any}}()],
        :variables => Any[Dict{Symbol, Dict{Symbol, Any}}()])
    get_conditional_dict!(conditional_y_dict, conditional_y_tuple)

    prune_conditional_dict!(conditional_y_dict)
    prune_conditional_dict!(conditional_dict)
    !isempty(conditional_y_dict) && for k in keys(conditional_y_dict[varclass][1])
        vd[k] = copy(conditional_y_dict[varclass][1][k])
        vd[k][:condition] = (:if, condition, conditional_dict, conditional_y_tuple)
    end
end

function parse_variables!(exprs, vs, dict, mod, body, varclass, kwargs)
    expr = Expr(:block)
    push!(exprs, expr)
    for arg in body.args
        arg isa LineNumberNode && continue
        MLStyle.@match arg begin
            Expr(:if, condition, x) => begin
                conditional_expr = Expr(:if, condition, Expr(:block))
                conditional_dict = handle_conditional_vars!(x,
                    conditional_expr.args[2],
                    mod,
                    varclass,
                    kwargs)
                push!(expr.args, conditional_expr)
                push_conditional_dict!(dict, condition, conditional_dict, nothing, varclass)
            end
            Expr(:if, condition, x, y) => begin
                conditional_expr = Expr(:if, condition, Expr(:block))
                conditional_dict = handle_conditional_vars!(x,
                    conditional_expr.args[2],
                    mod,
                    varclass,
                    kwargs)
                conditional_y_expr, conditional_y_tuple = handle_y_vars(y,
                    conditional_dict,
                    mod,
                    varclass,
                    kwargs)
                push!(conditional_expr.args, conditional_y_expr)
                push!(expr.args, conditional_expr)
                push_conditional_dict!(dict,
                    condition,
                    conditional_dict,
                    conditional_y_tuple,
                    varclass)
            end
            _ => parse_variable_arg!(exprs, vs, dict, mod, arg, varclass, kwargs)
        end
    end
end

function handle_y_vars(y, dict, mod, varclass, kwargs)
    conditional_dict = if Meta.isexpr(y, :elseif)
        conditional_y_expr = Expr(:elseif, y.args[1], Expr(:block))
        conditional_dict = handle_conditional_vars!(y.args[2],
            conditional_y_expr.args[2],
            mod,
            varclass,
            kwargs)
        _y_expr, _conditional_dict = handle_y_vars(y.args[end], dict, mod, varclass, kwargs)
        push!(conditional_y_expr.args, _y_expr)
        (:elseif, y.args[1], conditional_dict, _conditional_dict)
    else
        conditional_y_expr = Expr(:block)
        handle_conditional_vars!(y, conditional_y_expr, mod, varclass, kwargs)
    end
    conditional_y_expr, conditional_dict
end

function handle_if_x_equations!(condition, dict, ifexpr, x)
    if Meta.isexpr(x, :block)
        push!(ifexpr.args, condition, :(push!(equations, $(x.args...))))
        return readable_code.(x.args)
    else
        push!(ifexpr.args, condition, :(push!(equations, $x)))
        return readable_code(x)
    end
    # push!(dict[:equations], [:if, readable_code(condition), readable_code.(x.args)])
end

function handle_if_y_equations!(ifexpr, y, dict)
    if y.head == :elseif
        elseifexpr = Expr(:elseif)
        eq_entry = [:elseif, readable_code.(y.args[1].args)...]
        push!(eq_entry, handle_if_x_equations!(y.args[1], dict, elseifexpr, y.args[2]))
        get(y.args, 3, nothing) !== nothing &&
            push!(eq_entry, handle_if_y_equations!(elseifexpr, y.args[3], dict))
        push!(ifexpr.args, elseifexpr)
        (eq_entry...,)
    else
        if Meta.isexpr(y, :block)
            push!(ifexpr.args, :(push!(equations, $(y.args...))))
        else
            push!(ifexpr.args, :(push!(equations, $(y))))
        end
        readable_code.(y.args)
    end
end

function parse_equations!(exprs, eqs, dict, body)
    dict[:equations] = []
    Base.remove_linenums!(body)
    for arg in body.args
        MLStyle.@match arg begin
            Expr(:if, condition, x) => begin
                ifexpr = Expr(:if)
                eq_entry = handle_if_x_equations!(condition, dict, ifexpr, x)
                push!(exprs, ifexpr)
                push!(dict[:equations], (:if, condition, eq_entry))
            end
            Expr(:if, condition, x, y) => begin
                ifexpr = Expr(:if)
                xeq_entry = handle_if_x_equations!(condition, dict, ifexpr, x)
                yeq_entry = handle_if_y_equations!(ifexpr, y, dict)
                push!(exprs, ifexpr)
                push!(dict[:equations], (:if, condition, xeq_entry, yeq_entry))
            end
            _ => begin
                push!(eqs, arg)
                push!(dict[:equations], readable_code.(eqs)...)
            end
        end
    end
end

function parse_icon!(body::String, dict, icon, mod)
    icon_dir = get(ENV, "MTK_ICONS_DIR", joinpath(DEPOT_PATH[1], "mtk_icons"))
    dict[:icon] = icon[] = if isfile(body)
        URI("file:///" * abspath(body))
    elseif (iconpath = joinpath(icon_dir, body); isfile(iconpath))
        URI("file:///" * abspath(iconpath))
    elseif try
        Base.isvalid(URI(body))
    catch e
        false
    end
        URI(body)
    elseif (_body = lstrip(body); startswith(_body, r"<\?xml|<svg"))
        String(_body) # With Julia-1.10 promoting `SubString{String}` to `String` can be dropped.
    else
        error("\n$body is not a valid icon")
    end
end

function parse_icon!(body::Symbol, dict, icon, mod)
    parse_icon!(getfield(mod, body), dict, icon, mod)
end

### Parsing Components:

function component_args!(a, b, expr, varexpr, kwargs)
    # Whenever `b` is a function call, skip the first arg aka the function name.
    # Whenever it is a kwargs list, include it.
    start = b.head == :call ? 2 : 1
    for i in start:lastindex(b.args)
        arg = b.args[i]
        arg isa LineNumberNode && continue
        MLStyle.@match arg begin
            x::Symbol || Expr(:kw, x) => begin
                _v = _rename(a, x)
                b.args[i] = Expr(:kw, x, _v)
                push!(varexpr.args, :((@isdefined $x) && ($_v = $x)))
                push!(kwargs, Expr(:kw, _v, nothing))
                # dict[:kwargs][_v] = nothing
            end
            Expr(:parameters, x...) => begin
                component_args!(a, arg, expr, varexpr, kwargs)
            end
            Expr(:kw, x, y) => begin
                _v = _rename(a, x)
                b.args[i] = Expr(:kw, x, _v)
                push!(varexpr.args, :($_v = $_v === nothing ? $y : $_v))
                push!(kwargs, Expr(:kw, _v, nothing))
                # dict[:kwargs][_v] = nothing
            end
            _ => error("Could not parse $arg of component $a")
        end
    end
end

function _parse_components!(exprs, body, kwargs)
    expr = Expr(:block)
    varexpr = Expr(:block)
    # push!(exprs, varexpr)
    comps = Vector{Union{Symbol, Expr}}[]
    comp_names = []

    for arg in body.args
        arg isa LineNumberNode && continue
        MLStyle.@match arg begin
            Expr(:block) => begin
                # TODO: Do we need this?
                error("Multiple `@components` block detected within a single block")
            end
            Expr(:(=), a, b) => begin
                arg = deepcopy(arg)
                b = deepcopy(arg.args[2])

                component_args!(a, b, expr, varexpr, kwargs)

                arg.args[2] = b
                push!(expr.args, arg)
                push!(comp_names, a)
                if (isa(b.args[1], Symbol) || Meta.isexpr(b.args[1], :.))
                    push!(comps, [a, b.args[1]])
                end
            end
            _ => error("Couldn't parse the component body: $arg")
        end
    end
    return comp_names, comps, expr, varexpr
end

function push_conditional_component!(ifexpr, expr_vec, comp_names, varexpr)
    blk = Expr(:block)
    push!(blk.args, varexpr)
    push!(blk.args, :(@named begin
        $(expr_vec.args...)
    end))
    push!(blk.args, :($push!(systems, $(comp_names...))))
    push!(ifexpr.args, blk)
end

function handle_if_x!(mod, exprs, ifexpr, x, kwargs, condition = nothing)
    push!(ifexpr.args, condition)
    comp_names, comps, expr_vec, varexpr = _parse_components!(ifexpr, x, kwargs)
    push_conditional_component!(ifexpr, expr_vec, comp_names, varexpr)
    comps
end

function handle_if_y!(exprs, ifexpr, y, kwargs)
    Base.remove_linenums!(y)
    if Meta.isexpr(y, :elseif)
        comps = [:elseif, y.args[1]]
        elseifexpr = Expr(:elseif)
        push!(comps, handle_if_x!(mod, exprs, elseifexpr, y.args[2], kwargs, y.args[1]))
        get(y.args, 3, nothing) !== nothing &&
            push!(comps, handle_if_y!(exprs, elseifexpr, y.args[3], kwargs))
        push!(ifexpr.args, elseifexpr)
        (comps...,)
    else
        comp_names, comps, expr_vec, varexpr = _parse_components!(exprs, y, kwargs)
        push_conditional_component!(ifexpr, expr_vec, comp_names, varexpr)
        comps
    end
end

function handle_conditional_components(condition, dict, exprs, kwargs, x, y = nothing)
    ifexpr = Expr(:if)
    comps = handle_if_x!(mod, exprs, ifexpr, x, kwargs, condition)
    ycomps = y === nothing ? [] : handle_if_y!(exprs, ifexpr, y, kwargs)
    push!(exprs, ifexpr)
    push!(dict[:components], (:if, condition, comps, ycomps))
end

function parse_components!(exprs, cs, dict, compbody, kwargs)
    dict[:components] = []
    Base.remove_linenums!(compbody)
    for arg in compbody.args
        MLStyle.@match arg begin
            Expr(:if, condition, x) => begin
                handle_conditional_components(condition, dict, exprs, kwargs, x)
            end
            Expr(:if, condition, x, y) => begin
                handle_conditional_components(condition, dict, exprs, kwargs, x, y)
            end
            Expr(:(=), a, b) => begin
                comp_names, comps, expr_vec, varexpr = _parse_components!(exprs,
                    :(begin
                        $arg
                    end),
                    kwargs)
                push!(cs, comp_names...)
                push!(dict[:components], comps...)
                push!(exprs, varexpr, :(@named begin
                    $(expr_vec.args...)
                end))
            end
            _ => error("Couldn't parse the component body $compbody")
        end
    end
end

function _rename(compname, varname)
    compname = Symbol(compname, :__, varname)
end

# Handle top level branching
push_something!(v, ::Nothing) = v
push_something!(v, x) = push!(v, x)
push_something!(v::Dict, x::Dict) = merge!(v, x)
push_something!(v, x...) = push_something!.(Ref(v), x)

define_blocks(branch) = [Expr(branch), Expr(branch), Expr(branch), Expr(branch)]

function parse_top_level_branch(condition, x, y = nothing, branch = :if)
    blocks::Vector{Union{Expr, Nothing}} = component_blk, equations_blk, parameter_blk, variable_blk = define_blocks(branch)

    for arg in x
        if arg.args[1] == Symbol("@components")
            push_something!(component_blk.args, condition, arg.args[end])
        elseif arg.args[1] == Symbol("@equations")
            push_something!(equations_blk.args, condition, arg.args[end])
        elseif arg.args[1] == Symbol("@variables")
            push_something!(variable_blk.args, condition, arg.args[end])
        elseif arg.args[1] == Symbol("@parameters")
            push_something!(parameter_blk.args, condition, arg.args[end])
        else
            error("$(arg.args[1]) isn't supported")
        end
    end

    if y !== nothing
        yblocks = if y.head == :elseif
            parse_top_level_branch(y.args[1],
                y.args[2].args,
                lastindex(y.args) == 3 ? y.args[3] : nothing,
                :elseif)
        else
            yblocks = parse_top_level_branch(nothing, y.args, nothing, :block)

            for i in 1:lastindex(yblocks)
                yblocks[i] !== nothing && (yblocks[i] = yblocks[i].args[end])
            end
            yblocks
        end
        for i in 1:lastindex(yblocks)
            if lastindex(blocks[i].args) == 1
                push_something!(blocks[i].args, Expr(:block), yblocks[i])
            elseif lastindex(blocks[i].args) == 0
                blocks[i] = yblocks[i]
            else
                push_something!(blocks[i].args, yblocks[i])
            end
        end
    end

    for i in 1:lastindex(blocks)
        blocks[i] !== nothing && isempty(blocks[i].args) && (blocks[i] = nothing)
    end

    return blocks
end

function parse_conditional_model_statements(comps, dict, eqs, exprs, kwargs, mod,
        ps, vs, component_blk, equations_blk, parameter_blk, variable_blk)
    parameter_blk !== nothing &&
        parse_variables!(exprs.args, ps, dict, mod, :(begin
                $parameter_blk
            end), :parameters, kwargs)

    variable_blk !== nothing &&
        parse_variables!(exprs.args, vs, dict, mod, :(begin
                $variable_blk
            end), :variables, kwargs)

    component_blk !== nothing &&
        parse_components!(exprs.args,
            comps, dict, :(begin
                $component_blk
            end), kwargs)

    equations_blk !== nothing &&
        parse_equations!(exprs.args, eqs, dict, :(begin
            $equations_blk
        end))
end
