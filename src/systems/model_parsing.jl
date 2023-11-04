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
    The dictionary with metadata like keyword arguements (:kwargs), base
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
    kwargs, ps, sps, vs, = [], [], [], []

    for arg in expr.args
        arg isa LineNumberNode && continue
        if arg.head == :macrocall
            parse_model!(exprs.args, comps, ext, eqs, icon, vs, ps,
                sps, dict, mod, arg, kwargs)
        elseif arg.head == :block
            push!(exprs.args, arg)
        elseif isconnector
            # Connectors can have variables listed without `@variables` prefix or
            # begin block.
            parse_variable_arg!(exprs, vs, dict, mod, arg, :variables, kwargs)
        else
            error("$arg is not valid syntax. Expected a macro call.")
        end
    end

    iv = get(dict, :independent_variable, nothing)
    if iv === nothing
        iv = dict[:independent_variable] = variable(:t)
    end

    gui_metadata = isassigned(icon) > 0 ? GUIMetadata(GlobalRef(mod, name), icon[]) :
                   GUIMetadata(GlobalRef(mod, name))

    sys = :($ODESystem($Equation[$(eqs...)], $iv, [$(vs...)], [$(ps...)];
        name, systems = [$(comps...)], gui_metadata = $gui_metadata))

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
                        dict[varclass][getname(var)][type] = mt
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
                        dict[varclass][getname(var)][type] = mt
                    end
                end
                var = set_var_metadata(var, meta)
            end
            (var, def)
        end
        Expr(:ref, a, b...) => begin
            parse_variable_def!(dict, mod, a, varclass, kwargs;
                def, indices = [eval.(b)...])
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
    @assert isequal(iv, prev_iv)
    vd = get!(dict, varclass) do
        Dict{Symbol, Dict{Symbol, Any}}()
    end
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
        _ => error("Cannot parse default $a")
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
        parse_icon!(icon, dict, body)
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

function parse_components!(exprs, cs, dict, body, kwargs)
    expr = Expr(:block)
    varexpr = Expr(:block)
    push!(exprs, varexpr)
    comps = Vector{Symbol}[]
    for arg in body.args
        arg isa LineNumberNode && continue
        MLStyle.@match arg begin
            Expr(:(=), a, b) => begin
                push!(cs, a)
                push!(comps, [a, b.args[1]])
                arg = deepcopy(arg)
                b = deepcopy(arg.args[2])

                component_args!(a, b, dict, expr, varexpr, kwargs)

                arg.args[2] = b
                push!(expr.args, arg)
            end
            _ => error("`@components` only takes assignment expressions. Got $arg")
        end
    end

    push!(exprs, :(@named $expr))

    dict[:components] = comps
end

function _rename(compname, varname)
    compname = Symbol(compname, :__, varname)
end

function component_args!(a, b, dict, expr, varexpr, kwargs)
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
                dict[:kwargs][_v] = nothing
            end
            Expr(:parameters, x...) => begin
                component_args!(a, arg, dict, expr, varexpr, kwargs)
            end
            Expr(:kw, x, y) => begin
                _v = _rename(a, x)
                b.args[i] = Expr(:kw, x, _v)
                push!(varexpr.args, :($_v = $_v === nothing ? $y : $_v))
                push!(kwargs, Expr(:kw, _v, nothing))
                dict[:kwargs][_v] = nothing
            end
            _ => error("Could not parse $arg of component $a")
        end
    end
end

function extend_args!(a, b, dict, expr, kwargs, varexpr, has_param = false)
    # Whenever `b` is a function call, skip the first arg aka the function name.
    # Whenver it is a kwargs list, include it.
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

function parse_variable_arg!(expr, vs, dict, mod, arg, varclass, kwargs)
    vv, def = parse_variable_def!(dict, mod, arg, varclass, kwargs)
    name = getname(vv)
    push!(expr.args,
        :($name = $name === nothing ?
                  $setdefault($vv, $def) :
                  $setdefault($vv, $name)))
    vv isa Num ? push!(vs, name) : push!(vs, :($name...))
end

function parse_variables!(exprs, vs, dict, mod, body, varclass, kwargs)
    expr = Expr(:block)
    push!(exprs, expr)
    for arg in body.args
        arg isa LineNumberNode && continue
        parse_variable_arg!(expr, vs, dict, mod, arg, varclass, kwargs)
    end
end

function parse_equations!(exprs, eqs, dict, body)
    for arg in body.args
        arg isa LineNumberNode && continue
        push!(eqs, arg)
    end
    # TODO: does this work with TOML?
    dict[:equations] = readable_code.(eqs)
end

function parse_icon!(icon, dict, body::String)
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

function parse_icon!(icon, dict, body::Expr)
    parse_icon!(icon, dict, eval(body))
end
