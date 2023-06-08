macro connector(name::Symbol, body)
    esc(connector_macro(__module__, name, body))
end

struct Model{F, S}
    f::F
    structure::S
end
(m::Model)(args...; kw...) = m.f(args...; kw...)

using MLStyle

function connector_macro(mod, name, body)
    if !Meta.isexpr(body, :block)
        err = """
        connector body must be a block! It should be in the form of
        ```
        @connector Pin begin
            v(t) = 1
            (i(t) = 1), [connect = Flow]
        end
        ```
        """
        error(err)
    end
    vs = Num[]
    icon = Ref{Union{String, URI}}()
    dict = Dict{Symbol, Any}()
    for arg in body.args
        arg isa LineNumberNode && continue
        if arg.head == :macrocall && arg.args[1] == Symbol("@icon")
            parse_icon!(icon, dict, dict, arg.args[end])
            continue
        end
        push!(vs, Num(parse_variable_def!(dict, mod, arg, :variables, kwargs)))
    end
    iv = get(dict, :independent_variable, nothing)
    if iv === nothing
        error("$name doesn't have a independent variable")
    end
    gui_metadata = isassigned(icon) ? GUIMetadata(GlobalRef(mod, name), icon[]) :
                   nothing
    quote
        $name = $Model((; name) -> begin
                var"#___sys___" = $ODESystem($(Equation[]), $iv, $vs, $([]);
                    name, gui_metadata = $gui_metadata)
                $Setfield.@set!(var"#___sys___".connector_type=$connector_type(var"#___sys___"))
            end, $dict)
    end
end

function parse_variable_def!(dict, mod, arg, varclass, kwargs)
    MLStyle.@match arg begin
        ::Symbol => generate_var!(dict, arg, varclass)
        Expr(:call, a, b) => generate_var!(dict, a, b, varclass)
        Expr(:(=), a, b) => begin
            var = parse_variable_def!(dict, mod, a, varclass, kwargs)
            def = parse_default(mod, b, kwargs)
            dict[varclass][getname(var)][:default] = def
            var = setdefault(var, def)
            var
        end
        Expr(:tuple, a, b) => begin
            var = parse_variable_def!(dict, mod, a, varclass, kwargs)
            meta = parse_metadata(mod, b)
            if (ct = get(meta, VariableConnectType, nothing)) !== nothing
                dict[varclass][getname(var)][:connection_type] = nameof(ct)
            end
            set_var_metadata(var, meta)
        end
        _ => error("$arg cannot be parsed")
    end
end

function generate_var(a, varclass)
    var = Symbolics.variable(a)
    if varclass == :parameters
        var = toparam(var)
    end
    var
end
function generate_var!(dict, a, varclass)
    var = generate_var(a, varclass)
    vd = get!(dict, varclass) do
        Dict{Symbol, Dict{Symbol, Any}}()
    end
    vd[a] = Dict{Symbol, Any}()
    var
end
function generate_var!(dict, a, b, varclass)
    iv = generate_var(b, :variables)
    prev_iv = get!(dict, :independent_variable) do
        iv
    end
    @assert isequal(iv, prev_iv)
    vd = get!(dict, varclass) do
        Dict{Symbol, Dict{Symbol, Any}}()
    end
    vd[a] = Dict{Symbol, Any}()
    var = Symbolics.variable(a, T = SymbolicUtils.FnType{Tuple{Real}, Real})(iv)
    if varclass == :parameters
        var = toparam(var)
    end
    var
end

function set_kwargs!(kwargs, a)
    push!(kwargs, a)
    return a
end

function parse_default(mod, a, kwargs)
    a = Base.remove_linenums!(deepcopy(a))
    MLStyle.@match a begin
        Expr(:block, a) => set_kwargs!(kwargs, a)
        ::Symbol => get_var(mod, a)
        ::Number => a
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
        a = setmetadata(a, m, v)
    end
    a
end
function get_var(mod::Module, b)
    b isa Symbol ? getproperty(mod, b) : b
end

macro model(name::Symbol, expr)
    esc(model_macro(__module__, name, expr))
end

function model_macro(mod, name, expr)
    exprs = Expr(:block)
    dict = Dict{Symbol, Any}()
    comps = Symbol[]
    ext = Ref{Any}(nothing)
    vs = Symbol[]
    ps = Symbol[]
    eqs = Expr[]
    icon = Ref{Union{String, URI}}()
    kwargs = []
    for arg in expr.args
        arg isa LineNumberNode && continue
        arg.head == :macrocall || error("$arg is not valid syntax. Expected a macro call.")
        parse_model!(exprs.args, comps, ext, eqs, vs, ps, icon, dict, mod, arg, kwargs)
    end
    iv = get(dict, :independent_variable, nothing)
    if iv === nothing
        iv = dict[:independent_variable] = variable(:t)
    end
    gui_metadata = isassigned(icon) > 0 ? GUIMetadata(GlobalRef(mod, name), icon[]) :
                   nothing
    sys = :($ODESystem($Equation[$(eqs...)], $iv, [$(vs...)], [$(ps...)];
        systems = [$(comps...)], name, gui_metadata = $gui_metadata))
    if ext[] === nothing
        push!(exprs.args, sys)
    else
        push!(exprs.args, :($extend($sys, $(ext[]))))
    end
    @info "Exprs: $exprs"
    :($name = $Model((; name, kwargs...) -> $exprs, $dict))
end

function parse_model!(exprs, comps, ext, eqs, vs, ps, icon, dict, mod, arg, kwargs)
    mname = arg.args[1]
    body = arg.args[end]
    if mname == Symbol("@components")
        parse_components!(exprs, comps, dict, body)
    elseif mname == Symbol("@extend")
        parse_extend!(exprs, ext, dict, body)
    elseif mname == Symbol("@variables")
        parse_variables!(exprs, vs, dict, mod, body, :variables, kwargs)
    elseif mname == Symbol("@parameters")
        parse_variables!(exprs, ps, dict, mod, body, :parameters, kwargs)
    elseif mname == Symbol("@equations")
        parse_equations!(exprs, eqs, dict, body)
    elseif mname == Symbol("@icon")
        parse_icon!(icon, dict, mod, body)
    else
        error("$mname is not handled.")
    end
end

function parse_components!(exprs, cs, dict, body)
    expr = Expr(:block)
    push!(exprs, expr)
    comps = Vector{String}[]
    for arg in body.args
        arg isa LineNumberNode && continue
        MLStyle.@match arg begin
            Expr(:(=), a, b) => begin
                push!(cs, a)
                push!(comps, [String(a), String(b.args[1])])
                arg = deepcopy(arg)
                b = deepcopy(arg.args[2])
                push!(b.args, Expr(:kw, :name, Meta.quot(a)))
                arg.args[2] = b
                push!(expr.args, arg)
            end
            _ => error("`@components` only takes assignment expressions. Got $arg")
        end
    end
    dict[:components] = comps
end

function parse_extend!(exprs, ext, dict, body)
    expr = Expr(:block)
    push!(exprs, expr)
    body = deepcopy(body)
    MLStyle.@match body begin
        Expr(:(=), a, b) => begin
            vars = nothing
            if Meta.isexpr(b, :(=))
                vars = a
                if !Meta.isexpr(vars, :tuple)
                    error("`@extend` destructuring only takes an tuple as LHS. Got $body")
                end
                a, b = b.args
                vars, a, b
            end
            ext[] = a
            push!(b.args, Expr(:kw, :name, Meta.quot(a)))
            dict[:extend] = [Symbol.(vars.args), a, b.args[1]]
            push!(expr.args, :($a = $b))
            if vars !== nothing
                push!(expr.args, :(@unpack $vars = $a))
            end
        end
        _ => error("`@extend` only takes an assignment expression. Got $body")
    end
end

function parse_variables!(exprs, vs, dict, mod, body, varclass, kwargs)
    expr = Expr(:block)
    push!(exprs, expr)
    for arg in body.args
        arg isa LineNumberNode && continue
        vv = parse_variable_def!(dict, mod, arg, varclass, kwargs)
        v = Num(vv)
        name = getname(v)
        push!(vs, name)
        push!(expr.args, :($name = $v))
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

function parse_icon!(icon, dict, mod, body::String)
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
    else
        error("$body is not a valid icon")
    end
end

function parse_icon!(icon, dict, mod, body::Expr)
    _icon = body.args[end]
    dict[:icon] = icon[] = MLStyle.@match _icon begin
        ::Symbol => get_var(mod, _icon)
        ::String => _icon
        Expr(:call, read, a...) => eval(_icon)
        _ => error("$_icon isn't a valid icon")
    end
end
