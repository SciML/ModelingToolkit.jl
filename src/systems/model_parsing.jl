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

Base.parentmodule(m::Model) = parentmodule(m.f)

for f in (:connector, :mtkmodel)
    isconnector = f == :connector ? true : false
    @eval begin
        macro $f(name::Symbol, body)
            esc($(:_model_macro)(__module__, name, body, $isconnector))
        end
    end
end

flatten_equations(eqs::Vector{Equation}, eq::Equation) = vcat(eqs, [eq])
flatten_equations(eq::Vector{Equation}, eqs::Vector{Equation}) = vcat(eq, eqs)
function flatten_equations(eqs::Vector{Union{Equation, Vector{Equation}}})
    foldl(flatten_equations, eqs; init = Equation[])
end

function _model_macro(mod, name, expr, isconnector)
    exprs = Expr(:block)
    dict = Dict{Symbol, Any}(
        :constants => Dict{Symbol, Dict}(),
        :defaults => Dict{Symbol, Any}(),
        :kwargs => Dict{Symbol, Dict}(),
        :structural_parameters => Dict{Symbol, Dict}()
    )
    comps = Union{Symbol, Expr}[]
    ext = Ref{Any}(nothing)
    eqs = Expr[]
    icon = Ref{Union{String, URI}}()
    ps, sps, vs, = [], [], []
    c_evts = []
    d_evts = []
    kwargs = OrderedCollections.OrderedSet()
    where_types = Union{Symbol, Expr}[]

    push!(exprs.args, :(variables = []))
    push!(exprs.args, :(parameters = []))
    push!(exprs.args, :(systems = ODESystem[]))
    push!(exprs.args, :(equations = Union{Equation, Vector{Equation}}[]))
    push!(exprs.args, :(defaults = Dict{Num, Union{Number, Symbol, Function}}()))

    Base.remove_linenums!(expr)
    for arg in expr.args
        if arg.head == :macrocall
            parse_model!(exprs.args, comps, ext, eqs, icon, vs, ps,
                sps, c_evts, d_evts, dict, mod, arg, kwargs, where_types)
        elseif arg.head == :block
            push!(exprs.args, arg)
        elseif arg.head == :if
            MLStyle.@match arg begin
                Expr(:if, condition, x) => begin
                    parse_conditional_model_statements(comps, dict, eqs, exprs, kwargs,
                        mod, ps, vs, where_types,
                        parse_top_level_branch(condition, x.args)...)
                end
                Expr(:if, condition, x, y) => begin
                    parse_conditional_model_statements(comps, dict, eqs, exprs, kwargs,
                        mod, ps, vs, where_types,
                        parse_top_level_branch(condition, x.args, y)...)
                end
                _ => error("Got an invalid argument: $arg")
            end
        elseif isconnector
            # Connectors can have variables listed without `@variables` prefix or
            # begin block.
            parse_variable_arg!(
                exprs.args, vs, dict, mod, arg, :variables, kwargs, where_types)
        else
            error("$arg is not valid syntax. Expected a macro call.")
        end
    end

    iv = get(dict, :independent_variable, nothing)
    if iv === nothing
        iv = dict[:independent_variable] = get_t(mod, :t)
    end

    push!(exprs.args, :(push!(equations, $(eqs...))))
    push!(exprs.args, :(push!(parameters, $(ps...))))
    push!(exprs.args, :(push!(systems, $(comps...))))
    push!(exprs.args, :(push!(variables, $(vs...))))

    gui_metadata = isassigned(icon) > 0 ? GUIMetadata(GlobalRef(mod, name), icon[]) :
                   GUIMetadata(GlobalRef(mod, name))

    @inline pop_structure_dict!.(
        Ref(dict), [:constants, :defaults, :kwargs, :structural_parameters])

    sys = :($ODESystem($(flatten_equations)(equations), $iv, variables, parameters;
        name, systems, gui_metadata = $gui_metadata, defaults))

    if ext[] === nothing
        push!(exprs.args, :(var"#___sys___" = $sys))
    else
        push!(exprs.args, :(var"#___sys___" = $extend($sys, $(ext[]))))
    end

    isconnector && push!(exprs.args,
        :($Setfield.@set!(var"#___sys___".connector_type=$connector_type(var"#___sys___"))))

    !isempty(c_evts) && push!(exprs.args,
        :($Setfield.@set!(var"#___sys___".continuous_events=$SymbolicContinuousCallback.([
            $(c_evts...)
        ]))))

    !isempty(d_evts) && push!(exprs.args,
        :($Setfield.@set!(var"#___sys___".discrete_events=$SymbolicDiscreteCallback.([
            $(d_evts...)
        ]))))

    f = if length(where_types) == 0
        :($(Symbol(:__, name, :__))(; name, $(kwargs...)) = $exprs)
    else
        f_with_where = Expr(:where)
        push!(f_with_where.args,
            :($(Symbol(:__, name, :__))(; name, $(kwargs...))), where_types...)
        :($f_with_where = $exprs)
    end

    :($name = $Model($f, $dict, $isconnector))
end

pop_structure_dict!(dict, key) = length(dict[key]) == 0 && pop!(dict, key)

function update_kwargs_and_metadata!(dict, kwargs, a, def, indices, type, var,
        varclass, where_types, meta)
    if indices isa Nothing
        if !isnothing(meta) && haskey(meta, VariableUnit)
            uvar = gensym()
            push!(where_types, uvar)
            push!(kwargs, Expr(:kw, :($a::Union{Nothing, $uvar}), nothing))
        else
            push!(kwargs, Expr(:kw, :($a::Union{Nothing, $type}), nothing))
        end
        dict[:kwargs][getname(var)] = Dict(:value => def, :type => type)
    else
        vartype = gensym(:T)
        push!(kwargs,
            Expr(:kw,
                Expr(:(::), a,
                    Expr(:curly, :Union, :Nothing, Expr(:curly, :AbstractArray, vartype))),
                nothing))
        if !isnothing(meta) && haskey(meta, VariableUnit)
            push!(where_types, vartype)
        else
            push!(where_types, :($vartype <: $type))
        end
        dict[:kwargs][getname(var)] = Dict(:value => def, :type => AbstractArray{type})
    end
    if dict[varclass] isa Vector
        dict[varclass][1][getname(var)][:type] = AbstractArray{type}
    else
        dict[varclass][getname(var)][:type] = type
    end
end

function parse_variable_def!(dict, mod, arg, varclass, kwargs, where_types;
        def = nothing, indices::Union{Vector{UnitRange{Int}}, Nothing} = nothing,
        type::Type = Real, meta = Dict{DataType, Expr}())
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
        (:dist, VariableDistribution)]

    arg isa LineNumberNode && return
    MLStyle.@match arg begin
        a::Symbol => begin
            var = generate_var!(dict, a, varclass; indices, type)
            update_kwargs_and_metadata!(dict, kwargs, a, def, indices, type, var,
                varclass, where_types, meta)
            return var, def, Dict()
        end
        Expr(:(::), a, type) => begin
            type = getfield(mod, type)
            parse_variable_def!(
                dict, mod, a, varclass, kwargs, where_types; def, type, meta)
        end
        Expr(:(::), Expr(:call, a, b), type) => begin
            type = getfield(mod, type)
            def = _type_check!(def, a, type, varclass)
            parse_variable_def!(
                dict, mod, a, varclass, kwargs, where_types; def, type, meta)
        end
        Expr(:call, a, b) => begin
            var = generate_var!(dict, a, b, varclass, mod; indices, type)
            update_kwargs_and_metadata!(dict, kwargs, a, def, indices, type, var,
                varclass, where_types, meta)
            return var, def, Dict()
        end
        Expr(:(=), a, b) => begin
            Base.remove_linenums!(b)
            def, meta = parse_default(mod, b)
            var, def, _ = parse_variable_def!(
                dict, mod, a, varclass, kwargs, where_types; def, type, meta)
            if dict[varclass] isa Vector
                dict[varclass][1][getname(var)][:default] = def
            else
                dict[varclass][getname(var)][:default] = def
            end
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
                var, metadata_with_exprs = set_var_metadata(var, meta)
                return var, def, metadata_with_exprs
            end
            return var, def, Dict()
        end
        Expr(:tuple, a, b) => begin
            meta = parse_metadata(mod, b)
            var, def, _ = parse_variable_def!(
                dict, mod, a, varclass, kwargs, where_types; type, meta)
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
                var, metadata_with_exprs = set_var_metadata(var, meta)
                return var, def, metadata_with_exprs
            end
            return var, def, Dict()
        end
        Expr(:ref, a, b...) => begin
            indices = map(i -> UnitRange(i.args[2], i.args[end]), b)
            parse_variable_def!(dict, mod, a, varclass, kwargs, where_types;
                def, indices, type, meta)
        end
        _ => error("$arg cannot be parsed")
    end
end

function generate_var(a, varclass;
        indices::Union{Vector{UnitRange{Int}}, Nothing} = nothing,
        type = Real)
    var = indices === nothing ? Symbolics.variable(a; T = type) :
          first(@variables $a[indices...]::type)
    if varclass == :parameters
        var = toparam(var)
    elseif varclass == :independent_variables
        var = toiv(var)
    end
    var
end

singular(sym) = last(string(sym)) == 's' ? Symbol(string(sym)[1:(end - 1)]) : sym

function check_name_uniqueness(dict, a, newvarclass)
    for varclass in [:variables, :parameters, :structural_parameters, :constants]
        dvarclass = get(dict, varclass, nothing)
        if dvarclass !== nothing && a in keys(dvarclass)
            error("Cannot create a $(singular(newvarclass)) `$(a)` because there is already a $(singular(varclass)) with that name")
        end
    end
end

function generate_var!(dict, a, varclass;
        indices::Union{Vector{UnitRange{Int}}, Nothing} = nothing,
        type = Real)
    check_name_uniqueness(dict, a, varclass)
    vd = get!(dict, varclass) do
        Dict{Symbol, Dict{Symbol, Any}}()
    end
    vd isa Vector && (vd = first(vd))
    vd[a] = Dict{Symbol, Any}()
    indices !== nothing && (vd[a][:size] = Tuple(lastindex.(indices)))
    generate_var(a, varclass; indices, type)
end

function generate_var!(dict, a, b, varclass, mod;
        indices::Union{Vector{UnitRange{Int}}, Nothing} = nothing,
        type = Real)
    iv = b == :t ? get_t(mod, b) : generate_var(b, :independent_variables)
    prev_iv = get!(dict, :independent_variable) do
        iv
    end
    @assert isequal(iv, prev_iv) "Multiple independent variables are used in the model"
    check_name_uniqueness(dict, a, varclass)
    vd = get!(dict, varclass) do
        Dict{Symbol, Dict{Symbol, Any}}()
    end
    vd isa Vector && (vd = first(vd))
    vd[a] = Dict{Symbol, Any}()
    var = if indices === nothing
        Symbolics.variable(a, T = SymbolicUtils.FnType{Tuple{Any}, type})(iv)
    else
        vd[a][:size] = Tuple(lastindex.(indices))
        first(@variables $a(iv)[indices...]::type)
    end
    if varclass == :parameters
        var = toparam(var)
    end
    var
end

# Use the `t` defined in the `mod`. When it is unavailable, generate a new `t` with a warning.
function get_t(mod, t)
    try
        get_var(mod, t)
    catch e
        if e isa UndefVarError
            @warn("Could not find a predefined `t` in `$mod`; generating a new one within this model.\nConsider defining it or importing `t` (or `t_nounits`, `t_unitful` as `t`) from ModelingToolkit.")
            variable(:t)
        else
            throw(e)
        end
    end
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
        Expr(:vect, x...) => begin
            (a, nothing)
        end
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

function _set_var_metadata!(metadata_with_exprs, a, m, v::Expr)
    push!(metadata_with_exprs, m => v)
    a
end
function _set_var_metadata!(metadata_with_exprs, a, m, v)
    wrap(set_scalar_metadata(unwrap(a), m, v))
end

function set_var_metadata(a, ms)
    metadata_with_exprs = Dict{DataType, Expr}()
    for (m, v) in ms
        if m == VariableGuess && v isa Symbol
            v = quote
                $v
            end
        end
        a = _set_var_metadata!(metadata_with_exprs, a, m, v)
    end
    a, metadata_with_exprs
end

function get_var(mod::Module, b)
    if b isa Symbol
        isdefined(mod, b) && return getproperty(mod, b)
        isdefined(@__MODULE__, b) && return getproperty(@__MODULE__, b)
    end
    b
end

function parse_model!(exprs, comps, ext, eqs, icon, vs, ps, sps, c_evts, d_evts,
        dict, mod, arg, kwargs, where_types)
    mname = arg.args[1]
    body = arg.args[end]
    if mname == Symbol("@components")
        parse_components!(exprs, comps, dict, body, kwargs)
    elseif mname == Symbol("@extend")
        parse_extend!(exprs, ext, dict, mod, body, kwargs)
    elseif mname == Symbol("@variables")
        parse_variables!(exprs, vs, dict, mod, body, :variables, kwargs, where_types)
    elseif mname == Symbol("@parameters")
        parse_variables!(exprs, ps, dict, mod, body, :parameters, kwargs, where_types)
    elseif mname == Symbol("@structural_parameters")
        parse_structural_parameters!(exprs, sps, dict, mod, body, kwargs)
    elseif mname == Symbol("@equations")
        parse_equations!(exprs, eqs, dict, body)
    elseif mname == Symbol("@constants")
        parse_constants!(exprs, dict, body, mod)
    elseif mname == Symbol("@continuous_events")
        parse_continuous_events!(c_evts, dict, body)
    elseif mname == Symbol("@discrete_events")
        parse_discrete_events!(d_evts, dict, body)
    elseif mname == Symbol("@icon")
        isassigned(icon) && error("This model has more than one icon.")
        parse_icon!(body, dict, icon, mod)
    elseif mname == Symbol("@defaults")
        parse_system_defaults!(exprs, arg, dict)
    else
        error("$mname is not handled.")
    end
end

function parse_constants!(exprs, dict, body, mod)
    Base.remove_linenums!(body)
    for arg in body.args
        MLStyle.@match arg begin
            Expr(:(=), Expr(:(::), a, type), Expr(:tuple, b, metadata)) || Expr(:(=), Expr(:(::), a, type), b) => begin
                type = getfield(mod, type)
                b = _type_check!(get_var(mod, b), a, type, :constants)
                push!(exprs,
                    :($(Symbolics._parse_vars(
                        :constants, type, [:($a = $b), metadata], toconstant))))
                dict[:constants][a] = Dict(:value => b, :type => type)
                if @isdefined metadata
                    for data in metadata.args
                        dict[:constants][a][data.args[1]] = data.args[2]
                    end
                end
            end
            Expr(:(=), a, Expr(:tuple, b, metadata)) => begin
                push!(exprs,
                    :($(Symbolics._parse_vars(
                        :constants, Real, [:($a = $b), metadata], toconstant))))
                dict[:constants][a] = Dict{Symbol, Any}(:value => get_var(mod, b))
                for data in metadata.args
                    dict[:constants][a][data.args[1]] = data.args[2]
                end
            end
            Expr(:(=), a, b) => begin
                push!(exprs,
                    :($(Symbolics._parse_vars(
                        :constants, Real, [:($a = $b)], toconstant))))
                dict[:constants][a] = Dict(:value => get_var(mod, b))
            end
            _ => error("""Malformed constant definition `$arg`. Please use the following syntax:
                ```
                @constants begin
                    var = value, [description = "This is an example constant."]
                end
                ```
            """)
        end
    end
end

push_additional_defaults!(dict, a, b::Number) = dict[:defaults][a] = b
push_additional_defaults!(dict, a, b::QuoteNode) = dict[:defaults][a] = b.value
function push_additional_defaults!(dict, a, b::Expr)
    dict[:defaults][a] = readable_code(b)
end

function parse_system_defaults!(exprs, defaults_body, dict)
    for default_arg in defaults_body.args[end].args
        # for arg in default_arg.args
        MLStyle.@match default_arg begin
            # For cases like `p => 1` and `p => f()`. In both cases the definitions of
            # `a`, here `p` and when `b` is a function, here `f` are available while
            # defining the model
            Expr(:call, :(=>), a, b) => begin
                push!(exprs, :(defaults[$a] = $b))
                push_additional_defaults!(dict, a, b)
            end
            _ => error("Invalid `@defaults` entry: `$default_arg`")
        end
    end
end

function parse_structural_parameters!(exprs, sps, dict, mod, body, kwargs)
    Base.remove_linenums!(body)
    for arg in body.args
        MLStyle.@match arg begin
            Expr(:(=), Expr(:(::), a, type), b) => begin
                type = getfield(mod, type)
                b = _type_check!(get_var(mod, b), a, type, :structural_parameters)
                push!(sps, a)
                push!(kwargs, Expr(:kw, Expr(:(::), a, type), b))
                dict[:structural_parameters][a] = dict[:kwargs][a] = Dict(
                    :value => b, :type => type)
            end
            Expr(:(=), a, b) => begin
                push!(sps, a)
                push!(kwargs, Expr(:kw, a, b))
                dict[:structural_parameters][a] = dict[:kwargs][a] = Dict(:value => b)
            end
            a => begin
                push!(sps, a)
                push!(kwargs, a)
                dict[:structural_parameters][a] = dict[:kwargs][a] = Dict(:value => nothing)
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
                dict[:kwargs][x] = Dict(:value => nothing)
            end
            Expr(:kw, x) => begin
                push!(kwargs, Expr(:kw, x, nothing))
                dict[:kwargs][x] = Dict(:value => nothing)
            end
            Expr(:kw, x, y) => begin
                b.args[i] = Expr(:kw, x, x)
                push!(varexpr.args, :($x = $x === nothing ? $y : $x))
                push!(kwargs, Expr(:kw, x, nothing))
                dict[:kwargs][x] = Dict(:value => nothing)
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

function parse_variable_arg!(exprs, vs, dict, mod, arg, varclass, kwargs, where_types)
    name, ex = parse_variable_arg(dict, mod, arg, varclass, kwargs, where_types)
    push!(vs, name)
    push!(exprs, ex)
end

function convert_units(varunits::DynamicQuantities.Quantity, value)
    DynamicQuantities.ustrip(DynamicQuantities.uconvert(
        DynamicQuantities.SymbolicUnits.as_quantity(varunits), value))
end

function convert_units(
        varunits::DynamicQuantities.Quantity, value::AbstractArray{T}) where {T}
    DynamicQuantities.ustrip.(DynamicQuantities.uconvert.(
        DynamicQuantities.SymbolicUnits.as_quantity(varunits), value))
end

function convert_units(varunits::Unitful.FreeUnits, value)
    Unitful.ustrip(varunits, value)
end

function convert_units(varunits::Unitful.FreeUnits, value::AbstractArray{T}) where {T}
    Unitful.ustrip.(varunits, value)
end

function convert_units(varunits::Unitful.FreeUnits, value::Num)
    value
end

function convert_units(varunits::DynamicQuantities.Quantity, value::Num)
    value
end

function parse_variable_arg(dict, mod, arg, varclass, kwargs, where_types)
    vv, def, metadata_with_exprs = parse_variable_def!(
        dict, mod, arg, varclass, kwargs, where_types)
    name = getname(vv)

    varexpr = if haskey(metadata_with_exprs, VariableUnit)
        unit = metadata_with_exprs[VariableUnit]
        quote
            $name = if $name === nothing
                $setdefault($vv, $def)
            else
                try
                    $setdefault($vv, $convert_units($unit, $name))
                catch e
                    if isa(e, $(DynamicQuantities.DimensionError)) ||
                       isa(e, $(Unitful.DimensionError))
                        error("Unable to convert units for \'" * string(:($$vv)) * "\'")
                    elseif isa(e, MethodError)
                        error("No or invalid units provided for \'" * string(:($$vv)) *
                              "\'")
                    else
                        rethrow(e)
                    end
                end
            end
        end
    else
        quote
            $name = if $name === nothing
                $setdefault($vv, $def)
            else
                $setdefault($vv, $name)
            end
        end
    end

    metadata_expr = Expr(:block)
    for (k, v) in metadata_with_exprs
        push!(metadata_expr.args,
            :($name = $wrap($set_scalar_metadata($unwrap($name), $k, $v))))
    end

    push!(varexpr.args, metadata_expr)
    return name, varexpr
end

function handle_conditional_vars!(
        arg, conditional_branch, mod, varclass, kwargs, where_types)
    conditional_dict = Dict(:kwargs => Dict(),
        :parameters => Any[Dict{Symbol, Dict{Symbol, Any}}()],
        :variables => Any[Dict{Symbol, Dict{Symbol, Any}}()])
    for _arg in arg.args
        name, ex = parse_variable_arg(
            conditional_dict, mod, _arg, varclass, kwargs, where_types)
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

function parse_variables!(exprs, vs, dict, mod, body, varclass, kwargs, where_types)
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
                    kwargs,
                    where_types)
                push!(expr.args, conditional_expr)
                push_conditional_dict!(dict, condition, conditional_dict, nothing, varclass)
            end
            Expr(:if, condition, x, y) => begin
                conditional_expr = Expr(:if, condition, Expr(:block))
                conditional_dict = handle_conditional_vars!(x,
                    conditional_expr.args[2],
                    mod,
                    varclass,
                    kwargs,
                    where_types)
                conditional_y_expr, conditional_y_tuple = handle_y_vars(y,
                    conditional_dict,
                    mod,
                    varclass,
                    kwargs, where_types)
                push!(conditional_expr.args, conditional_y_expr)
                push!(expr.args, conditional_expr)
                push_conditional_dict!(dict,
                    condition,
                    conditional_dict,
                    conditional_y_tuple,
                    varclass)
            end
            _ => parse_variable_arg!(
                exprs, vs, dict, mod, arg, varclass, kwargs, where_types)
        end
    end
end

function handle_y_vars(y, dict, mod, varclass, kwargs, where_types)
    conditional_dict = if Meta.isexpr(y, :elseif)
        conditional_y_expr = Expr(:elseif, y.args[1], Expr(:block))
        conditional_dict = handle_conditional_vars!(y.args[2],
            conditional_y_expr.args[2],
            mod,
            varclass,
            kwargs,
            where_types)
        _y_expr, _conditional_dict = handle_y_vars(
            y.args[end], dict, mod, varclass, kwargs, where_types)
        push!(conditional_y_expr.args, _y_expr)
        (:elseif, y.args[1], conditional_dict, _conditional_dict)
    else
        conditional_y_expr = Expr(:block)
        handle_conditional_vars!(y, conditional_y_expr, mod, varclass, kwargs, where_types)
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

function parse_continuous_events!(c_evts, dict, body)
    dict[:continuous_events] = []
    Base.remove_linenums!(body)
    for arg in body.args
        push!(c_evts, arg)
        push!(dict[:continuous_events], readable_code.(c_evts)...)
    end
end

function parse_discrete_events!(d_evts, dict, body)
    dict[:discrete_events] = []
    Base.remove_linenums!(body)
    for arg in body.args
        push!(d_evts, arg)
        push!(dict[:discrete_events], readable_code.(d_evts)...)
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

function component_args!(a, b, varexpr, kwargs; index_name = nothing)
    # Whenever `b` is a function call, skip the first arg aka the function name.
    # Whenever it is a kwargs list, include it.
    start = b.head == :call ? 2 : 1
    for i in start:lastindex(b.args)
        arg = b.args[i]
        arg isa LineNumberNode && continue
        MLStyle.@match arg begin
            x::Symbol || Expr(:kw, x) => begin
                varname, _varname = _rename(a, x)
                b.args[i] = Expr(:kw, x, _varname)
                push!(varexpr.args, :((if $varname !== nothing
                    $_varname = $varname
                elseif @isdefined $x
                    # Allow users to define a var in `structural_parameters` and set
                    # that as positional arg of subcomponents; it is useful for cases
                    # where it needs to be passed to multiple subcomponents.
                    $_varname = $x
                end)))
                push!(kwargs, Expr(:kw, varname, nothing))
                # dict[:kwargs][varname] = nothing
            end
            Expr(:parameters, x...) => begin
                component_args!(a, arg, varexpr, kwargs)
            end
            Expr(:kw, x, y) => begin
                varname, _varname = _rename(a, x)
                b.args[i] = Expr(:kw, x, _varname)
                if isnothing(index_name)
                    push!(varexpr.args, :($_varname = $varname === nothing ? $y : $varname))
                else
                    push!(varexpr.args,
                        :($_varname = $varname === nothing ? $y : $varname[$index_name]))
                end
                push!(kwargs, Expr(:kw, varname, nothing))
                # dict[:kwargs][varname] = nothing
            end
            _ => error("Could not parse $arg of component $a")
        end
    end
end

model_name(name, range) = Symbol.(name, :_, collect(range))

function _parse_components!(body, kwargs)
    local expr
    varexpr = Expr(:block)
    comps = Vector{Union{Union{Expr, Symbol}, Expr}}[]
    comp_names = []

    Base.remove_linenums!(body)
    arg = body.args[end]

    MLStyle.@match arg begin
        Expr(:(=), a, Expr(:comprehension, Expr(:generator, b, Expr(:(=), c, d)))) => begin
            array_varexpr = Expr(:block)

            push!(comp_names, :($a...))
            push!(comps, [a, b.args[1], d])
            b = deepcopy(b)

            component_args!(a, b, array_varexpr, kwargs; index_name = c)

            expr = _named_idxs(a, d, :($c -> $b); extra_args = array_varexpr)
        end
        Expr(:(=), a, Expr(:comprehension, Expr(:generator, b, Expr(:filter, e, Expr(:(=), c, d))))) => begin
            error("List comprehensions with conditional statements aren't supported.")
        end
        Expr(:(=), a, Expr(:comprehension, Expr(:generator, b, Expr(:(=), c, d), e...))) => begin
            # Note that `e` is of the form `Tuple{Expr(:(=), c, d)}`
            error("More than one index isn't supported while building component array")
        end
        Expr(:block) => begin
            # TODO: Do we need this?
            error("Multiple `@components` block detected within a single block")
        end
        Expr(:(=), a, Expr(:for, Expr(:(=), c, d), b)) => begin
            Base.remove_linenums!(b)
            array_varexpr = Expr(:block)
            push!(array_varexpr.args, b.args[1:(end - 1)]...)
            push!(comp_names, :($a...))
            push!(comps, [a, b.args[end].args[1], d])
            b = deepcopy(b)

            component_args!(a, b.args[end], array_varexpr, kwargs; index_name = c)

            expr = _named_idxs(a, d, :($c -> $(b.args[end])); extra_args = array_varexpr)
        end
        Expr(:(=), a, b) => begin
            arg = deepcopy(arg)
            b = deepcopy(arg.args[2])

            component_args!(a, b, varexpr, kwargs)

            arg.args[2] = b
            expr = :(@named $arg)
            push!(comp_names, a)
            if (isa(b.args[1], Symbol) || Meta.isexpr(b.args[1], :.))
                push!(comps, [a, b.args[1]])
            end
        end
        _ => error("Couldn't parse the component body: $arg")
    end

    return comp_names, comps, expr, varexpr
end

function push_conditional_component!(ifexpr, expr_vec, comp_names, varexpr)
    blk = Expr(:block)
    push!(blk.args, varexpr)
    push!(blk.args, expr_vec)
    push!(blk.args, :($push!(systems, $(comp_names...))))
    push!(ifexpr.args, blk)
end

function handle_if_x!(mod, exprs, ifexpr, x, kwargs, condition = nothing)
    push!(ifexpr.args, condition)
    comp_names, comps, expr_vec, varexpr = _parse_components!(x, kwargs)
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
        comp_names, comps, expr_vec, varexpr = _parse_components!(y, kwargs)
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
            # Either the arg is top level component declaration or an invalid cause - both are handled by `_parse_components`
            _ => begin
                comp_names, comps, expr_vec, varexpr = _parse_components!(:(begin
                        $arg
                    end),
                    kwargs)
                push!(cs, comp_names...)
                push!(dict[:components], comps...)
                push!(exprs, varexpr, expr_vec)
            end
        end
    end
end

function _rename(compname, varname)
    compname = Symbol(compname, :__, varname)
    (compname, Symbol(:_, compname))
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
        ps, vs, where_types, component_blk, equations_blk, parameter_blk, variable_blk)
    parameter_blk !== nothing &&
        parse_variables!(
            exprs.args, ps, dict, mod, :(begin
                $parameter_blk
            end), :parameters, kwargs, where_types)

    variable_blk !== nothing &&
        parse_variables!(
            exprs.args, vs, dict, mod, :(begin
                $variable_blk
            end), :variables, kwargs, where_types)

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

function _type_check!(val, a, type, class)
    if val isa type
        return val
    else
        try
            return convert(type, val)
        catch e
            throw(TypeError(Symbol("`@mtkmodel`"),
                "`$class`, while assigning to `$a`", type, typeof(val)))
        end
    end
end
