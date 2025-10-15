""""""
function eval_or_rgf(expr::Expr; eval_expression = false, eval_module = @__MODULE__)
    if eval_expression
        return eval_module.eval(expr)
    else
        return drop_expr(RuntimeGeneratedFunction(eval_module, eval_module, expr))
    end
end
eval_or_rgf(::Nothing; kws...) = nothing
""""""
function generated_argument_name(i::Int)
    return Symbol(:__mtk_arg_, i)
end
""""""
function array_variable_assignments(args...; argument_name = generated_argument_name)
    var_to_arridxs = Dict{BasicSymbolic, Array{Tuple{Int, Int}}}()
    for (i, arg) in enumerate(args)
        symbolic_type(arg) == NotSymbolic() || continue
        arg isa AbstractArray || continue
        for (j, var) in enumerate(arg)
            var = unwrap(var)
            iscall(var) || continue
            operation(var) == getindex || continue
            arrvar = arguments(var)[1]
            idxbuffer = get!(
                () -> map(Returns((0, 0)), eachindex(arrvar)), var_to_arridxs, arrvar)
            Origin(first.(axes(arrvar))...)(idxbuffer)[arguments(var)[2:end]...] = (i, j)
        end
    end
    assignments = Assignment[]
    for (arrvar, idxs) in var_to_arridxs
        any(iszero ∘ first, idxs) && continue
        if allequal(Iterators.map(first, idxs))
            buffer_idx = first(first(idxs))
            idxs = map(last, idxs)
            if first(idxs) < last(idxs) && vec(idxs) == first(idxs):last(idxs)
                idxs = first(idxs):last(idxs)
            elseif vec(idxs) == last(idxs):-1:first(idxs)
                idxs = last(idxs):-1:first(idxs)
            else
                idxs = SArray{Tuple{size(idxs)...}}(idxs)
            end
            expr = term(reshape, term(view, argument_name(buffer_idx), idxs),
                size(arrvar))
        else
            elems = map(idxs) do idx
                i, j = idx
                term(getindex, argument_name(i), j)
            end
            expr = term(SymbolicUtils.Code.create_array, SArray, nothing,
                Val(ndims(arrvar)), Val(length(arrvar)), elems...)
        end
        if any(x -> !isone(first(x)), axes(arrvar))
            expr = term(Origin(first.(axes(arrvar))...), expr)
        end
        push!(assignments, arrvar ← expr)
    end
    return assignments
end
""""""
function isdelay(var, iv)
    iv === nothing && return false
    if iscall(var) && ModelingToolkit.isoperator(var, Differential)
        return isdelay(arguments(var)[1], iv)
    end
    isvariable(var) || return false
    isparameter(var) && return false
    if iscall(var) && !ModelingToolkit.isoperator(var, Symbolics.Operator)
        args = arguments(var)
        length(args) == 1 || return false
        arg = args[1]
        isequal(arg, iv) && return false
        iscall(arg) || return true
        issym(operation(arg)) && !iscalledparameter(arg) && return false
        return true
    end
    return false
end
""""""
const DDE_HISTORY_FUN = SSym(:___history___; type = SU.FnType{Tuple{Any, <:Real}, Vector{Real}}, shape = SU.Unknown(1))
const BVP_SOLUTION = SSym(:__sol__; type = Symbolics.FnType{Tuple{<:Real}, Vector{Real}}, shape = SU.Unknown(1))
""""""
function delay_to_function(
        sys::AbstractSystem, eqs = full_equations(sys); param_arg = MTKPARAMETERS_ARG, histfn = DDE_HISTORY_FUN)
    delay_to_function(eqs,
        get_iv(sys),
        Dict{Any, Int}(operation(s) => i for (i, s) in enumerate(unknowns(sys))),
        parameters(sys),
        histfn; param_arg)
end
function delay_to_function(eqs::Vector, iv, sts, ps, h; param_arg = MTKPARAMETERS_ARG)
    delay_to_function.(eqs, (iv,), (sts,), (ps,), (h,); param_arg)
end
function delay_to_function(eq::Equation, iv, sts, ps, h; param_arg = MTKPARAMETERS_ARG)
    delay_to_function(eq.lhs, iv, sts, ps, h; param_arg) ~ delay_to_function(
        eq.rhs, iv, sts, ps, h; param_arg)
end
function delay_to_function(expr, iv, sts, ps, h; param_arg = MTKPARAMETERS_ARG)
    if isdelay(expr, iv)
        v = operation(expr)
        time = arguments(expr)[1]
        idx = sts[v]
        return term(getindex, h(param_arg, time), idx, type = Real)
    elseif iscall(expr)
        return maketerm(typeof(expr),
            operation(expr),
            map(x -> delay_to_function(x, iv, sts, ps, h; param_arg), arguments(expr)),
            metadata(expr))
    else
        return expr
    end
end
""""""
function build_function_wrapper(sys::AbstractSystem, expr, args...; p_start = 2,
        p_end = is_time_dependent(sys) ? length(args) - 1 : length(args),
        wrap_delays = is_dde(sys), histfn = DDE_HISTORY_FUN, histfn_symbolic = histfn, wrap_code = identity,
        add_observed = true, filter_observed = Returns(true),
        create_bindings = false, output_type = nothing, mkarray = nothing,
        wrap_mtkparameters = true, extra_assignments = Assignment[], cse = true, kwargs...)
    isscalar = !(expr isa AbstractArray || symbolic_type(expr) == ArraySymbolic())
    obs = filter(filter_observed, observed(sys))
    if wrap_delays
        param_arg = is_split(sys) ? MTKPARAMETERS_ARG : generated_argument_name(p_start)
        obs = map(obs) do eq
            delay_to_function(sys, eq; param_arg, histfn)
        end
        expr = delay_to_function(sys, expr; param_arg, histfn)
        args = (args[1:(p_start - 1)]..., histfn_symbolic, args[p_start:end]...)
        p_start += 1
        p_end += 1
    end
    pdeps = get_parameter_dependencies(sys)
    if add_observed && !isempty(obs)
        obsidxs = observed_equations_used_by(sys, expr; obs)
    else
        obsidxs = Int[]
    end
    pdepidxs = observed_equations_used_by(sys, expr; obs = pdeps)
    for i in obsidxs
        union!(pdepidxs, observed_equations_used_by(sys, obs[i].rhs; obs = pdeps))
    end
    assignments = array_variable_assignments(args...)
    for eq in Iterators.flatten((pdeps[pdepidxs], obs[obsidxs]))
        push!(assignments, eq.lhs ← eq.rhs)
    end
    append!(assignments, extra_assignments)
    args = ntuple(Val(length(args))) do i
        arg = args[i]
        if symbolic_type(arg) == NotSymbolic() && arg isa AbstractArray
            DestructuredArgs(arg, generated_argument_name(i); create_bindings)
        else
            arg
        end
    end
    if is_split(sys) && wrap_mtkparameters
        if p_start > p_end
            args = (args[1:(p_start - 1)]..., MTKPARAMETERS_ARG, args[(p_end + 1):end]...)
        else
            args = (args[1:(p_start - 1)]...,
                DestructuredArgs(collect(args[p_start:p_end]), MTKPARAMETERS_ARG),
                args[(p_end + 1):end]...)
        end
    end
    if has_preface(sys) && (pref = preface(sys)) !== nothing
        append!(assignments, pref)
    end
    wrap_code = wrap_code .∘ wrap_assignments(isscalar, assignments)
    similarto = nothing
    if output_type === Tuple
        expr = MakeTuple(Tuple(expr))
        wrap_code = wrap_code[1]
    elseif mkarray === nothing
        similarto = output_type
    else
        expr = mkarray(expr, output_type)
        wrap_code = wrap_code[2]
    end
    if wrap_code isa Tuple && symbolic_type(expr) == ScalarSymbolic()
        wrap_code = wrap_code[1]
    end
    return build_function(expr, args...; wrap_code, similarto, cse, kwargs...)
end
""""""
struct GeneratedFunctionWrapper{P, O, I} <: Function
    f_oop::O
    f_iip::I
end
function GeneratedFunctionWrapper{P}(foop::O, fiip::I) where {P, O, I}
    GeneratedFunctionWrapper{P, O, I}(foop, fiip)
end
function GeneratedFunctionWrapper{P}(::Type{Val{true}}, foop, fiip; kwargs...) where {P}
    :($(GeneratedFunctionWrapper{P})($foop, $fiip))
end
function GeneratedFunctionWrapper{P}(::Type{Val{false}}, foop, fiip; kws...) where {P}
    GeneratedFunctionWrapper{P}(eval_or_rgf(foop; kws...), eval_or_rgf(fiip; kws...))
end
function (gfw::GeneratedFunctionWrapper)(args...)
    _generated_call(gfw, args...)
end
@generated function _generated_call(gfw::GeneratedFunctionWrapper{P}, args...) where {P}
    paramidx, nargs, issplit = P
    iip = false
    if length(args) == nargs + 1
        nargs += 1
        paramidx += 1
        iip = true
    end
    if length(args) != nargs
        throw(ArgumentError("Expected $nargs arguments, got $(length(args))."))
    end
    f = iip ? :(gfw.f_iip) : :(gfw.f_oop)
    if !issplit
        return :($f(args...))
    end
    if args[paramidx] <: Union{Tuple, MTKParameters} &&
       !(args[paramidx] <: Tuple{Vararg{Number}})
        return :($f(args...))
    else
        fargs = ntuple(Val(length(args))) do i
            i == paramidx ? :((args[$i], nothing)) : :(args[$i])
        end
        return :($f($(fargs...)))
    end
end
""""""
function maybe_compile_function(expression, wrap_gfw::Type{Val{true}},
        gfw_args::Tuple{Int, Int, Bool}, f::NTuple{2, Expr}; kwargs...)
    GeneratedFunctionWrapper{gfw_args}(expression, f...; kwargs...)
end
function maybe_compile_function(expression::Type{Val{false}}, wrap_gfw::Type{Val{false}},
        gfw_args::Tuple{Int, Int, Bool}, f::NTuple{2, Expr}; kwargs...)
    eval_or_rgf.(f; kwargs...)
end
function maybe_compile_function(expression::Type{Val{true}}, wrap_gfw::Type{Val{false}},
        gfw_args::Tuple{Int, Int, Bool}, f::Union{Expr, NTuple{2, Expr}}; kwargs...)
    return f
end
function maybe_compile_function(expression, wrap_gfw::Type{Val{true}},
        gfw_args::Tuple{Int, Int, Bool}, f::Expr; kwargs...)
    GeneratedFunctionWrapper{gfw_args}(expression, f, nothing; kwargs...)
end
function maybe_compile_function(expression::Type{Val{false}}, wrap_gfw::Type{Val{false}},
        gfw_args::Tuple{Int, Int, Bool}, f::Expr; kwargs...)
    eval_or_rgf(f; kwargs...)
end
