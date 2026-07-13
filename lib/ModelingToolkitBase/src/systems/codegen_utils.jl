const _COMPILE_OPTIONS = Dict{Symbol, Int}(
    :off => 0,
    :on => 1,
    :all => 2,
    :min => 3,
)

"""
    CompilerOptions(; optlevel = -1, compile = :default, infer = :default)

Options controlling the Julia compiler for generated functions.

Note that this feature is considered experimental.

# Fields
- `optlevel::Int`: LLVM optimization level (0-3), or -1 (default) to inherit from the module.
- `compile::Int`: Compilation mode as an integer (0=off, 1=on, 2=all, 3=min), or -1 (default)
  to inherit. Can also be specified as a symbol: `:off`, `:on`, `:all`, `:min`, or `:default`.
- `infer::Int`: Type inference (0=off, 1=on), or -1 (default) to inherit. Can also be specified
  as a Bool or `:default`.
"""
struct CompilerOptions
    optlevel::Int
    compile::Int
    infer::Int
end

function CompilerOptions(;
        optlevel::Int = -1, compile::Union{Int, Symbol} = :default,
        infer::Union{Int, Bool, Symbol} = :default
    )
    _compile = if compile isa Symbol
        compile === :default ? -1 : get(_COMPILE_OPTIONS, compile) do
                throw(ArgumentError("Invalid compile option: $compile. Valid options: :default, :off, :on, :all, :min"))
        end
    else
        compile
    end
    _infer = if infer isa Symbol
        infer === :default ? -1 :
            throw(ArgumentError("Invalid infer option: $infer. Valid options: :default, true, false"))
    elseif infer isa Bool
        Int(infer)
    else
        infer
    end
    return CompilerOptions(optlevel, _compile, _infer)
end

_has_options(co::CompilerOptions) = co.optlevel != -1 || co.compile != -1 || co.infer != -1

const _COMPILER_OPTIONS_SUPPORTED = isdefined(Base.Experimental, :set_compile!)

# Fallback modules with module-level compiler options for Julia versions
# that don't support per-method compiler options (Base.Experimental.set_compile!).
# We pre-create modules for the most common combinations:
#   - optimize=0
#   - optimize=1
#   - optimize=0, infer=false
# On Julia versions with per-method support, these modules are still defined
# but without any special compiler options (they are never used in that case).
module _EvalModuleOpt0
    @static if !isdefined(Base.Experimental, :set_compile!)
        Base.Experimental.@compiler_options optimize = 0
        using RuntimeGeneratedFunctions
        RuntimeGeneratedFunctions.init(@__MODULE__)
    end
end
module _EvalModuleOpt1
    @static if !isdefined(Base.Experimental, :set_compile!)
        Base.Experimental.@compiler_options optimize = 1
        using RuntimeGeneratedFunctions
        RuntimeGeneratedFunctions.init(@__MODULE__)
    end
end
module _EvalModuleOpt0NoInfer
    @static if !isdefined(Base.Experimental, :set_compile!)
        Base.Experimental.@compiler_options optimize = 0 infer = false
        using RuntimeGeneratedFunctions
        RuntimeGeneratedFunctions.init(@__MODULE__)
    end
end

"""
    _resolve_fallback_eval_module(co::CompilerOptions)

Return a pre-made module with the appropriate module-level compiler options for the
given `CompilerOptions`. This is used as a fallback on Julia versions that do not support
per-method compiler options. Only a limited set of option combinations is supported;
unsupported combinations will throw an error.
"""
function _resolve_fallback_eval_module(co::CompilerOptions)
    if co.compile != -1
        throw(
            ArgumentError(
                "Non-default `compile` compiler option is not supported on Julia $VERSION. " *
                    "Per-method compiler options require Julia with `Base.Experimental.set_compile!` support."
            )
        )
    end
    if co.optlevel == 0 && co.infer == -1
        return _EvalModuleOpt0
    elseif co.optlevel == 1 && co.infer == -1
        return _EvalModuleOpt1
    elseif co.optlevel == 0 && co.infer == 0
        return _EvalModuleOpt0NoInfer
    else
        throw(
            ArgumentError(
                "The compiler option combination (optlevel=$(co.optlevel), infer=$(co.infer)) " *
                    "is not supported on Julia $VERSION without per-method compiler options. " *
                    "Supported fallback combinations: (optlevel=0), (optlevel=1), (optlevel=0, infer=false)."
            )
        )
    end
end

"""
    $(TYPEDSIGNATURES)

Given a function expression `expr`, return a callable version of it.

# Keyword arguments
- `eval_expression`: Whether to use `eval` to make `expr` callable. If `false`, uses
  RuntimeGeneratedFunctions.jl.
- `eval_module`: The module to `eval` the expression `expr` in. If `!eval_expression`,
  this is the cache and context module for the `RuntimeGeneratedFunction`.
- `compiler_options`: A [`CompilerOptions`](@ref) controlling LLVM optimization level,
  compilation mode, and type inference for the generated function.
"""
function eval_or_rgf(
        expr::Expr; eval_expression = false, eval_module = @__MODULE__,
        compiler_options::CompilerOptions = CompilerOptions()
    )
    if _has_options(compiler_options)
        if _COMPILER_OPTIONS_SUPPORTED
            expr = _inject_compiler_options_meta(expr, compiler_options)
        else
            if eval_module !== @__MODULE__
                throw(
                    ArgumentError(
                        "Cannot use both a non-default `eval_module` and non-default `compiler_options` " *
                            "on Julia $VERSION. Per-method compiler options require `Base.Experimental.set_compile!` support."
                    )
                )
            end
            eval_module = _resolve_fallback_eval_module(compiler_options)
        end
    end
    if eval_expression
        return eval_module.eval(expr)
    else
        return drop_expr(RuntimeGeneratedFunction(eval_module, eval_module, expr))
    end
end

eval_or_rgf(::Nothing; kws...) = nothing

function _inject_compiler_options_meta(expr::Expr, co::CompilerOptions)
    if expr.head === :function || expr.head === :->
        body = expr.args[2]
        metas = Expr[]
        co.optlevel != -1 && push!(metas, Expr(:meta, Expr(:optlevel, co.optlevel)))
        co.compile != -1 && push!(metas, Expr(:meta, Expr(:compile, co.compile)))
        co.infer != -1 && push!(metas, Expr(:meta, Expr(:infer, co.infer)))
        if body isa Expr && body.head === :block
            new_body = Expr(:block, metas..., body.args...)
        else
            new_body = Expr(:block, metas..., body)
        end
        return Expr(expr.head, expr.args[1], new_body)
    end
    return expr
end

"""
    $(TYPEDSIGNATURES)

Return the name for the `i`th argument in a function generated by `build_function_wrapper`.
"""
function generated_argument_name(i::Int)
    return Symbol(:__mtk_arg_, i)
end

function compute_array_variable_buffer_idxs(@nospecialize(args); ignore_vars = Set{SymbolicT}(), ignore_arg_idxs = nothing)
    return _compute_array_variable_buffer_idxs(args isa Vector ? args : collect(args), ignore_vars, ignore_arg_idxs)
end

function _compute_array_variable_buffer_idxs(args::Vector, ignore_vars, ignore_arg_idxs = nothing)
    # map array symbolic to an identically sized array where each element is (buffer_idx, idx_in_buffer)
    var_to_arridxs = Dict{SymbolicT, Vector{Tuple{Int, Int}}}()
    for (i, arg) in enumerate(args)
        # filter out non-arrays
        # any element of args which is not an array is assumed to not contain a
        # scalarized array symbolic. This works because the only non-array element
        # is the independent variable
        arg isa Vector{SymbolicT} || continue
        # entire arg-vectors whose decomposition is already accounted for (e.g. the
        # parameter slice, handled via the cached `param_var_to_arridxs`) are skipped
        # here so we never re-run `split_indexed_var`/`get_stable_index` on them.
        ignore_arg_idxs === nothing || i in ignore_arg_idxs && continue

        # go through symbolics
        for (j, var) in enumerate(arg)
            var = unwrap(var)
            # filter out non-array-symbolics
            arrvar, isarr = split_indexed_var(var)
            isarr || continue
            arrvar in ignore_vars && continue
            # get and/or construct the buffer storing indexes
            idxbuffer = get(var_to_arridxs, arrvar, nothing)
            if idxbuffer === nothing
                idxbuffer = map(Returns((0, 0)), SU.stable_eachindex(arrvar))
            end
            idxbuffer[SU.as_linear_idx(SU.shape(arrvar), get_stable_index(var))] = (i, j)
            var_to_arridxs[arrvar] = idxbuffer
        end
    end

    return var_to_arridxs
end

function array_variable_buffer_idxs_to_assignments(
        var_to_arridxs::Dict{SymbolicT, Vector{Tuple{Int, Int}}};
        argument_name = generated_argument_name, buffer_offset = 0,
        filter_vars = nothing,
    )
    assignments = Assignment[]
    for (arrvar, idxs) in var_to_arridxs
        filter_vars === nothing || arrvar in filter_vars || continue
        # all elements of the array need to be present in `args` to form the
        # reconstructing assignment
        any(iszero ∘ first, idxs) && continue

        var_size_expr = Expr(:tuple)
        sharrvar = SU.shape(arrvar)
        for ax in sharrvar
            push!(var_size_expr.args, length(ax))
        end

        # if they are all in the same buffer, we can take a shortcut and `view` into it
        if allequal(Iterators.map(first, idxs))
            buffer_idx = first(first(idxs)) + buffer_offset
            idxs = map(last, idxs)
            shape = SU.ShapeVecT((1:length(idxs),))
            # if all the elements are contiguous and ordered, turn the array of indexes into a range
            # to help reduce allocations
            if first(idxs) < last(idxs) && vec(idxs) == first(idxs):last(idxs)
                expr = Expr(:call, view, argument_name(buffer_idx), first(idxs):last(idxs))
            elseif vec(idxs) == first(idxs):-1:last(idxs)
                expr = Expr(:call, view, argument_name(buffer_idx), first(idxs):-1:last(idxs))
            elseif length(idxs) <= 16
                expr = Expr(:call, SVector)
                append!(expr.args, idxs)
                expr = Expr(:call, view, argument_name(buffer_idx), expr)
            else
                expr = Expr(:call, view, argument_name(buffer_idx), idxs)
            end
        else
            if length(idxs) <= 16
                expr = Expr(:call, SVector)
            else
                expr = Expr(:vect)
            end
            for idx in idxs
                i, j = idx
                push!(expr.args, Expr(:ref, argument_name(i), j))
            end
        end
        expr = Expr(:call, reshape, expr, var_size_expr)
        if any(!isone ∘ first, sharrvar)
            expr_origin = Expr(:call, Origin)
            for ax in sharrvar
                push!(expr_origin.args, first(ax))
            end
            expr = Expr(:call, expr_origin, expr)
        end
        push!(assignments, Assignment(arrvar, expr))
    end

    return assignments
end

"""
    $(TYPEDSIGNATURES)

Given the arguments to `build_function_wrapper`, return a list of assignments which
reconstruct array variables if they are present scalarized in `args`.

# Keyword Arguments

- `argument_name` a function of the form `(::Int) -> Symbol` which takes the index of
  an argument to the generated function and returns the name of the argument in the
  generated function.
"""
function array_variable_assignments(
        @nospecialize(args...); ignore_vars = Set{SymbolicT}(), filter_vars = nothing,
        argument_name = generated_argument_name, buffer_offset = 0, ignore_arg_idxs = nothing
    )
    var_to_arridxs = compute_array_variable_buffer_idxs(args; ignore_vars, ignore_arg_idxs)
    return array_variable_buffer_idxs_to_assignments(var_to_arridxs; argument_name, buffer_offset, filter_vars)
end

"""
    $(TYPEDSIGNATURES)

Check if the variable `var` is a delayed variable, where `iv` is the independent
variable.
"""
function isdelay(var, iv)
    iv === nothing && return false
    if iscall(var) && ModelingToolkitBase.isoperator(var, Differential)
        return isdelay(arguments(var)[1], iv)
    end
    isvariable(var) || return false
    isparameter(var) && return false
    if iscall(var) && !ModelingToolkitBase.isoperator(var, Symbolics.Operator)
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

"""
The argument of generated functions corresponding to the history function.
"""
const DDE_HISTORY_FUN = SSym(:___history___; type = SU.FnType{Tuple{Any, <:Real}, Vector{Real}, Nothing}, shape = SU.Unknown(1))
const BVP_SOLUTION = SSym(:__sol__; type = Symbolics.FnType{Tuple{<:Real}, Vector{Real}, Nothing}, shape = SU.Unknown(1))
const DDE_AT_IDX_SYM = SSym(:__delayvar_idxₘₜₖ; type = Int, shape = UnitRange{Int}[])
const DDE_DELAY_SYM = SSym(:__delayxₘₜₖ; type = Real, shape = UnitRange{Int}[])

"""
    $(TYPEDSIGNATURES)

Turn delayed unknowns in `eqs` into calls to `DDE_HISTORY_FUNCTION`.

# Arguments

- `sys`: The system of DDEs.
- `eqs`: The equations to convert.

# Keyword Arguments

- `param_arg`: The name of the variable containing the parameter object.
- `histfn`: The history function to use for codegen, called as `histfn(p, t)`
"""
function delay_to_function(
        sys::AbstractSystem, eqs = full_equations(sys); param_arg = MTKPARAMETERS_ARG, histfn = DDE_HISTORY_FUN
    )
    return delay_to_function(
        eqs,
        get_iv(sys),
        Dict{Any, Int}(operation(s) => i for (i, s) in enumerate(unknowns(sys))),
        parameters(sys),
        histfn; param_arg
    )
end
function delay_to_function(eqs::Vector, iv, sts, ps, h; param_arg = MTKPARAMETERS_ARG)
    return delay_to_function.(eqs, (iv,), (sts,), (ps,), (h,); param_arg)
end
function delay_to_function(eq::Equation, iv, sts, ps, h; param_arg = MTKPARAMETERS_ARG)
    return delay_to_function(eq.lhs, iv, sts, ps, h; param_arg) ~ delay_to_function(
        eq.rhs, iv, sts, ps, h; param_arg
    )
end
function delay_to_function(expr, iv, sts, ps, h; param_arg = MTKPARAMETERS_ARG)
    if isdelay(expr, iv)
        v = operation(expr)
        time = arguments(expr)[1]
        idx = sts[v]
        return term(getindex, h(param_arg, time), idx, type = Real)
    elseif iscall(expr)
        return maketerm(
            typeof(expr),
            operation(expr),
            map(x -> delay_to_function(x, iv, sts, ps, h; param_arg), arguments(expr)),
            metadata(expr)
        )
    else
        return expr
    end
end

function __search_dervars_recurse(x::SymbolicT)
    return iscall(x) && Moshi.Match.@match x begin
        BSImpl.Term(; f) && if f isa Operator end => false
        _ => true
    end
end

struct ParameterArrayAssignments
    var_to_arridxs::Dict{SymbolicT, Vector{Tuple{Int, Int}}}
end

function should_invalidate_mutable_cache_entry(::Type{ParameterArrayAssignments}, patch::NamedTuple)
    return haskey(patch, :ps)
end

function find_arrvars_is_atomic(ex::SymbolicT)
    return SU.default_is_atomic(ex) && Symbolics.isarraysymbolic(ex)
end

"""
    $(TYPEDSIGNATURES)

A wrapper around `build_function` which performs the necessary transformations for
code generation of all types of systems. `expr` is the expression returned from the
generated functions, and `args` are the arguments.

# Keyword Arguments

- `p_start`, `p_end`: Denotes the indexes in `args` where the buffers of the splatted
  `MTKParameters` object are present. These are collapsed into a single argument and
  destructured inside the function. `p_start` must also be provided for non-split systems
  since it is used by `wrap_delays`.
- `u_arg`: The index in `args` of the argument corresponding to `unknowns(sys)` (the `u`
  vector). If `-1` (the default), the u vector is not treated specially. Otherwise, the
  argument must be a `Vector` and is wrapped in a `DestructuredArgs` with the common
  identifier `MTKUNKNOWNS_ARG`, giving it the predictable name `___mtkunknowns___`.
- `compress_args`: A list of argument ranges that end before `p_start`.
  Each range will be compressed into a single argument to the function. For example,
  If there are 5 elements in `args` and `compress_args = [2:3]`, then the generated function
  will take 4 arguments, where the second should be an indexable collection of the second
  and third elements in `args`.
- `wrap_delays`: Whether to transform delayed unknowns of `sys` present in `expr` into
  calls to a history function. The history function is added to the list of arguments
  right before parameters, at the index `p_start`.
- `histfn`: The history function to use for transforming delayed terms. For any delayed
  term `x(expr)`, this is called as `histfn(p, expr)` where `p` is the parameter object.
- `histfn_symbolic`: The symbolic history function variable to add as an argument to the
  generated function.
- `wrap_code`: Forwarded to `build_function`.
- `create_bindings`: Whether to explicitly destructure arrays of symbolics present in
  `args` in the generated code. If `false`, all usages of the individual symbolics will
  instead call `getindex` on the relevant argument. This is useful if the generated
  function writes to one of its arguments and expects subsequent code to use the new
  values. Note that the collapsed `MTKParameters` argument will always be explicitly
  destructured regardless of this keyword argument.
- `output_type`: The type of the output buffer. If `mkarray` (see below) is `nothing`,
  this will be passed to the `similarto` argument of `build_function`. If `output_type`
  is `Tuple`, `expr` will be wrapped in `SymbolicUtils.Code.MakeTuple` (regardless of
  whether it is scalar or an array).
- `mkarray`: A function which accepts `expr` and `output_type` and returns a code
  generation object similar to `MakeArray` or `MakeTuple` to be used to generate
  code for `expr`.
- `wrap_mtkparameters`: Whether to collapse parameter buffers for a split system into a
  argument.
- `extra_assignments`: Extra `Assignment` statements to prefix to `expr`, after all other
  assignments.

All other keyword arguments are forwarded to `build_function`.
"""
Base.@nospecializeinfer function build_function_wrapper(
        sys::AbstractSystem, @nospecialize(expr), @nospecialize(args...); p_start = 2,
        p_end = is_time_dependent(sys) ? length(args) - 1 : length(args), compress_args = UnitRange{Int}[],
        non_standard_param_layout = false, u_arg::Integer = -1,
        wrap_delays = is_dde(sys), histfn = DDE_HISTORY_FUN, histfn_symbolic = histfn, wrap_code = identity,
        create_bindings = false, @nospecialize(output_type::Union{Nothing, Type} = nothing), mkarray = nothing,
        wrap_mtkparameters = true, extra_assignments = Assignment[],
        n_param_buffers::Int = -1, optimize = nothing, @nospecialize(kwargs...)
    )
    isscalar = !(expr isa AbstractArray || symbolic_type(expr) == ArraySymbolic())
    args = Vector{Any}(collect(args))
    assignments = Assignment[]

    if u_arg != -1
        args[u_arg] isa AbstractVector ||
            throw(ArgumentError("argument at u_arg = $u_arg must be a Vector, got $(typeof(args[u_arg]))"))
    end

    u_argument_name = if u_arg == -1
        generated_argument_name
    else
        push!(assignments, Assignment(generated_argument_name(u_arg), :___mtkunknowns___))
        i -> i == u_arg ? :___mtkunknowns___ : generated_argument_name(i)
    end
    # turn delayed unknowns into calls to the history function
    if wrap_delays
        param_arg = is_split(sys) ? MTKPARAMETERS_ARG : generated_argument_name(p_start)
        delayfn_expr = STerm(
            Base.Fix1, SArgsT((histfn_symbolic, param_arg));
            type = SU.FnType{Tuple, Vector{Real}, Any}, shape = SU.Unknown(1)
        )
        for (i, var) in enumerate(unknowns(sys))
            Moshi.Match.@match var begin
                BSImpl.Term(; f) && if f isa SymbolicT end => begin
                    rhs = STerm(
                        ComposedFunction, SArgsT((Base.Fix2(getindex, i), delayfn_expr));
                        type = SU.FnType{Tuple, Real, Any}, shape = UnitRange{Int}[]
                    )
                    push!(assignments, Assignment(f, rhs))
                end
                _ => nothing
            end
        end

        # add extra argument
        insert!(args, p_start, histfn_symbolic)
        p_start += 1
        p_end += 1
    end

    ir_info = get_ir_info(sys)
    ir = get_irstructure(sys)
    expr = ir_info.obs_subber(expr)
    stmts = SU.IRStructureSearchBuffer(ir, Set{SymbolicT}())
    SU.search_variables!(stmts, expr; is_atomic = Returns(true))
    for subexpr in stmts
        SU.populate_ir!(ir, subexpr)
    end

    # Subset the `IRStructure` so that optimization passes run on the minimal set of
    # expressions
    ir = SU.subset_ir(ir, stmts)
    # We don't use `replace_node!` here to insert the observed, and instead prefer to
    # substitute as done above. The issue arises due to array variables passed scalarized
    # to the function. Consider a case where `p1` and `p2` are array parameters, `p1 => p2`
    # is a binding, and `p2` is passed scalarized to the function. Codegen for `Func` adds
    # `p2[1] => Symbol("...")` to `rewrites`, expecting that `SU.Code.codegen_ir!` checks
    # this and doesn't actually try to codegen the `getindex`. However, if an expression
    # has `p1[1]`, then `replace_node!(ir, p1, p2)` won't update this expression. As a result,
    # when codegen hits `p1[1]`, it won't use the replacement and codegen the `getindex`, which
    # generates `p2[1]`. `p2` unscalarized is not defined because `array_variable_assignments`
    # doesn't find it in the expression. Thus, where codegen should generate something like
    # `var"##cse#1234" = __mtk_arg_2[13]`, it will generate
    # ```julia
    # var"##cse#1234" = var"p2"
    # var"##cse#1235" = var"##cse#1234"[1]
    # ```

    required_arrvars = Set{SymbolicT}()
    search_buffer = SU.IRStructureSearchBuffer(ir, required_arrvars)
    SU.search_variables!(search_buffer, expr; is_atomic = find_arrvars_is_atomic, recurse = !SU.default_is_atomic)

    # assignments for reconstructing scalarized array symbolics
    if non_standard_param_layout
        append!(assignments, array_variable_assignments(args...; filter_vars = required_arrvars, argument_name = u_argument_name))
    else
        # `n_param_buffers` marks the boundary between the `reorder_parameters` buffers and
        # any `cachesyms` the caller appended to the parameter slice (cachesyms are always a
        # suffix). When it is `-1` (unspecified) the whole slice is parameters (no cachesym
        # tail). The cache always stores the PARAM-ONLY decomposition (shared, propagated to
        # SCC subsystems); the per-call cachesym tail is decomposed fresh and merged below.
        pbuf_end = n_param_buffers == -1 ? p_end : (p_start + n_param_buffers - 1)
        cached = check_mutable_cache(sys, ParameterArrayAssignments, ParameterArrayAssignments, nothing)
        if cached isa ParameterArrayAssignments
            param_var_to_arridxs = cached.var_to_arridxs
        else
            param_var_to_arridxs = compute_array_variable_buffer_idxs(args[p_start:pbuf_end])
            store_to_mutable_cache!(sys, ParameterArrayAssignments, ParameterArrayAssignments(param_var_to_arridxs))
        end
        # Merge the cachesym tail, shifting its slice-relative buffer indices past the
        # parameter buffers. cachesyms are prior-SCC unknowns / CSE temporaries -- never array
        # parameters -- so their array-variable keyset is disjoint from the params'. The
        # merged dict is identical to decomposing the whole `args[p_start:p_end]` slice.
        if pbuf_end < p_end
            merged = copy(param_var_to_arridxs)
            tail = compute_array_variable_buffer_idxs(
                args[(pbuf_end + 1):p_end]; ignore_vars = keys(merged)
            )
            shift = pbuf_end - p_start + 1
            for (arrvar, idxs) in tail
                merged[arrvar] = [(i + shift, j) for (i, j) in idxs]
            end
            param_var_to_arridxs = merged
        end
        append!(
            assignments, array_variable_buffer_idxs_to_assignments(
                param_var_to_arridxs; buffer_offset = p_start - 1, filter_vars = required_arrvars
            )
        )
        other_assigns = array_variable_assignments(args...; ignore_vars = keys(param_var_to_arridxs), filter_vars = required_arrvars, argument_name = u_argument_name, ignore_arg_idxs = p_start:p_end)
        append!(assignments, other_assigns)
    end
    append!(assignments, extra_assignments)

    for (i, arg) in enumerate(args)
        # Make sure to use the proper names for arguments
        args[i] = if u_arg != -1 && i == u_arg
            DestructuredArgs(arg, MTKUNKNOWNS_ARG; create_bindings)
        elseif symbolic_type(arg) == NotSymbolic() && arg isa AbstractArray
            DestructuredArgs(arg, generated_argument_name(i); create_bindings)
        else
            arg
        end
    end

    # wrap into a single MTKParameters argument
    if is_split(sys) && wrap_mtkparameters
        if p_start > p_end
            # In case there are no parameter buffers, still insert an argument
            args = [args[1:(p_start - 1)]; MTKPARAMETERS_ARG; args[(p_end + 1):end]]
        else
            # cannot apply `create_bindings` here since it doesn't nest
            args = [args[1:(p_start - 1)]; DestructuredArgs(collect(args[p_start:p_end]), MTKPARAMETERS_ARG); args[(p_end + 1):end]]
        end
    end

    sort!(compress_args; by = first)
    reverse!(compress_args)
    for (i, range) in enumerate(compress_args)
        compressed = DestructuredArgs(args[range], Symbol(:__compressed, i))
        deleteat!(args, range)
        insert!(args, first(range), compressed)
    end

    # add preface assignments
    if has_preface(sys) && (pref = preface(sys)) !== nothing
        append!(assignments, pref)
    end

    wrap_code = wrap_code .∘ wrap_assignments(false, assignments)

    # handling of `output_type` and `mkarray`
    similarto = nothing
    if output_type === Tuple
        expr = MakeTuple(Tuple(expr))
    elseif mkarray === nothing
        similarto = output_type
    else
        expr = mkarray(expr, output_type)
    end

    optimize = resolve_optimize_option(optimize)
    # Bundle the code-generation options into a single `Symbolics.CodegenFunctionOptions` rather than
    # splatting them (and any stray forwarded keyword arguments) as `kwargs...`. Because
    # `CodegenFunctionOptions` is one concrete type regardless of the option values, `codegen_function`
    # only needs to be compiled once instead of once per distinct set of keyword arguments.
    # Unknown keyword arguments are ignored by the `CodegenFunctionOptions` constructor, matching the
    # previous behaviour where `codegen_function` silently dropped them.
    codegen_options = Symbolics.CodegenFunctionOptions(; wrap_code, similarto, optimize, kwargs...)
    return Symbolics.codegen_function(ir, expr, args, codegen_options)
end

resolve_optimize_option(x) = x
resolve_optimize_option(::Nothing) = nothing

"""
    $(TYPEDEF)

A wrapper around a generated in-place and out-of-place function. The type-parameter `P`
must be a 3-tuple where the first element is the index of the parameter object in the
arguments, the second is the expected number of arguments in the out-of-place variant
of the function, and the third is a boolean indicating whether the generated functions
are for a split system. For scalar functions, the inplace variant can be `nothing`.
"""
struct GeneratedFunctionWrapper{P, O, I} <: Function
    f_oop::O
    f_iip::I
end

function GeneratedFunctionWrapper{P}(foop::O, fiip::I) where {P, O, I}
    return GeneratedFunctionWrapper{P, O, I}(foop, fiip)
end

function GeneratedFunctionWrapper{P}(::Type{Val{true}}, foop, fiip; kwargs...) where {P}
    return :($(GeneratedFunctionWrapper{P})($foop, $fiip))
end

function GeneratedFunctionWrapper{P}(
        ::Type{Val{false}}, foop, fiip;
        compiler_options::CompilerOptions = CompilerOptions(), kws...
    ) where {P}
    return GeneratedFunctionWrapper{P}(
        eval_or_rgf(foop; compiler_options, kws...),
        eval_or_rgf(fiip; compiler_options, kws...)
    )
end

function (gfw::GeneratedFunctionWrapper)(args...)
    return _generated_call(gfw, args...)
end

function SciMLBase.numargs(::GeneratedFunctionWrapper{P}) where {P}
    n_oop = P[2]
    return (n_oop, n_oop + 1)
end

@generated function _generated_call(gfw::GeneratedFunctionWrapper{P}, args...) where {P}
    paramidx, nargs, issplit = P
    iip = false
    # IIP case has one more argument
    if length(args) == nargs + 1
        nargs += 1
        paramidx += 1
        iip = true
    end
    if length(args) != nargs
        throw(ArgumentError("Expected $nargs arguments, got $(length(args))."))
    end

    # the function to use
    f = iip ? :(gfw.f_iip) : :(gfw.f_oop)
    # non-split systems just call it as-is
    if !issplit
        return :($f(args...))
    end
    if args[paramidx] <: Union{Tuple, MTKParameters} &&
            !(args[paramidx] <: Tuple{Vararg{Number}})
        # for split systems, call it as-is if the parameter object is a tuple or MTKParameters
        # but not if it is a tuple of numbers
        return :($f(args...))
    else
        # The user provided a single buffer/tuple for the parameter object, so wrap that
        # one in a tuple
        fargs = ntuple(Val(length(args))) do i
            i == paramidx ? :((args[$i], nothing)) : :(args[$i])
        end
        return :($f($(fargs...)))
    end
end

"""
    $(TYPEDSIGNATURES)

Optionally compile a method and optionally wrap it in a `GeneratedFunctionWrapper` on the
basis of `expression` `wrap_gfw`, both of type `Union{Type{Val{true}}, Type{Val{false}}}`.
`gfw_args` is the first type parameter of `GeneratedFunctionWrapper`. `f` is a tuple of
function expressions of the form `(oop, iip)` or a single out-of-place function expression.

# Keyword Arguments

- `compiler_options`: A [`CompilerOptions`](@ref) controlling LLVM optimization level,
  compilation mode, and type inference for the generated function.

All other keyword arguments are forwarded to `eval_or_rgf`.
"""
function maybe_compile_function(
        expression, wrap_gfw::Type{Val{true}},
        gfw_args::Tuple{Int, Int, Bool}, f::NTuple{2, Expr};
        compiler_options::CompilerOptions = CompilerOptions(), kwargs...
    )
    return GeneratedFunctionWrapper{gfw_args}(expression, f...; compiler_options, kwargs...)
end

function maybe_compile_function(
        expression::Type{Val{false}}, wrap_gfw::Type{Val{false}},
        gfw_args::Tuple{Int, Int, Bool}, f::NTuple{2, Expr};
        compiler_options::CompilerOptions = CompilerOptions(), kwargs...
    )
    return eval_or_rgf.(f; compiler_options, kwargs...)
end

function maybe_compile_function(
        expression::Type{Val{true}}, wrap_gfw::Type{Val{false}},
        gfw_args::Tuple{Int, Int, Bool}, f::Union{Expr, NTuple{2, Expr}};
        compiler_options::CompilerOptions = CompilerOptions(), kwargs...
    )
    return f
end

function maybe_compile_function(
        expression, wrap_gfw::Type{Val{true}},
        gfw_args::Tuple{Int, Int, Bool}, f::Expr;
        compiler_options::CompilerOptions = CompilerOptions(), kwargs...
    )
    return GeneratedFunctionWrapper{gfw_args}(expression, f, nothing; compiler_options, kwargs...)
end

function maybe_compile_function(
        expression::Type{Val{false}}, wrap_gfw::Type{Val{false}},
        gfw_args::Tuple{Int, Int, Bool}, f::Expr;
        compiler_options::CompilerOptions = CompilerOptions(), kwargs...
    )
    return eval_or_rgf(f; compiler_options, kwargs...)
end
