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

# Normalize an `expression`/`wrap_gfw`-style option (`Val{true}`/`Val{false}` as a type or
# instance, or a `Bool`) to a `Bool`.
_gfo_bool(::Type{Val{B}}) where {B} = B::Bool
_gfo_bool(::Val{B}) where {B} = B::Bool
_gfo_bool(b::Bool) = b

"""
    GeneratedFunctionOptions{expression, wrap_gfw}(; kwargs...)

Options for the code-generation entry points (`generate_rhs`, `generate_jacobian`, ...):
the "output/compile" layer sitting one level above [`BuildFunctionWrapperOptions`](@ref).
It controls how the generated code is realized (returned as an `Expr` vs compiled to a
callable, and whether wrapped in a `GeneratedFunctionWrapper`) and holds a nested
`Symbolics.CodegenFunctionOptions` (`codegen`) with the low-level code-generation options
threaded down to `build_function_wrapper`.

`expression` and `wrap_gfw` are `Bool` type parameters (each accepts `Val{true}`/`Val{false}`
or a `Bool` in the keyword constructor). Both gate the return type of the generators —
`expression` selects `Expr` vs a compiled callable, `wrap_gfw` selects wrapping in a
`GeneratedFunctionWrapper` — so making them type parameters lets the generators branch on
them statically. All other options are fields, so for a fixed `(expression, wrap_gfw)` the
type is invariant to their values.

# Keyword arguments

- `expression`, `wrap_gfw`: as above (lifted into the type parameters).
- `eval_expression`: whether to `eval` the generated code (vs `RuntimeGeneratedFunctions`).
- `eval_module`: the module used to realize the generated code.
- `compiler_options`: a [`CompilerOptions`](@ref).
- `codegen_function_options`: a `Symbolics.CodegenFunctionOptions` with the code-generation
  options forwarded down to `build_function_wrapper`.
"""
struct GeneratedFunctionOptions{expression, wrap_gfw}
    eval_expression::Bool
    eval_module::Module
    compiler_options::CompilerOptions
    codegen::Symbolics.CodegenFunctionOptions
end

function GeneratedFunctionOptions(;
        expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__,
        compiler_options::CompilerOptions = CompilerOptions(),
        codegen_function_options::Symbolics.CodegenFunctionOptions = Symbolics.CodegenFunctionOptions()
    )
    return GeneratedFunctionOptions{_gfo_bool(expression), _gfo_bool(wrap_gfw)}(
        eval_expression, eval_module, compiler_options, codegen_function_options
    )
end

# Recover the `expression`/`wrap_gfw` switches as `Val` types for dispatch into
# `maybe_compile_function` / for branching in generator bodies.
expression_val(::GeneratedFunctionOptions{E}) where {E} = Val{E}
wrap_gfw_val(::GeneratedFunctionOptions{E, W}) where {E, W} = Val{W}

"""
    ObservedFunctionCache(sys, opts::GeneratedFunctionOptions; steady_state = false)

Build an `ObservedFunctionCache` for `sys`, taking `eval_expression`, `eval_module`,
`checkbounds`, and `optimize` from `opts` (its own `expression` type parameter decides
whether an `Expr` or a live `ObservedFunctionCache` is returned) rather than as separate
keywords. This is the form every `*Function` constructor should use, since they already
have a `GeneratedFunctionOptions` (`opts.codegen` off a `SciMLFunctionOptions`) on hand.
"""
function ObservedFunctionCache(
        sys, opts::GeneratedFunctionOptions{E}; steady_state::Bool = false
    ) where {E}
    (; eval_expression, eval_module) = opts
    checkbounds = opts.codegen.checkbounds
    optimize = opts.codegen.optimize
    return if E
        :(
            $ObservedFunctionCache(
                $sys, Dict(), $steady_state, $eval_expression,
                $eval_module, $checkbounds, $optimize
            )
        )
    else
        ObservedFunctionCache(
            sys, Dict(), steady_state, eval_expression, eval_module, checkbounds,
            optimize,
        )
    end
end

"""
    ObservedFunctionCache(sys; expression, steady_state, eval_expression, eval_module, checkbounds, optimize)

Backward-compatible keyword constructor, kept for callers that don't have a
`GeneratedFunctionOptions` on hand. Assembles one from the given keywords and forwards to
[`ObservedFunctionCache(sys, ::GeneratedFunctionOptions)`](@ref).
"""
function ObservedFunctionCache(
        sys; expression = Val{false}, steady_state::Bool = false, eval_expression = false,
        eval_module = @__MODULE__, checkbounds = true, optimize = nothing,
    )
    opts = GeneratedFunctionOptions(;
        expression, eval_expression, eval_module,
        codegen_function_options = Symbolics.CodegenFunctionOptions(; checkbounds, optimize)
    )
    return ObservedFunctionCache(sys, opts; steady_state)
end

"""
    SciMLFunctionOptions{expression}(; kwargs...)

Bundle of options shared by the `SciMLBase.*Function` constructors (`ODEFunction`,
`DAEFunction`, `SDEFunction`, `NonlinearFunction`, `OptimizationFunction`,
`DDEFunction`, `SDDEFunction`, `DiscreteFunction`, `ImplicitDiscreteFunction`,
`IntervalNonlinearFunction`, `SemilinearODEFunction`, `ODEInputFunction`). Holds a
nested [`GeneratedFunctionOptions`](@ref) (`codegen`) for the code-generation
options, plus the options specific to assembling a `*Function` (Jacobian/tgrad
requests, sparsity, initialization data, the `u0`/`p`/`t` specialization inputs).

`expression` is a `Bool` **type parameter** (not a field), matching
`GeneratedFunctionOptions`: `SciMLFunctionOptions{true}` requests code generation
that returns an `Expr`, `SciMLFunctionOptions{false}` requests a runtime-callable
function. Everything else is a field, so for a fixed `expression` the type is
invariant to the option *values* — constructors that thread this struct through
are not recompiled per keyword-set, but they *are* specialized (deliberately) on
`expression`. `wrap_gfw` is not exposed here: every `*Function` constructor
builds its `codegen` with `wrap_gfw = Val{true}`, so it is fixed rather than
threaded through as an option.

Not every field is consumed by every constructor (e.g. `OptimizationFunction`
does not consume `analytic`/`initialization_data`; `ODEInputFunction` does not
consume `expression`/`check_compatibility`) — these are documented no-ops on
constructors that don't have a use for them, rather than a reason to keep the
field out of the shared struct. `t` in particular is only relevant for
time-dependent functions; constructors for time-independent problems (e.g.
`NonlinearFunction`, `OptimizationFunction`, `IntervalNonlinearFunction`) simply
don't consume it.

# Keyword arguments

- `u0`, `p`, `t`: the operating point / independent variable inputs used for
  `iip`/`spec` specialization.
- `jac`, `tgrad`, `sparsity`: whether to generate a Jacobian / time-gradient, and
  whether to report Jacobian sparsity.
- `sparse`: whether generated Jacobians/mass matrices should be sparse.
- `analytic`: an optional analytic solution function.
- `simplify`: whether to run `SymbolicUtils.simplify` on the symbolic
  derivative-related quantities (Jacobians, time-gradients, mass matrices, and,
  for `OptimizationFunction`, cost/constraint gradients and Hessians) after
  differentiation and before code generation. It never affects RHS codegen.
- `initialization_data`: optional `OverrideInitData` for problem initialization.
- `check_compatibility`: whether to validate the system against the target
  `*Function`/`*Problem` type before assembling it.
- `expression`, `eval_expression`, `eval_module`, `compiler_options`, and the
  `Symbolics.CodegenFunctionOptions` fields (`checkbounds`, `optimize`,
  `nanmath`, `wrap_code`, `iip_config`, `sort_addmul`, `similarto`,
  `outputidxs`, `skipzeros`): forwarded into the nested `codegen`.

Any other keyword is accepted and silently dropped (forwarded into, and
discarded by, `Symbolics.CodegenFunctionOptions`'s own trailing `kwargs...`).
This is the backward-compatibility escape hatch for the `kwargs...` wrappers
each `*Function` constructor keeps around this struct: unknown keywords —
whether genuinely irrelevant or leaked in from unrelated call sites (e.g.
MTKBase's initialization-problem plumbing, which forwards a problem's keywords
indiscriminately to every constructor it touches) — are ignored rather than
raising a `MethodError`.
"""
struct SciMLFunctionOptions{expression}
    codegen::GeneratedFunctionOptions{expression, true}
    u0::Any
    p::Any
    t::Any
    jac::Bool
    tgrad::Bool
    sparse::Bool
    sparsity::Bool
    analytic::Any
    simplify::Bool
    initialization_data::Any
    check_compatibility::Bool
end

function SciMLFunctionOptions(;
        u0 = nothing, p = nothing, t = nothing,
        jac::Bool = false, tgrad::Bool = false, sparse::Bool = false,
        sparsity::Bool = false, analytic = nothing, simplify::Bool = false,
        initialization_data = nothing, expression = Val{false},
        check_compatibility::Bool = true,
        eval_expression::Bool = false, eval_module::Module = @__MODULE__,
        compiler_options::CompilerOptions = CompilerOptions(),
        checkbounds::Bool = false, optimize = nothing,
        nanmath::Bool = true, wrap_code::Tuple = (identity, identity),
        iip_config::NTuple{2, Bool} = (true, true), sort_addmul::Bool = false,
        similarto = nothing, outputidxs = nothing, skipzeros::Bool = false,
        kwargs...,
    )
    E = _gfo_bool(expression)
    # The leftover `kwargs...` here is the backward-compatibility escape hatch: callers of
    # the `*Function` constructors' `kwargs...` wrappers may pass keywords that belong to
    # neither `SciMLFunctionOptions` nor `Symbolics.CodegenFunctionOptions` (e.g. options
    # that only apply to a sibling constructor, or that leak in through MTKBase's
    # initialization-problem plumbing). `CodegenFunctionOptions` has its own trailing
    # `kwargs...` that silently discards anything it doesn't recognize, so nothing needs
    # to be named here for this to be safe — unknown keywords are dropped, not errored on.
    codegen = GeneratedFunctionOptions(;
        expression = Val{E}, wrap_gfw = Val{true}, eval_expression, eval_module,
        compiler_options,
        codegen_function_options = Symbolics.CodegenFunctionOptions(;
            nanmath, wrap_code, checkbounds, iip_config, sort_addmul, optimize,
            similarto, outputidxs, skipzeros, kwargs...,
        )
    )
    return SciMLFunctionOptions{E}(
        codegen, u0, p, t, jac, tgrad, sparse, sparsity, analytic, simplify,
        initialization_data, check_compatibility,
    )
end

# Recover the `expression` switch as a `Val` type for dispatch/return in method bodies.
expression_val(::SciMLFunctionOptions{E}) where {E} = Val{E}

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
    BuildFunctionWrapperOptions(; kwargs...)

Options for [`build_function_wrapper`](@ref). Bundles the wrapper-specific options
together with a nested `Symbolics.CodegenFunctionOptions` (the `codegen` field)
holding the code-generation options that are ultimately forwarded to
`codegen_function`.

Threading a single `BuildFunctionWrapperOptions` (one concrete type) instead of a
`kwargs...` bundle means the body of `build_function_wrapper` is compiled once
regardless of which options are set. Keeping the code-generation options in a
nested `CodegenFunctionOptions` means that if `codegen_function` gains new options
this struct does not need to change — they ride along inside `codegen`.

The keyword constructor mirrors the historical keyword arguments of
`build_function_wrapper`. The code-generation options (including `wrap_code` and
`optimize`) are collected into the nested `codegen`. `build_function_wrapper` reads
`wrap_code`/`optimize` from `codegen` and writes the transformed values back — the
composed `wrap_code`, the resolved `optimize`, and the `output_type`-derived
`similarto` — before calling `codegen_function`. Keywords not recognized by
`CodegenFunctionOptions` are ignored, matching the previous behaviour.

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
- `codegen_function_options`: A [`Symbolics.CodegenFunctionOptions`](@ref) holding the
  code-generation options forwarded to `codegen_function`. `build_function_wrapper` reads
  `wrap_code` and `optimize` from it and writes the transformed values back (the composed
  `wrap_code`, the resolved `optimize`, and the `output_type`-derived `similarto`) before
  generating code. The backwards-compatible keyword form of `build_function_wrapper`
  assembles this from loose keyword arguments; direct callers of this constructor must pass
  a fully-formed `CodegenFunctionOptions` (the constructor accepts no other keywords).
"""
struct BuildFunctionWrapperOptions
    p_start::Int
    # `nothing` means "resolve from `sys`/`args` in the struct-based method"
    p_end::Union{Nothing, Int}
    compress_args::Vector{UnitRange{Int}}
    non_standard_param_layout::Bool
    u_arg::Int
    # `nothing` means "resolve to `is_dde(sys)` in the struct-based method"
    wrap_delays::Union{Nothing, Bool}
    histfn::Any
    histfn_symbolic::SymbolicT
    create_bindings::Bool
    output_type::Any
    mkarray::Any
    wrap_mtkparameters::Bool
    extra_assignments::Vector{Assignment}
    n_param_buffers::Union{Nothing, Int}
    codegen::Symbolics.CodegenFunctionOptions
end

function BuildFunctionWrapperOptions(;
        p_start::Integer = 2, p_end = nothing, compress_args = UnitRange{Int}[],
        non_standard_param_layout = false, u_arg::Integer = -1, wrap_delays = nothing,
        histfn = DDE_HISTORY_FUN, histfn_symbolic = histfn,
        create_bindings = false, output_type = nothing, mkarray = nothing,
        wrap_mtkparameters = true, extra_assignments = Assignment[],
        n_param_buffers = nothing,
        codegen_function_options::Symbolics.CodegenFunctionOptions = Symbolics.CodegenFunctionOptions()
    )
    # The constructor only accepts valid `build_function_wrapper` options; there is no
    # `kwargs...` sink. Code-generation options are supplied as a fully-formed
    # `CodegenFunctionOptions` via `codegen_function_options` (`build_function_wrapper` reads
    # `wrap_code`/`optimize` from it and writes the transformed values back before calling
    # `codegen_function`). The backwards-compatible keyword form of `build_function_wrapper`
    # is responsible for assembling that `CodegenFunctionOptions` from loose keywords.
    return BuildFunctionWrapperOptions(
        p_start, p_end, compress_args, non_standard_param_layout, u_arg, wrap_delays,
        histfn, histfn_symbolic, create_bindings, output_type, mkarray,
        wrap_mtkparameters, extra_assignments, n_param_buffers, codegen_function_options
    )
end

"""
    build_function_wrapper(sys::AbstractSystem, expr, args...; kwargs...)

Backwards-compatibility keyword-argument form of `build_function_wrapper`. The keyword
arguments (documented on [`BuildFunctionWrapperOptions`](@ref)) are bundled into a
`BuildFunctionWrapperOptions` and forwarded to the primary method,
[`build_function_wrapper(sys, expr, args, opts::BuildFunctionWrapperOptions)`](@ref). This
method exists only for backwards compatibility; new code should construct a
`BuildFunctionWrapperOptions` and call that method directly.
"""
Base.@nospecializeinfer function build_function_wrapper(
        sys::AbstractSystem, @nospecialize(expr), @nospecialize(args...); p_start = 2,
        p_end = nothing, compress_args = UnitRange{Int}[],
        non_standard_param_layout = false, u_arg::Integer = -1,
        wrap_delays = is_dde(sys), histfn = DDE_HISTORY_FUN, histfn_symbolic = histfn,
        wrap_code = (identity, identity),
        create_bindings = false, @nospecialize(output_type::Union{Nothing, Type} = nothing), mkarray = nothing,
        wrap_mtkparameters = true, extra_assignments = Assignment[],
        n_param_buffers::Int = -1, optimize = nothing, @nospecialize(kwargs...)
    )
    # Thin, backwards-compatible keyword entry point: bundle the options into a single
    # `BuildFunctionWrapperOptions` (computing the `sys`-dependent defaults here, where
    # `sys`/`args` are in scope) and dispatch to the struct-based method below, whose body
    # is compiled once regardless of the options. This shim is the one place that accepts
    # loose code-generation keywords: `wrap_code` and `optimize` (which the struct-based
    # method reads back out of `codegen`) together with any other keywords in `kwargs...`
    # ride into a `CodegenFunctionOptions` (which drops any it does not recognize) passed
    # via `codegen_function_options`.
    opts = BuildFunctionWrapperOptions(;
        p_start, p_end, compress_args, non_standard_param_layout, u_arg, wrap_delays,
        histfn, histfn_symbolic, create_bindings, output_type, mkarray,
        wrap_mtkparameters, extra_assignments, n_param_buffers,
        codegen_function_options = Symbolics.CodegenFunctionOptions(; wrap_code, optimize, kwargs...)
    )
    return build_function_wrapper(sys, expr, collect(Any, args), opts)
end

"""
    build_function_wrapper(sys::AbstractSystem, expr, args, opts::BuildFunctionWrapperOptions)

A wrapper around `build_function` which performs the necessary transformations for
code generation of all types of systems. `expr` is the expression returned from the
generated functions, and `args` is the `Vector{Any}` of arguments.

Options are supplied as a [`BuildFunctionWrapperOptions`](@ref); see its docstring for the
available options. This is the primary method — the keyword-argument form of
`build_function_wrapper` is a backwards-compatibility shim that bundles its keywords into a
`BuildFunctionWrapperOptions` and calls this method.
"""
Base.@nospecializeinfer function build_function_wrapper(
        sys::AbstractSystem, @nospecialize(expr), args::Vector{Any},
        opts::BuildFunctionWrapperOptions
    )
    (;
        p_start, p_end, compress_args, non_standard_param_layout, u_arg, wrap_delays,
        histfn_symbolic, create_bindings, output_type, mkarray, wrap_mtkparameters,
        extra_assignments, n_param_buffers,
    ) = opts
    # Resolve the `sys`/`args`-dependent defaults (left as `nothing` when constructing the
    # options without a `sys` in scope). `length(args)` here is the original argument count,
    # before `wrap_delays` may insert the history-function argument below.
    p_end = p_end === nothing ? (is_time_dependent(sys) ? length(args) - 1 : length(args)) : p_end
    wrap_delays = wrap_delays === nothing ? is_dde(sys) : wrap_delays
    # `-1` is the "unset" sentinel used below (`pbuf_end`). Direct callers may leave
    # `n_param_buffers` as `nothing`; normalize it to the sentinel here.
    n_param_buffers = n_param_buffers === nothing ? -1 : n_param_buffers
    isscalar = !(expr isa AbstractArray || symbolic_type(expr) == ArraySymbolic())
    # Copy so the in-place mutations below (`insert!`, `deleteat!`, `setindex!`) do not
    # affect the caller's vector.
    args = copy(args)
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

    wrap_code = opts.codegen.wrap_code .∘ wrap_assignments(false, assignments)

    # handling of `output_type` and `mkarray`. A `similarto` already present on
    # `opts.codegen` (set directly via `codegen_function_options` by callers such as
    # `generate_history`, bypassing the `output_type` wrapper option) takes precedence over
    # the `output_type`-derived value, matching the pre-`BuildFunctionWrapperOptions`
    # behaviour where such an explicit `similarto` was passed as a later keyword argument.
    similarto = opts.codegen.similarto
    if output_type === Tuple
        expr = MakeTuple(Tuple(expr))
    elseif mkarray === nothing
        similarto = similarto === nothing ? output_type : similarto
    else
        expr = mkarray(expr, output_type)
    end

    optimize = resolve_optimize_option(opts.codegen.optimize)
    # Write the wrapper-derived code-generation options (the composed `wrap_code`, the
    # `output_type`-derived `similarto`, and the resolved `optimize`) into the nested
    # `CodegenFunctionOptions`. `setproperties` overrides only these fields and preserves the
    # rest (`nanmath`, `checkbounds`, `iip_config`, ...), so `codegen_function` still sees a
    # single concrete `CodegenFunctionOptions` type and is compiled once.
    codegen_options = setproperties(opts.codegen, (; wrap_code, similarto, optimize))
    return Symbolics.codegen_function(ir, expr, args, codegen_options)
end

resolve_optimize_option(x) = x
resolve_optimize_option(::Nothing) = nothing

"""
    $(TYPEDEF)

A wrapper around a generated in-place and out-of-place function. The type-parameter `P`
must be a `Tuple` type `Tuple{A, B, C}` where the first element `A` is the index of the
parameter object in the arguments, the second `B` is the expected number of arguments in
the out-of-place variant of the function, and the third `C` is a boolean indicating whether
the generated functions are for a split system. Encoding these as a `Tuple` type (rather
than a tuple value) allows dispatching on the individual elements. For scalar functions,
the inplace variant can be `nothing`.

For backward compatibility, `P` may also be provided as the old 3-tuple *value* form
`(A, B, C)`, which is converted to the `Tuple{A, B, C}` type form.
"""
struct GeneratedFunctionWrapper{P, O, I} <: Function
    f_oop::O
    f_iip::I
end

# Normalize the `P` type parameter to a `Tuple{...}` type. The old form passed a tuple
# *value* (e.g. `(2, 3, true)`); the new form is a `Tuple` *type* (e.g. `Tuple{2, 3, true}`).
_gfw_params_type(P::Tuple) = Tuple{P...}
_gfw_params_type(::Type{P}) where {P <: Tuple} = P

function GeneratedFunctionWrapper{P}(foop::O, fiip::I) where {P, O, I}
    return GeneratedFunctionWrapper{_gfw_params_type(P), O, I}(foop, fiip)
end

# The wrapped functions are stateless generated functions, so there is nothing to
# deep-copy. Overriding this as the identity avoids `deepcopy` recursing into the
# function internals, which improves `juliac` trimmability.
Base.deepcopy_internal(gfw::GeneratedFunctionWrapper, ::IdDict) = gfw

function GeneratedFunctionWrapper{P}(::Type{Val{true}}, foop, fiip; kwargs...) where {P}
    return :($(GeneratedFunctionWrapper{_gfw_params_type(P)})($foop, $fiip))
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

function (gfw::GeneratedFunctionWrapper{Tuple{PIdx, NArgs, Split}})(args::Vararg{Any, NArgs}) where {PIdx, NArgs, Split}
    # non-split systems just call it as-is
    Split || return gfw.f_oop(args...)
    if args[PIdx] isa Union{Tuple, MTKParameters} && !(args[PIdx] isa Tuple{Vararg{Number}})
        # for split systems, call it as-is if the parameter object is a tuple or MTKParameters
        # but not if it is a tuple of numbers
        return gfw.f_oop(args...)
    end
    # The user provided a single buffer/tuple for the parameter object, so wrap that
    # one in a tuple
    @set! args[PIdx] = (args[PIdx], nothing)
    return gfw.f_oop(args...)
end

function (gfw::GeneratedFunctionWrapper{Tuple{PIdx, NArgs, Split}})(args::Vararg{Any, N}) where {PIdx, NArgs, Split, N}
    # IIP case has one more argument
    if NArgs + 1 != N
        throw(MethodError(gfw, args))
    end
    Split || return gfw.f_iip(args...)
    if args[PIdx + 1] isa Union{Tuple, MTKParameters} && !(args[PIdx + 1] isa Tuple{Vararg{Number}})
        return gfw.f_iip(args...)
    end
    @set! args[PIdx + 1] = (args[PIdx + 1], nothing)
    return gfw.f_iip(args...)
end

function SciMLBase.numargs(::GeneratedFunctionWrapper{Tuple{P, N, S}}) where {P, N, S}
    return (N, N + 1)
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
