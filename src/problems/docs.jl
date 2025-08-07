const U0_P_DOCS = """
The order of unknowns is determined by `unknowns(sys)`. If the system is split
[`is_split`](@ref) create an [`MTKParameters`](@ref) object. Otherwise, a parameter vector.
Initial values provided in terms of other variables will be symbolically evaluated.
The type of `op` will be used to determine the type of the containers. For example, if
given as an `SArray` of key-value pairs, `u0` will be an appropriately sized `SVector`
and the parameter object will be an `MTKParameters` object with `SArray`s inside.
"""

const EVAL_EXPR_MOD_KWARGS = """
- `eval_expression`: Whether to compile any functions via `eval` or
  `RuntimeGeneratedFunctions`.
- `eval_module`: If `eval_expression == true`, the module to `eval` into. Otherwise, the
  module in which to generate the `RuntimeGeneratedFunction`.
"""

const INITIALIZEPROB_KWARGS = """
- `guesses`: The guesses for variables in the system, used as initial values for the
  initialization problem.
- `warn_initialize_determined`: Warn if the initialization system is under/over-determined.
- `initialization_eqs`: Extra equations to use in the initialization problem.
- `fully_determined`: Override whether the initialization system is fully determined.
- `use_scc`: Whether to use `SCCNonlinearProblem` for initialization if the system is fully
  determined.
"""

const PROBLEM_KWARGS = """
$EVAL_EXPR_MOD_KWARGS
$INITIALIZEPROB_KWARGS
- `check_initialization_units`: Enable or disable unit checks when constructing the
  initialization problem.
- `tofloat`: Passed to [`varmap_to_vars`](@ref) when building the parameter vector of
  a non-split system.
- `u0_eltype`: The `eltype` of the `u0` vector. If `nothing`, finds the promoted floating point
  type from `op`.
- `u0_constructor`: A function to apply to the `u0` value returned from
  [`varmap_to_vars`](@ref).
  to construct the final `u0` value.
- `p_constructor`: A function to apply to each array buffer created when constructing the
  parameter object.
- `warn_cyclic_dependency`: Whether to emit a warning listing out cycles in initial
  conditions provided for unknowns and parameters.
- `circular_dependency_max_cycle_length`: Maximum length of cycle to check for. Only
  applicable if `warn_cyclic_dependency == true`.
- `circular_dependency_max_cycles`: Maximum number of cycles to check for. Only applicable
  if `warn_cyclic_dependency == true`.
- `substitution_limit`: The number times to substitute initial conditions into each other
  to attempt to arrive at a numeric value.
"""

const TIME_DEPENDENT_PROBLEM_KWARGS = """
- `callback`: An extra callback or `CallbackSet` to add to the problem, in addition to the
  ones defined symbolically in the system.
"""

const PROBLEM_INTERNALS_HEADER = """
# Extended docs

The following API is internal and may change or be removed without notice. Its usage is
highly discouraged.
"""

const INTERNAL_INITIALIZEPROB_KWARGS = """
- `time_dependent_init`: Whether to build a time-dependent initialization for the problem. A
  time-dependent initialization solves for a consistent `u0`, whereas a time-independent one
  only runs parameter initialization.
- `algebraic_only`: Whether to build the initialization problem using only algebraic equations.
- `allow_incomplete`: Whether to allow incomplete initialization problems.
"""

const PROBLEM_INTERNAL_KWARGS = """
- `build_initializeprob`: If `false`, avoids building the initialization problem.
- `check_length`: Whether to check the number of equations along with number of unknowns and
  length of `u0` vector for consistency. If `false`, do not check with equations. This is
  forwarded to `check_eqs_u0`.
$INTERNAL_INITIALIZEPROB_KWARGS
"""

function problem_ctors(prob, istd)
    if istd
        """
            SciMLBase.$prob(sys::System, op, tspan::NTuple{2}; kwargs...)
            SciMLBase.$prob{iip}(sys::System, op, tspan::NTuple{2}; kwargs...)
            SciMLBase.$prob{iip, specialize}(sys::System, op, tspan::NTuple{2}; kwargs...)
        """
    else
        """
            SciMLBase.$prob(sys::System, op; kwargs...)
            SciMLBase.$prob{iip}(sys::System, op; kwargs...)
            SciMLBase.$prob{iip, specialize}(sys::System, op; kwargs...)
        """
    end
end

function prob_fun_common_kwargs(T, istd)
    return """
    - `check_compatibility`: Whether to check if the given system `sys` contains all the
      information necessary to create a `$T` and no more. If disabled, assumes that `sys`
      at least contains the necessary information.
    - `expression`: `Val{true}` to return an `Expr` that constructs the corresponding
      problem instead of the problem itself. `Val{false}` otherwise.
      $(istd ? " Constructing the expression does not support callbacks" : "")
    """
end

function problem_docstring(prob, func, istd; init = true, extra_body = "")
    if func isa DataType
        func = "`$func`"
    end
    return """
    $(problem_ctors(prob, istd))

    Build a `$prob` given a system `sys` and operating point `op`
    $(istd ? " and timespan `tspan`" : ""). `iip` is a boolean indicating whether the
    problem should be in-place. `specialization` is a `SciMLBase.AbstractSpecalize` subtype
    indicating the level of specialization of the $func. The operating point should be an
    iterable collection of key-value pairs mapping variables/parameters in the system to the
    (initial) values they should take in `$prob`. Any values not provided will fallback to
    the corresponding default (if present).

    $(init ? istd ? TIME_DEPENDENT_INIT : TIME_INDEPENDENT_INIT : "")

    $extra_body

    # Keyword arguments

    $PROBLEM_KWARGS
    $(istd ? TIME_DEPENDENT_PROBLEM_KWARGS : "")
    $(prob_fun_common_kwargs(prob, istd))

    All other keyword arguments are forwarded to the $func constructor.

    $PROBLEM_INTERNALS_HEADER

    $PROBLEM_INTERNAL_KWARGS
    """
end

const TIME_DEPENDENT_INIT = """
ModelingToolkit will build an initialization problem where all initial values for
unknowns or observables of `sys` (either explicitly provided or in defaults) will
be constraints. To remove an initial condition in the defaults (without providing
a replacement) give the corresponding variable a value of `nothing` in the operating
point. The initialization problem will also run parameter initialization. See the
[Initialization](@ref initialization) documentation for more information.
"""

const TIME_INDEPENDENT_INIT = """
ModelingToolkit will build an initialization problem that will run parameter
initialization. Since it does not solve for initial values of unknowns, observed
equations will not be initialization constraints. If an initialization equation
of the system must involve the initial value of an unknown `x`, it must be used as
`Initial(x)` in the equation. For example, an equation to be used to solve for parameter
`p` in terms of unknowns `x` and `y` must be provided as `Initial(x) + Initial(y) ~ p`
instead of `x + y ~ p`. See the [Initialization](@ref initialization) documentation
for more information.
"""

const BV_EXTRA_BODY = """
Boundary value conditions are supplied to Systems in the form of a list of constraints.
These equations  should specify values that state variables should take at specific points,
as in `x(0.5) ~ 1`). More general constraints that  should hold over the entire solution,
such as `x(t)^2 + y(t)^2`, should be  specified as one of the equations used to build the
`System`.

If a `System` without `constraints` is specified, it will be treated as an initial value problem. 

```julia
    @parameters g t_c = 0.5
    @variables x(..) y(t) λ(t)
    eqs = [D(D(x(t))) ~ λ * x(t)
           D(D(y)) ~ λ * y - g
           x(t)^2 + y^2 ~ 1]
    cstr = [x(0.5) ~ 1]
    @mtkcompile pend = System(eqs, t; constraints = cstrs)

    tspan = (0.0, 1.5)
    u0map = [x(t) => 0.6, y => 0.8]
    parammap = [g => 1]
    guesses = [λ => 1]

    bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; guesses, check_length = false)
```

If the `System` has algebraic equations, like `x(t)^2 + y(t)^2`, the resulting 
`BVProblem` must be solved using BVDAE solvers, such as Ascher.
"""

for (mod, prob, func, istd, kws) in [
    (SciMLBase, :ODEProblem, ODEFunction, true, (;)),
    (SciMLBase, :SteadyStateProblem, ODEFunction, false, (;)),
    (SciMLBase, :BVProblem, ODEFunction, true,
        (; init = false, extra_body = BV_EXTRA_BODY)),
    (SciMLBase, :DAEProblem, DAEFunction, true, (;)),
    (SciMLBase, :DDEProblem, DDEFunction, true, (;)),
    (SciMLBase, :SDEProblem, SDEFunction, true, (;)),
    (SciMLBase, :SDDEProblem, SDDEFunction, true, (;)),
    (JumpProcesses, :JumpProblem, "inner SciMLFunction", true, (; init = false)),
    (SciMLBase, :DiscreteProblem, DiscreteFunction, true, (;)),
    (SciMLBase, :ImplicitDiscreteProblem, ImplicitDiscreteFunction, true, (;)),
    (SciMLBase, :NonlinearProblem, NonlinearFunction, false, (;)),
    (SciMLBase, :NonlinearLeastSquaresProblem, NonlinearFunction, false, (;)),
    (SciMLBase, :SCCNonlinearProblem, NonlinearFunction, false, (; init = false)),
    (SciMLBase, :OptimizationProblem, OptimizationFunction, false, (; init = false))
]
    @eval @doc problem_docstring($mod.$prob, $func, $istd) $mod.$prob
end

function function_docstring(func, istd, optionals)
    return """
        $func(sys::System; kwargs...)
        $func{iip}(sys::System; kwargs...)
        $func{iip, specialize}(sys::System; kwargs...)

    Create a `$func` from the given `sys`. `iip` is a boolean indicating whether the
    function should be in-place. `specialization` is a `SciMLBase.AbstractSpecalize`
    subtype indicating the level of specialization of the $func.

    Beyond the arguments listed below, this constructor accepts all keyword arguments
    supported by the DifferentialEquations.jl `solve` function. For a complete list
    and detailed descriptions, see the [DifferentialEquations.jl solve documentation](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/).

    # Keyword arguments

    - `u0`: The `u0` vector for the corresponding problem, if available. Can be obtained
      using [`ModelingToolkit.get_u0`](@ref).
    - `p`: The parameter object for the corresponding problem, if available. Can be obtained
      using [`ModelingToolkit.get_p`](@ref).
    $(istd ? TIME_DEPENDENT_FUNCTION_KWARGS : "")
    $EVAL_EXPR_MOD_KWARGS
    - `checkbounds`: Whether to enable bounds checking in the generated code.
    - `simplify`: Whether to `simplify` any symbolically computed jacobians/hessians/etc.
    - `cse`: Whether to enable Common Subexpression Elimination (CSE) on the generated code.
      This typically improves performance of the generated code but reduces readability.
    - `sparse`: Whether to generate jacobian/hessian/etc. functions that return/operate on
      sparse matrices. Also controls whether the mass matrix is sparse, wherever applicable.
    $(prob_fun_common_kwargs(func, istd))
    $(process_optional_function_kwargs(optionals))
    - `kwargs...`: Additional keyword arguments passed to the solver

    All other keyword arguments are forwarded to the `$func` struct constructor.
    """
end

const TIME_DEPENDENT_FUNCTION_KWARGS = """
- `t`: The initial time for the corresponding problem, if available.
"""

const JAC_KWARGS = """
- `jac`: Whether to symbolically compute and generate code for the jacobian function.
"""

const TGRAD_KWARGS = """
- `tgrad`: Whether to symbolically compute and generate code for the `tgrad` function.
"""

const SPARSITY_KWARGS = """
- `sparsity`: Whether to provide symbolically compute and provide sparsity patterns for the
  jacobian/hessian/etc.
"""

const RESID_PROTOTYPE_KWARGS = """
- `resid_prototype`: The prototype of the residual function `f` for a problem involving a
  nonlinear solve where the residual and `u0` have different sizes.
"""

const GRAD_KWARGS = """
- `grad`: Whether the symbolically compute and generate code for the gradient of the cost
  function with respect to unknowns.
"""

const HESS_KWARGS = """
- `hess`: Whether to symbolically compute and generate code for the hessian function.
"""

const CONSH_KWARGS = """
- `cons_h`: Whether to symbolically compute and generate code for the hessian function of
  constraints. Since the constraint function is vector-valued, the hessian is a vector
  of hessian matrices.
"""

const CONSJ_KWARGS = """
- `cons_j`: Whether to symbolically compute and generate code for the jacobian function of
  constraints.
"""

const CONSSPARSE_KWARGS = """
- `cons_sparse`: Identical to the `sparse` keyword, but specifically for jacobian/hessian
  functions of the constraints.
"""

const INPUTFN_KWARGS = """
- `inputs`: The variables in the input vector. The system must have been simplified using
  `mtkcompile` with these variables passed as `inputs`.
- `disturbance_inputs`: The disturbance input variables. The system must have been
  simplified using `mtkcompile` with these variables passed as `disturbance_inputs`.
"""

const CONTROLJAC_KWARGS = """
- `controljac`: Whether to symbolically compute and generate code for the jacobian of
  the ODE with respect to the inputs.
"""

const OPTIONAL_FN_KWARGS_DICT = Dict(
    :jac => JAC_KWARGS,
    :tgrad => TGRAD_KWARGS,
    :sparsity => SPARSITY_KWARGS,
    :resid_prototype => RESID_PROTOTYPE_KWARGS,
    :grad => GRAD_KWARGS,
    :hess => HESS_KWARGS,
    :cons_h => CONSH_KWARGS,
    :cons_j => CONSJ_KWARGS,
    :cons_sparse => CONSSPARSE_KWARGS,
    :inputfn => INPUTFN_KWARGS,
    :controljac => CONTROLJAC_KWARGS
)

const SPARSITY_OPTIONALS = Set([:jac, :hess, :cons_h, :cons_j, :controljac])

const CONS_SPARSITY_OPTIONALS = Set([:cons_h, :cons_j])

function process_optional_function_kwargs(choices::Vector{Symbol})
    if !isdisjoint(choices, SPARSITY_OPTIONALS)
        push!(choices, :sparsity)
    end
    if !isdisjoint(choices, CONS_SPARSITY_OPTIONALS)
        push!(choices, :cons_sparse)
    end
    join(map(Base.Fix1(getindex, OPTIONAL_FN_KWARGS_DICT), choices), "\n")
end

for (mod, func, istd, optionals) in [
    (SciMLBase, :ODEFunction, true, [:jac, :tgrad]),
    (SciMLBase, :ODEInputFunction, true, [:inputfn, :jac, :tgrad, :controljac]),
    (SciMLBase, :DAEFunction, true, [:jac, :tgrad]),
    (SciMLBase, :DDEFunction, true, Symbol[]),
    (SciMLBase, :SDEFunction, true, [:jac, :tgrad]),
    (SciMLBase, :SDDEFunction, true, Symbol[]),
    (SciMLBase, :DiscreteFunction, true, Symbol[]),
    (SciMLBase, :ImplicitDiscreteFunction, true, Symbol[]),
    (SciMLBase, :NonlinearFunction, false, [:resid_prototype, :jac]),
    (SciMLBase, :IntervalNonlinearFunction, false, Symbol[]),
    (SciMLBase, :OptimizationFunction, false, [:jac, :grad, :hess, :cons_h, :cons_j])
]
    @eval @doc function_docstring($mod.$func, $istd, $optionals) $mod.$func
end

@doc """
    SciMLBase.HomotopyNonlinearFunction(sys::System; kwargs...)
    SciMLBase.HomotopyNonlinearFunction{iip}(sys::System; kwargs...)
    SciMLBase.HomotopyNonlinearFunction{iip, specialize}(sys::System; kwargs...)

Create a `HomotopyNonlinearFunction` from the given `sys`. `iip` is a boolean indicating
whether the function should be in-place. `specialization` is a `SciMLBase.AbstractSpecalize`
subtype indicating the level of specialization of the $func.

# Keyword arguments

- `u0`: The `u0` vector for the corresponding problem, if available. Can be obtained
  using [`ModelingToolkit.get_u0`](@ref).
- `p`: The parameter object for the corresponding problem, if available. Can be obtained
  using [`ModelingToolkit.get_p`](@ref).
$EVAL_EXPR_MOD_KWARGS
- `checkbounds`: Whether to enable bounds checking in the generated code.
- `simplify`: Whether to `simplify` any symbolically computed jacobians/hessians/etc.
- `cse`: Whether to enable Common Subexpression Elimination (CSE) on the generated code.
  This typically improves performance of the generated code but reduces readability.
- `fraction_cancel_fn`: The function to use to simplify fractions in the polynomial
  expression. A more powerful function can increase processing time but be able to
  eliminate more rational functions, thus improving solve time. Should be a function that
  takes a symbolic expression containing zero or more fraction expressions and returns the
  simplified expression. While this defaults to `SymbolicUtils.simplify_fractions`, a viable
  alternative is `SymbolicUtils.quick_cancel`

All keyword arguments are forwarded to the wrapped `NonlinearFunction` constructor.
""" SciMLBase.HomotopyNonlinearFunction

@doc """
    SciMLBase.IntervalNonlinearProblem(sys::System, uspan::NTuple{2}, parammap = SciMLBase.NullParameters(); kwargs...)

Create an `IntervalNonlinearProblem` from the given `sys`. This is only valid for a system
of nonlinear equations with a single equation and unknown. `uspan` is the interval in which
the root is to be found, and `parammap` is an iterable collection of key-value pairs
providing values for the parameters in the system.

$TIME_INDEPENDENT_INIT

# Keyword arguments

$PROBLEM_KWARGS
$(prob_fun_common_kwargs(IntervalNonlinearProblem, false))

All other keyword arguments are forwarded to the `IntervalNonlinearFunction` constructor.

$PROBLEM_INTERNALS_HEADER

$PROBLEM_INTERNAL_KWARGS
""" SciMLBase.IntervalNonlinearProblem

@doc """
    SciMLBase.LinearProblem(sys::System, op; kwargs...)
    SciMLBase.LinearProblem{iip}(sys::System, op; kwargs...)

Build a `LinearProblem` given a system `sys` and operating point `op`. `iip` is a boolean
indicating whether the problem should be in-place. The operating point should be an
iterable collection of key-value pairs mapping variables/parameters in the system to the
(initial) values they should take in `LinearProblem`. Any values not provided will
fallback to the corresponding default (if present).

Note that since `u0` is optional for `LinearProblem`, values of unknowns do not need to be
specified in `op` to create a `LinearProblem`. In such a case, `prob.u0` will be `nothing`
and attempting to symbolically index the problem with an unknown, observable, or expression
depending on unknowns/observables will error.

Updating the parameters automatically updates the `A` and `b` arrays.

# Keyword arguments

$PROBLEM_KWARGS
$(prob_fun_common_kwargs(LinearProblem, false))

All other keyword arguments are forwarded to the $func constructor.

$PROBLEM_INTERNALS_HEADER

$PROBLEM_INTERNAL_KWARGS
""" SciMLBase.LinearProblem
