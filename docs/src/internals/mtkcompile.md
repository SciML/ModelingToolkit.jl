# The `mtkcompile` pipeline

!!! warning "Internal API"
    The content on this page describes internal implementation details of ModelingToolkit.
    These APIs may change without notice.

`mtkcompile` is one of the most complicated pieces in ModelingToolkit.jl. This document
aims to go over the important parts of how it works.

## The basic version

The `mtkcompile` function itself is defined in ModelingToolkitBase.jl. It is pretty barebones,
and the only things it really implements are:
- `additional_passes`.
- Calling `complete`.

`mtkcompile` hands off to `_mtkcompile`, which filters out systems needing special handling:
ones with poissonians, jumps, brownians, or noise, etc. This then calls the special
`__mtkcompile`. ModelingToolkitBase defines a very barebones version of this function, which
aims to implement the essentials. It handles `inputs`, `outputs` and `disturbance_inputs`. It
supports explicit differential equations of the form `D(x) ~ ...` and performs a very basic
version of tearing to eliminate observed equations.

## The complex version

ModelingToolkit.jl then pirates `__mtkcompile` to define the advanced simplification. It is
the only package allowed to do so. First, this identifies the state machine subsystems written
using the experimental state machine syntax. These are not compiled in this package, simply
identified and tracked. Next, the `expand_connections` pass expands `connect` equations using
the Modelica connection semantics. In this process, we also retain source information about
which subcomponent each non-`connect` equation came from.

A pass called `discover_maybe_zeros` uses `bindings` and `initial_conditions` to identify
variables and parameters that might be zero in the system. This helps the rest of `mtkcompile`
make better simplification decisions and avoid singularities. `TearingState` is a data structure
that contains symbolic and structural information about the system. It represents the system
as a list of equations, a list of variables, and some graphs to describe them. SDEs are identified
here. They are simplified by removing brownian variables, simplifying the resultant DAE and
introducing them again.

DAE simplification is then handled by `mtkcompile!`. This mutates the `state::TearingState` which
contains a reference to the system. At the end of this process, `state.sys` will be the simplified
system. The `TearingState` constructor is passed `defer_scalarization = true`, which means it does not
scalarize the equations of the system eagerly. This is useful because `mtkcompile!` performs clock
inference and handles hybrid systems, where deferred scalarization is useful to preserve information.

For non-time-dependent information, the concept of a hybrid system does not apply. They are immediately
scalarized and simplified like normal. Nonlinear systems of equations hit this path. They are
equivalent to DAEs with no differential variables, so the simplification path is identical.
Clock inference happens via the `ClockInference` data structure. This and most of the rest of
the passes are now implemented in StateSelection.jl and ModelingToolkitTearing.jl. Refer to the
clock inference section for details on how it works.

Discrete partitions are handled via the `additional_passes` machanism. The first pass implementing
the `discrete_compile_pass` trait is used as the discrete compiler. It is called here and provided
the results of clock inference along with the simplified continuous subsystem. This pass should
then compile the discrete partitions and return a modified version of the system it is passed in.

`mtkcompile!` calls `_mtkcompile!`. This function starts the structural passes. It first validates
the input and output variables, and converts input variables into parameters. This also applies to
variables that are input to the continuous partition from the discrete partitions. We then have
a sequence of passes that aim to eliminate as much from the system as possible before the more
complicated index reduction and tearing passes.

### `eliminate_perfect_aliases!`

This pass identifies equatins of the form `x ~ y` and `x ~ -y` where `x` and `y` are variables.
`x` and `y` are thus identified as aliases, and chains of aliases are grouped into an alias group.
One variable in each group is chosen as the "target". Each variable in the group is then rewritten
as that target (or its negation). In case a group contains both `x ~ y` and `x ~ -y`, the entire
group must be `0` and is treated as such. Irreducible variables are not eliminated.

`connect` equations primarily form these sorts of aliasing equations. This pass thus greatly
reduces the number of equations and variables that the subsequent simplification has to handle.

### `trivial_tearing!`

This pass identifies variables already written by the user as observed equations. The equations
should have been originally written as explicit equations of the form `x ~ ...`, where `x` is
not differentiated and not present in any other equation. This can then easily be removed. This
sort of elimination can cascade, allowing removal of additional variables one after the other.

### `alias_elimination!`

Alias elimination builds an integer-linear-coefficient matrix. Every equation that can be
written as `Σ(a_i * x_i) ~ 0` where `a_i` are integer coefficients and `x_i` are variables
is a row in this matrix. A handwritten implementation of the Bareiss algorithm is then
applied to this matrix. In doing so, the equations are written down in a simpler form more
useful for the later tearing pass. The Bareiss implementation uses a custom pivot selection
algorithm. Variables present only in such integer-linear-coefficient equations are preferred
as pivots, since they must be solved by one of these equations. The next preference is
for highest order derivative variables, and then any remaning variables.

This pass often rewrites more equations into perfect alias form. Thus, we rerun
`eliminate_perfect_aliases!` after it.

### `eliminate_zero_variables_fixpoint!`

`alias_elimination!` also very frequently eliminates equations into `x ~ 0` form. Retaining
such equations is typically unnecessary. The information that a variable is identically
zero is very useful for subsequent simplification. This pass eliminates such variables and
substitutes the zero into all other equations. Irreducible variables are not eliminated, but
are substituted into other equations. When a variable `x` is eliminated, all derivatives of
`x` present in the system are also similarly eliminated. All antiderivatives of `x` are
analytically integrated. For example, if `D(D(x)) ~ 0` is found, we will also eliminate
`D(x) ~ var"x_t#0"` and `x ~ var"x_t#0" * t + var"x#0"` where `var"x_t#0"` and `var"x#0"`
are additional parameters introduced here. They are automatically solved for by initialization
and do not need special treatment by the user as long as initial conditions for `x` and `D(x)`
can be obtained. Often eliminating such zero variables can identify more zero variables,
so the pass is run in a fixpoint loop for a maximum number of iterations.

### Consistency checks

The system is checked for consistency. The `fully_determined` flag defaults to `true`
and requires that the system is not structurally singular. If it is set to `false`,
under- and over-determined systems are allowed. If it is set to `nothing`, the
pipeline will identify whether the system is structurally singular or not and simplify
accordingly.

Fully determined systems go through index reduction followed by tearing. Non-fully-determined
systems only go through tearing, since the index reduction algorithms will not terminate
on such inputs.

### Index reduction

For a consistent (solvable) system of equations

```math
F(\dot{x}, x, t) = 0
```

The index of the system is defined as "the minimum number of times that all or part of [the system]
must be differentiated with respect to `t` in order to determine $\dot{x}$ as a continuous function
of `x` and `t`" [1]. Most numerical solvers can only solve DAEs of index `<= 1`. ModelingToolkit.jl
uses a version of the dummy derivatives algorithm introduced in [1] for index reduction. This reduces
the system to index 1 (or index 0 if it is already index 0). It is recommended to read the paper [1]
as an introduction to this method. ModelingToolkit.jl calls `StateSelection.dummy_derivative_graph!`
which performs this algorithm on the structural part of `state` (stored in `state.structure`). This
updates the graph data structures, but does not update the symbolic information. This is done later
in what we term a "reassemble" stage. This section does not attempt to describe in detail
the algorithm. The aforementioned paper [1] and [2] are recommended reading, along with the
implementation of the algorithm itself.

ModelingToolkit.jl allows users to specify a priority for each variables. All scalarized elements
of an array variable get the same priority. This priority is used to influence choices made by
the state selection algorithm and is critical for the user to help guide the algorithm to simplify
the system more optimally. Ties in the state priority are resolved using the "canonical rank" of
the variable, which is a numerical value assigned to the variable based on its name. State
priorities are analyzed along the derivative chain for a variable. If any higher order derivative
of `x` has a negative state priority, the minimum such priority is assigned to `x`. Otherwise, the
maximum priority among its higher order derivatives is used. The state selection algorithm also
uses symbolic jacobian information if possible, specifically when all jacobian coefficients are
small integers.

The output of index reduction via the dummy derivatives algorithm consists of:
- A new list of equations, where some new equations are differentiated version of other equations.
- A new list of variables, where some are derivatives or other variables.
- A list of "dummy" derivatives, which are differentiated variables that aren't really derivatives.
  For example, The equation `D(x) ~ x_t` might be added where `x_t` is just an alternative name for
  `D(x)` because `D(x)` can't be solved analytically, and must be obtained from an algebraic equation.
  The algebraic equation then solves for `x_t`, and that indirectly defines the dynamics for `x`. Here,
  `x_t` is the dummy derivative.

### Tearing

The last part of `StateSelection.dummy_derivative_graph!` is running a tearing algorithm (the steps
are merged together a bit). The job of tearing is to causalize the system. Theoretically, the output
of dummy derivatives can be used directly if materialized as a `System`. However, it contains a lot
of redundant equations which can explicitly be solved for certain variables. The actual set of variables
that need to be integrated are usually a very small subset of the total number of variables in the model.
Tearing figures out this small subset. Note that tearing is NP-complete, and as such the tearing
algorithms used by ModelingToolkit are heuristics. Special care is taken in the primary tearing algorithm
(`CarpanzanoTearing`) to make sure the output is reproducible, since it depends on the order of equations
and variables in the system. This is part of why ModelingToolkit.jl does not guarantee the subset
of variables that a given system will reduce to.

Consider the bipartite incidence graph (`state.structure.graph`) of the system as an incidence
matrix. Tearing essentially takes this binary matrix $A$ and reorders it into the structure:

```math
A = \begin{bmatrix}
L & P \\
Q & R
\end{bmatrix}
```

Where $L$ is lower triangular. The goal typically is for $L$ to be as large as possible, so that the
integrator has to handle fewer variables.

This process starts by building a maximal matching of the incidence graph. This matching will eventually
define which equations are solved for which variables. This matching is coupled with the incidence graph
to build a `BipartiteGraphs.DiCMOBiGraph`. The documentation there is a good description of this data
structure. The strongly connected components (SCCs) of this graph identify variables (and, via the matching,
the corresponding equations) that must be solved simultaneously. Tearing then permutes/tears each SCC
individually. All variables in an SCC are unmatched, and gradually tearing finds a way to build the
assignment in each SCC in a causal order.

The output of tearing is:
- A matching mapping each variable to the equation used to solve it. This includes algebraic variables
  matched to the corresponding algebraic variables.
- The same matching, but algebraic variables and equations are unmatched. Selected state variables are
  marked with the `SelectedState` singleton.
- The variable SCCs.

After tearing, the results of all these structural (graph) operations are materialized symbolically
into a `System` in a step known as reassembly.

### Reassembly

Reassembly has many moving parts. The default algorithm used here is
`ModelingToolkitTearing.DefaultReassembleAlgorithm`. I'll cover the different parts in the order
they are run.

Firstly, we run `get_extra_eqs_vars` to get the extra equations/variables for over-/under-determined
systems. These need to be handled separately in some cases.

Discrete systems run `add_additional_history!`. The docstring here explains what it does with an
example. It is necessary to generate convenient sets of equations for discrete systems while
keeping the states minimal.

`substitute_derivative_algevars!` is necessary to handle dummy derivatives. A derivative variable
`D(x)` is identified as a dummy derivative if it is not matched to `SelectedState`. The variable
is renamed to `xˍt` (with more `t` appended for higher order dummy derivatives) and this is
substituted into all equations where it is present. Higher order derivatives of this variable
are also substituted to `D(xˍt)` (and so on).

`generate_derivative_variables!` is used to handle derivative values used in equations where
it is not the selected state. For example, the system

```julia
eqs = [
    D(x) ~ t
    D(x) ^ 2 + y^3 ~ 3
]
```

Has `D(x)` selected in the first equation, and not in the second. This is handled by
introducing a level of indirection by introducing a new variable similar to the one in
`substitute_derivative_algevars!`. `D(x)` is replaced by `xˍt`, and `D(x) ~ t` is replaced
with `D(x) ~ xˍt` and the equation `xˍt ~ t` is added. Doing this requires multiple steps,
since we're adding a new variable in this case instead of just renaming an existing one.

For each relevant derivative variable, we first check if an equation `D(x) ~ xˍt` already
exists to avoid adding duplicates. Apart from adding the new variables/equations and updating
the matching, we track the `D(x)` variable and the newly added `xˍt`. Note that this doesn't
actually substitute `xˍt` into the required equations. This will be handled in a later part
of the pipeline.

For each `(D(x), xˍt)` pair added, we find the SCC that `D(x)` is in and replace it with `xˍt`
to ensure that the correct variable is used and we generate equations in the correct order.
The singleton SCC consisting of the equation `D(x) ~ xˍt` and variable `xˍt` is added after
the SCC containing `D(x)`. There is an edge case here, where the new SCC must occur before
the SCC for `D(D(x))`. We check for this and appropriately adjust the index. More information
is in the inline comments in the function.

`generate_system_equations!` is the next step, and it does what it says on the tin. This
ends up as a fairly complicated function with several important invariants it tries to maintain.
First off, `total_sub` is a dictionary (`ExpandDerivativeDict`) mapping the LHS each differential
equations to its RHS. `DerivativeDict` is an `AbstractDict` wrapper in MTK which is useful since
nested derivatives (`(D(D(x)))`) are actually just one derivative operator with a higher `.order`.
So substituting `D(x) => y` into `D(D(x))` should turn into `D(y)` but typically it won't.
`DerivativeDict` handles this. `ExpandDerivativeDict` is an extension of this which also runs
`expand_derivatives` on the RHS. For example, if `D(x) => 2y + t` is a mapping then `DerivativeDict`
would return `D(2y + t)` for the key `D(D(x))` but `ExpandDerivativeDict` will return `2D(y) + 1`.
`total_sub` is substituted into every generated equation so that ones dependent on derivatives of
differential variables work correctly. This is what handles the `D(x) ~ xˍt` equations added in
`substitute_derivative_algevars!`.

`generate_system_equations!` generates equations on an SCC-by-SCC basis. This allows us to
generate equations in BLT sorted form. It also offers stronger guarantees that are useful during
debugging. An example is that `total_sub` is very intentionally populated and used for substitution
in parallel. If a bug is found where a derivative in an equation is not substituted by `total_sub`,
this is _never_ an issue with `total_sub`. Instead, it is most likely an SCC ordering or equation
ordering issue. AI tools such as Claude _love_ to blame `total_sub` and suggest substituting it
into all equations at the end. This is not something that should be necessary and is only a
symptomatic fix that hides the true issue.

We loop over all `var_sccs` and codegen them in turn. The first step is to obtain the sorted
`vscc` and `escc` (the latter using `var_eq_matching` and `full_var_eq_matching`). Ordinarily,
each equation in `escc` is generated in turn using the `EquationGenerator` struct. We check
`eq_var_matching` to obtain the corresponding variable (if exists) and hand off to `codegen_equation!`.
Codegen here happens in a very straightforward manner described extensively in the source
code comments. The unique edge case is that when a differential equation `D(x) ~ ...` is generated,
the incidence graph is updated such that each occurrence of `D(x)` is replaced by the incidence
of its RHS. The new symbolic equations are stored in `EquationGenerator`, along with:
- The indices of the equations they came from in `state`.
- The indices of variables they are matched to (0 for algebraic equations).
- The SCC index of each equation.

When `inline_linear_sccs = true`, we check if each SCC is a valid inline linear SCC before
the loop that calls `codegen_equation!`. `get_linear_scc_linsol` returns `nothing` if it
is not a valid inline linear SCC. Otherwise, it returns a 3-tuple of the linear solve
expression, and masks for `escc` and `vscc`. Equations/variables present in the mask are
the ones involved in the linear solve expression. The rest should be generated as
ordinary equations. This allows the inline linear SCCs to be as small as possible.

`get_linear_scc_linsol` considers an SCC valid if it has algebraic equations linear
in the variables of the SCC. Irreducible variable are never involved in the linear solve,
since the entire linear solve ends up as a sequence of observed or differential equations.
The symbolic variables in the SCC are not simply `fullvars[vscc]`. Differential variables
need to be processed as they would when generating the corresponding differential equation.
After doing this, `get_linear_scc_linsol` builds the symbolic `A` and `b` arrays. These matrices
currently involve all equations and variables in the SCC. Typically, a large portion of the SCC
is already torn. We can reuse this information to reduce the size of the linear solve. This
is handled by `__reduce_linear_system!`. It simply eliminates all torn equations in `A`.
The coefficients of the torn variables are back substituted into the remaining equations.
This turns a large sparse linear solve into a small dense solve. The implementation details
of this process are well-documented in the source code.

After this process, we perform the check for `analytical_linear_scc_limit`. SCCs smaller than
this limit are run through a symbolic LU factorization instead of emitting a runtime solve.
Single equation SCCs mandatorily hit this path, since they cannot be symbolically singular.
SCCs below the threshold must also satisfy the `allow_symbolic` and `allow_parameter`
requirements.

All other SCCs are generated as a symbolic expression with `MTKTearing.INLINE_LINEAR_SCC_OP`
as the head. To avoid allocating `A` and `b` in each `f` call, this uses the DiffCache API
in ModelingToolkitBase, enabling non-allocating construction that is robust to automatic
differentiation. `PreallocationTools.get_tmp` requires a "representative" or "reference"
value that defines whether it should return a standard floating point buffer or a buffer of
duals. This is the `reference` construction logic here. The expressions for `A` and `b`
are not `array_literal`. They use `ArrayMaker` with `fill!` to fill the buffer with zeros and
only write the non-zero elements. This can drastically reduce the size of the emitted code
for `prob.f.f`.

The equations in the linear solve stand independently. It is important that they are
generated first prior to the ones masked out and not present in the linear solve. The
latter group can refer to and require `total_sub` information from the former. The
loop that generates linear solve equations has to imitate the internals of
`codegen_equation!` a little. It performs broadly the same operations.

The zeros used for algebraic variables in `codegen_equation!` are now filled in
appropriately to maintain BLT ordering. After `generate_system_equations!`, we
now run `reorder_vars!` which uses the BLT ordering of equations and variables
returned from the previous step to subset `state`. Then, `update_simplified_system!`
propagates these changes to `state.sys`, which is the result of reassembly.

Back in ModelingToolkit.jl's `systemstructure.jl`, the observed equations are
topologically sorted. In ModelingToolkitBase.jl's `systems.jl`, the `additional_passes`
are run and the system passed through `complete`. This concludes `mtkcompile`.

## Clock Inference

Clock inference begins with the `ClockInference` struct. This takes `state` and
builds the data structures required for inference. The core data structure here is a
hypergraph. Nodes are of type `ClockVertex`, a Moshi.jl ADT with variants:
- `Variable`: A variable in the system.
- `Equation`: An equation in the system.
- `InitEquation`: An initialization equation in the system.
- `Clock`: A concrete, known clock in the system.
- `InferredClock`: An unknown clock in the system identified by UUID.
- `Expression`: An arbitrary expression we want to track the clock of.
- `IntegerSequence`: A compatibility node to support ModelingToolkit.jl's simple
  discrete systems.

Hyperedges denote clock equality: all vertices incident on an edge have the same clock.

Clocks are tracked using `ShiftIndex`. `k = ShiftIndex(clk)` creates a `ShiftIndex` with
the concrete clock `clk`. This can then be used to shift variables with syntax

```julia
@variables x(t)
fibonacci_eq = [x(k) ~ x(k-1) + x(k-2)]
```

When a shift is thus applied to a variable, the result is a `Shift(t, n)(x)` where the
`x` has a new metadata entry: `VariableTimeDomain` with the value `clk`. This embeds the
clock in the symbolic structure. For ModelingToolkit.jl's simple discrete systems, the
syntax `k = ShiftIndex(t)` is used where `t` is the independent variable in the system.
This hits the `IntegerSequence` case. When using `k = ShiftIndex()` where no clock is passed,
the clock is inferred. `k` internally gets a placeholder clock with a UUID to distinguish
uses of different inferred clocks.

The clock inference constructor goes over `state.fullvars` and checks if it has a
concrete known clock. The variable is then associated with that clock. This is likely not
strictly necessary due to the behavior of `infer_clocks!` (below). Often the same variable
may or may not have a clock in different equations. For example, a `connect` equation involving
a variable won't have a clock associated with it. The variable present in `fullvars` may or
may not have the `VariableTimeDomain` metadata depending on which equation it was picked up
from.

Every unique shift of a variable is a different entry in `fullvars`: `x(k)` and `x(k-1)`
are present separately. The first step in `infer_clocks!` is to add hyperedges marking
all shifts of the same variable as having the same clock. Inference then proceeds to
process all equations and initialization equations using the `InferEquationClosure`
callable struct. This struct finds all variables and operator applications in the equation.
It then builds an `InferVariableClosure` and calls it on each of the variables/operator calls.

The functioning of these closures is subtle and relies on carefully built aliasing guarantees.
`InferEquationClosure` (among other things) contains fields `hyperedge` and `arg_hyperedge`.
The former is the hyperedge for each equation, emptied between equations. The second is the
hyperedge for each argument to each operator call. It is cleared between operator arguments.
`InferEquationClosure` itself only populates `hyperedge` with the vertex for the equation
it is processing. `InferVariableClosure` contains (among others) the fields:
- `parent_hyperedge`: the hyperedge the variable/operator call currently being process
  exists in.
- `arg_hyperedge`: semantically identical to the same field in `InferEquationClosure`.

`parent_hyperedge` is critical here. Operator calls can be nested, and each argument to
each operator call can be in its own clock domain. When `InferEquationClosure` calls
`InferVariableClosure`, the `parent_hyperedge` of the latter is the `hyperedge` field
of the former since the clock of the variable/operator call must be the clock of the
equation. `InferVariableClosure` itself can recurse, since it might be processing
an operator call where an argument involves another operator call. The `parent_hyperedge`
for the inner operator call is the `arg_hyperedge` for the outer `InferVariableClosure`.

If `InferVariableClosure` gets a variable, it will simply add it to `parent_hyperedge`.
If the variable contains a `VariableTimeDomain` metadata entry, this will also be
added to `parent_hyperedge`. If `InferVariableClosure` gets an operator call, it will
do the same but also process each argument _if_ the operator satisfies the
`is_timevarying_operator` trait. Operators satisfying this trait implement `input_timedomain`
and `output_timedomain`. The former returns an array of clock domains for each of the
arguments, and the latter a single time domain for its output. They are not expeceted to
return concrete domains, but rather to use `InferredDiscrete`. `InferredDiscrete` takes
an integer index. For a given operator call, multiple `InferredDiscrete`s in `input_timedomain`
or `output_timedomain` with the same index have the same unknown clock. For example:

```julia
input_timedomain(::MyOp, _) = MTKTearing.InputTimeDomainElT[InferredDiscrete(1), InferredDiscrete(2)]
output_timedomain(::MyOp, _) = InferredDiscrete(1)
```

Here, this means that the two arguments of `MyOp` have no relation in clock (they may have the
same or different clocks). However, the first argument has the same clock as the output.
We use this information to add additional hyperedges to the graph. This information is tracked
per operator call using the `relative_hyperedges` field on the `Infer*Closure` structs.

Now, `InferVariableClosure` will do something very similar to `InferEquationClosure`. For
every argument of a timevarying operator call, it will use `search_variables!` to identify
all variables/nested operator calls. It will then build a new `InferVariableClosure` and
recurse on each of them. The `parent_hyperedge` for the inner `InferVariableClosure` is the
`arg_hyperedge` for the outer one. The `arg_hyperedge` is emptied before processing each
argument.

It is also often very useful to know the clock of each argument expression of operator
calls. Every argument thus constitutes an `Expression` vertex. An edge case exists to
avoid doing this for constant arguments. `InferEquationClosure` and `InferVariableClosure`
also track `must_be_discrete` - a list of vertices we know must be discrete due to
the presence of `InferredDiscrete`.

This behavior can be a bit confusing, but is made clearer by referring to the source
code and walking through it with a debugger. The test cases relevant to this in
ModelingToolkitTearing.jl are useful here.

After all this inference, the clock partitions are simply the connected components of the
hypergraph. `infer_clocks!` then processes each partition individually. First, it finds
all the `Clock` nodes in a partition. If it has no such `Clock` node and no `IntegerSequence`
node, the partition is marked as continuous. `ContinuousClock` nodes are added by the
`Differential` operator. However, an independent set of algebraic equations may not have
any such operator. Thus, an independent set of equations with no clock is assumed to be
continuous. If this partition contains any `must_be_discrete` nodes, we can infer that
this partition has to be discrete but does not have a discrete clock. This is an error.
If the partition has an explicit `ContinuousClock` and `must_be_discrete_nodes`, this
is also an error. If the partition has multiple clocks, it errors.

After all this validation, we process each variable, equation and init equation in
the system and mark its known time domain in the appropriate field of `ClockInference`
(`var_domain`/`eq_domain`/`init_eq_domain`). We also track information about
`Expression` vertices using the `expression_clocks` field.

## Bibliography

[1] Sven Erik Mattsson and Gustaf Söderlind. 1993. Index Reduction in Differential-Algebraic Equations Using Dummy Derivatives. SIAM J. Sci. Comput. 14, 3 (May 1993), 677–692. https://doi.org/10.1137/0914043

[2] Constantinos C. Pantelides. 1988. The Consistent Initialization of Differential-Algebraic Systems. SIAM J. Sci. Stat. Comput. 9, 2 (March 1988), 213–231. https://doi.org/10.1137/0909014
