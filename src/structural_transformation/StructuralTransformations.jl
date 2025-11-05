module StructuralTransformations

using Setfield: @set!, @set
using UnPack: @unpack

using Symbolics: unwrap, linear_expansion, fast_substitute
import Symbolics
using SymbolicUtils
using SymbolicUtils.Code
using SymbolicUtils.Rewriters
using SymbolicUtils: maketerm, iscall

using ModelingToolkit
using ModelingToolkit: System, AbstractSystem, var_from_nested_derivative, Differential,
                       unknowns, equations, vars, Symbolic, diff2term_with_unit,
                       shift2term_with_unit, value,
                       operation, arguments, Sym, Term, simplify, symbolic_linear_solve,
                       isdiffeq, isdifferential, isirreducible,
                       empty_substitutions, get_substitutions,
                       get_tearing_state, get_iv, independent_variables,
                       has_tearing_state, defaults, InvalidSystemException,
                       ExtraEquationsSystemException,
                       ExtraVariablesSystemException,
                       vars!, invalidate_cache!,
                       vars!, invalidate_cache!, Shift,
                       IncrementalCycleTracker, add_edge_checked!, topological_sort,
                       filter_kwargs, lower_varname_with_unit,
                       lower_shift_varname_with_unit, setio, SparseMatrixCLIL,
                       get_fullvars, has_equations, observed,
                       Schedule, schedule, iscomplete, get_schedule

using ModelingToolkit.BipartiteGraphs
import .BipartiteGraphs: invview, complete
import ModelingToolkit: var_derivative!, var_derivative_graph!
using Graphs
using ModelingToolkit: algeqs, EquationsView,
                       SystemStructure, TransformationState, TearingState,
                       mtkcompile!,
                       isdiffvar, isdervar, isalgvar, isdiffeq, algeqs, is_only_discrete,
                       dervars_range, diffvars_range, algvars_range,
                       DiffGraph, complete!,
                       get_fullvars, system_subset
using SymbolicIndexingInterface: symbolic_type, ArraySymbolic, NotSymbolic

using ModelingToolkit.DiffEqBase
using ModelingToolkit.StaticArrays
using RuntimeGeneratedFunctions: @RuntimeGeneratedFunction,
                                 RuntimeGeneratedFunctions,
                                 drop_expr

RuntimeGeneratedFunctions.init(@__MODULE__)

using SparseArrays

using SimpleNonlinearSolve

using DocStringExtensions

export tearing, dae_index_lowering, check_consistency
export dummy_derivative
export sorted_incidence_matrix,
       pantelides!, pantelides_reassemble, find_solvables!,
       linear_subsys_adjmat!
export tearing_substitution
export torn_system_jacobian_sparsity
export full_equations
export but_ordered_incidence, lowest_order_variable_mask, highest_order_variable_mask
export computed_highest_diff_variables
export shift2term, lower_shift_varname, simplify_shifts, distribute_shift

include("utils.jl")
include("tearing.jl")
include("pantelides.jl")
include("bipartite_tearing/modia_tearing.jl")
include("symbolics_tearing.jl")
include("partial_state_selection.jl")
include("codegen.jl")

end # module
