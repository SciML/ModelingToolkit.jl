module StructuralTransformations

using Setfield: @set!, @set
using UnPack: @unpack

using Symbolics: unwrap, linear_expansion, VartypeT, SymbolicT
import Symbolics
using SymbolicUtils
using SymbolicUtils: BSImpl
using SymbolicUtils.Code
using SymbolicUtils.Rewriters
using SymbolicUtils: maketerm, iscall, symtype
import SymbolicUtils as SU
import Moshi

using ModelingToolkit
using ModelingToolkitBase: System, AbstractSystem, var_from_nested_derivative, Differential,
                       unknowns, equations, diff2term_with_unit,
                       value,
                       operation, arguments, simplify, symbolic_linear_solve,
                       isdiffeq, isdifferential, isirreducible,
                       empty_substitutions, get_substitutions,
                       get_tearing_state, get_iv, independent_variables,
                       has_tearing_state, InvalidSystemException,
                       ExtraEquationsSystemException,
                       ExtraVariablesSystemException,
                       invalidate_cache!, Shift,
                       topological_sort,
                       filter_kwargs, lower_varname_with_unit,
                       setio,
                       has_equations, observed,
                       Schedule, schedule, iscomplete, get_schedule, VariableUnshifted,
                       VariableShift, DerivativeDict, shift2term, simplify_shifts,
                       distribute_shift

using BipartiteGraphs
import BipartiteGraphs: invview, complete, IncrementalCycleTracker, add_edge_checked!
using Graphs
using ModelingToolkit: mtkcompile!
using SymbolicIndexingInterface: symbolic_type, ArraySymbolic, NotSymbolic, getname

using ModelingToolkit.DiffEqBase
using ModelingToolkit.StaticArrays
import Symbolics: Num, Arr, CallAndWrap
import CommonSolve

using SparseArrays

using SimpleNonlinearSolve

using DocStringExtensions

import ModelingToolkitBase as MTKBase
import StateSelection
import StateSelection: CLIL, SelectedState
import ModelingToolkitTearing as MTKTearing
using ModelingToolkitTearing: TearingState, SystemStructure, ReassembleAlgorithm,
                              DefaultReassembleAlgorithm

export tearing, dae_index_lowering
export dummy_derivative
export sorted_incidence_matrix, pantelides_reassemble, find_solvables!
export tearing_substitution
export but_ordered_incidence, lowest_order_variable_mask, highest_order_variable_mask

include("utils.jl")
include("pantelides.jl")

function tearing_substitution(sys::AbstractSystem; kwargs...)
    neweqs = full_equations(sys::AbstractSystem; kwargs...)
    @set! sys.eqs = neweqs
    # @set! sys.substitutions = nothing
    @set! sys.schedule = nothing
end

include("symbolics_tearing.jl")

end # module
