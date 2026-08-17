"""
    StructuralTransformations

Developer-facing structural transformations used by ModelingToolkit's compiler and SciML
extension packages. The exported transformation functions are versioned developer API;
end-user applications should use [`ModelingToolkitBase.mtkcompile`](@ref) instead.
"""
module StructuralTransformations

using Setfield: @set!
using UnPack: @unpack

using Symbolics: SymbolicT
import Symbolics
import SymbolicUtils
using SymbolicUtils: BSImpl
import SymbolicUtils as SU
import Moshi

import ModelingToolkit
using ModelingToolkitBase: System, AbstractSystem, Differential,
    Equation, equations, full_equations, diff2term_with_unit,
    operation, arguments,
    isdiffeq, isdifferential,
    get_tearing_state, get_iv,
    invalidate_cache!,
    iscomplete, get_schedule

using SymbolicUtils: substitute

using BipartiteGraphs: maximal_matching, ndsts, unassigned, 𝑠neighbors
import BipartiteGraphs: complete
import Graphs
using Graphs: edges, inneighbors, nv, outneighbors
using Graphs.LinAlg: incidence_matrix
import CommonSolve

using SparseArrays: sparse

import DocStringExtensions
using DocStringExtensions: TYPEDSIGNATURES

import ModelingToolkitBase as MTKBase
import StateSelection
import StateSelection: find_solvables!
import ModelingToolkitTearing as MTKTearing
using ModelingToolkitTearing: TearingState, ReassembleAlgorithm,
    DefaultReassembleAlgorithm

export tearing, dae_index_lowering
export dummy_derivative
export sorted_incidence_matrix, pantelides_reassemble, find_solvables!
export tearing_substitution
export but_ordered_incidence, lowest_order_variable_mask, highest_order_variable_mask

include("utils.jl")
include("pantelides.jl")

"""
    tearing_substitution(sys::AbstractSystem; kwargs...)

Replace the equations of `sys` with its fully substituted equations.

This is a structural-transformation helper used by simplification passes. End-user code
should usually call [`ModelingToolkitBase.mtkcompile`](@ref).

# Arguments

- `sys`: system whose equations should be substituted.
- `kwargs...`: keyword arguments forwarded to `full_equations`.

# Returns

A copy of `sys` with substituted equations and no cached schedule.
"""
function tearing_substitution(sys::AbstractSystem; kwargs...)
    neweqs = full_equations(sys::AbstractSystem; kwargs...)
    @set! sys.eqs = neweqs
    # @set! sys.substitutions = nothing
    return @set! sys.schedule = nothing
end

include("symbolics_tearing.jl")

end # module
