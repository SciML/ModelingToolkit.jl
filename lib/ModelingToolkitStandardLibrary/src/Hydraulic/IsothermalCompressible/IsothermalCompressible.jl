"""
Library to model iso-thermal compressible liquid fluid flow
"""
module IsothermalCompressible

using ModelingToolkit, Symbolics
using ModelingToolkit: t_nounits as t, D_nounits as D

using ...Blocks: RealInput, RealOutput
using ...Mechanical.Translational: MechanicalPort, Mass

using IfElse: ifelse

export HydraulicPort, HydraulicFluid
include("utils.jl")

export Cap, Tube, FixedVolume, DynamicVolume, Open, FlowDivider, Valve, Volume, SpoolValve,
       SpoolValve2Way, Actuator
include("components.jl")

export MassFlow, Pressure, FixedPressure
include("sources.jl")

end
