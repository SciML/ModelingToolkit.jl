"""
Library of thermal system components to model heat transfer.
"""
module Thermal
using ModelingToolkit, Symbolics, IfElse
using ModelingToolkit: t_nounits as t, D_nounits as D
using ...Blocks: RealInput, RealOutput

export HeatPort, Element1D
include("utils.jl")

export BodyRadiation, ConvectiveConductor, ConvectiveResistor, HeatCapacitor,
       ThermalConductor,
       ThermalResistor, ThermalCollector
include("HeatTransfer/ideal_components.jl")

export RelativeTemperatureSensor, HeatFlowSensor, TemperatureSensor
include("HeatTransfer/sensors.jl")

export FixedHeatFlow, FixedTemperature, PrescribedHeatFlow, PrescribedTemperature
include("HeatTransfer/sources.jl")

# Simple components for 1-dimensional incompressible thermo-fluid flow models
# TODO:
# - FluidHeatFlow

end
