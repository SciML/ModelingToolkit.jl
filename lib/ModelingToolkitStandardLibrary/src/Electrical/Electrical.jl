"""
Library of electrical models.
This library contains electrical components to build up analog circuits.
"""
module Electrical

using ModelingToolkit, Symbolics, IfElse
using ModelingToolkit: t_nounits as t, D_nounits as D
using ..Thermal: HeatPort
using ..Mechanical.Rotational: Flange, Support
using ..Blocks: RealInput, RealOutput

export Pin, OnePort
include("utils.jl")

export Capacitor,
       Ground, Inductor, Resistor, Conductor, Short, IdealOpAmp, EMF,
       Diode, VariableResistor
include("Analog/ideal_components.jl")

export CurrentSensor, PotentialSensor, VoltageSensor, PowerSensor, MultiSensor
include("Analog/sensors.jl")

export Voltage, Current
include("Analog/sources.jl")

export NMOS, PMOS
include("Analog/mosfets.jl")

export NPN, PNP
include("Analog/transistors.jl")

# include("Digital/gates.jl")
# include("Digital/sources.jl")

# TODO:
# - digital
# - machines
# - multi-phase

export Logic
include("Digital/logic.jl")

export StdLogicVector, StdULogicVector,
       std_ulogic, UX01, UX01Z, X01, X01Z,
       get_logic_level
include("Digital/logic_vectors.jl")

export LogicTable,
       AndTable, OrTable, NotTable, XorTable
include("Digital/tables.jl")

end
