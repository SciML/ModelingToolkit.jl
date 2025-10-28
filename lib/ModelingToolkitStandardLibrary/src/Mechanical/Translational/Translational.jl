"""
Library to model 1-dimensional, translational mechanical systems
"""
module Translational

using ModelingToolkit, Symbolics
using ModelingToolkit: getdefault, t_nounits as t, D_nounits as D

using ModelingToolkitStandardLibrary.Blocks: RealInput, RealOutput
using IfElse: ifelse

export MechanicalPort
include("utils.jl")

export Mass, Spring, Damper, Fixed
include("components.jl")

export Force, Position, Velocity, Acceleration
include("sources.jl")

export ForceSensor, PositionSensor, AccelerationSensor
include("sensors.jl")

end
