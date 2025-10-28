"""
The module `Blocks` contains common input-output components, referred to as blocks.
"""
module Blocks
using ModelingToolkit, Symbolics
import IfElse: ifelse
import ..@symcheck
using ModelingToolkit: getdefault, t_nounits as t, D_nounits as D

export RealInput, RealInputArray, RealOutput, RealOutputArray, SISO
include("utils.jl")

export Gain, Sum, MatrixGain, Feedback, Add, Add3, Product, Division, Power, Modulo,
       UnaryMinus, Floor, Ceil
export Abs, Sign, Sqrt, Sin, Cos, Tan, Asin, Acos, Atan, Atan2, Sinh, Cosh, Tanh, Exp
export Log, Log10
include("math.jl")

export Constant, TimeVaryingFunction, Sine, Cosine, ContinuousClock, Ramp, Step, ExpSine,
       Square, Triangular, Parameter, SampledData,
       Interpolation, ParametrizedInterpolation
include("sources.jl")

export Limiter, DeadZone, SlewRateLimiter
include("nonlinear.jl")

export Integrator, Derivative, FirstOrder, SecondOrder, StateSpace, TransferFunction
export PI, LimPI, PID, LimPID
include("continuous.jl")

export AnalysisPoint, get_sensitivity, get_comp_sensitivity,
       get_looptransfer, open_loop

end
