module MTKDiffEqNoiseProcess

using ModelingToolkit: ModelingToolkit
using DiffEqNoiseProcess: WienerProcess

ModelingToolkit.scalar_noise() = WienerProcess(0.0, 0.0, 0.0)

end
