module MTKDiffEqNoiseProcessExt

using DiffEqNoiseProcess: WienerProcess
import ModelingToolkit

function ModelingToolkit.__default_wiener_process(::Nothing)
    return WienerProcess(0.0, 0.0, 0.0)
end

end
