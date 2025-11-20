module MTKDiffEqNoiseProcessExt

using DiffEqNoiseProcess: WienerProcess
import ModelingToolkitBase

function ModelingToolkitBase.__default_wiener_process(::Nothing)
    return WienerProcess(0.0, 0.0, 0.0)
end

end
