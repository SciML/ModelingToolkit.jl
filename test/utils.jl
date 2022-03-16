relerror(truth, pred) =
    maximum(abs.((truth .- pred)) ./ maximum(abs.(truth), dims = 2), dims = 2)[:, 1]

test_relerror(truth, pred, rerr) =
    all(relerror(truth, pred) .<= rerr)
