import math


def run(inputs, settings):
    tau_s = float(inputs["tau"]) * 1e-3
    ref = inputs["ref_curve"]
    grid = inputs["x_grid"]
    # Extract initial aspect ratio AR_0 at earliest time point (t = 0) of reference fit curve
    idx0 = min(range(len(ref["x"])), key=lambda i: abs(float(ref["x"][i])))
    t0 = float(ref["x"][idx0])
    ar0 = float(ref["y"][idx0])
    amp = ar0 - 1.0
    xs = [float(t) for t in grid["x"]]
    ys = [1.0 + amp * math.exp(-(t - t0) / tau_s) for t in xs]
    return {"AR": {"x": xs, "y": ys}}
