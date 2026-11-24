import math


def run(inputs, settings):
    n_steps = max(0, int(math.floor((inputs["t"] - inputs["t0"]) / inputs["dt_relax"])))
    # dz is increment in Z-axis span 2a (in nm -> 1e-3 um), so semi-axis a increases by 0.5*dz per step
    a = inputs["R_drop"] + 0.5 * n_steps * (inputs["dz"] * 1e-3)
    ar = (a / inputs["R_drop"]) ** 1.5
    return {"a": a, "aspect_ratio": ar}
