import numpy as np


def run(inputs, settings):
    ts = np.array([float(v) for v in inputs["R_curve"]["x"]])
    Rs = np.array([float(v) for v in inputs["R_curve"]["y"]])
    t_lo = float(settings.get("t_lo", 0.1))
    t_hi = float(settings.get("t_hi", 12.0))
    m = (ts >= t_lo) & (ts <= t_hi) & (Rs > 0.0)
    beta = float(np.polyfit(np.log(ts[m]), np.log(Rs[m]), 1)[0])
    return {"coarsening_exponent": beta}
