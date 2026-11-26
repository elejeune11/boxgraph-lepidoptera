import math
from scipy.optimize import brentq


def run(inputs, settings):
    a = inputs["a"]
    eps_b = inputs["eps_b"]
    rho = inputs["rho"]
    lam = inputs["lam"]
    E_c = inputs["E_c"]
    rtol = float(settings.get("rtol", 1e-10))
    h_rel = float(settings.get("h_rel", 1e-6))

    def energy(R):
        x = R / a
        e_bind = -math.pi * x * eps_b * (x - 2.0) / (x * x)
        e_def = rho * math.pi * ((R - lam) ** 2) * E_c / (x * x)
        return e_bind + e_def

    def dE(R):
        h = R * h_rel
        return (energy(R + h) - energy(R - h)) / (2.0 * h)

    lo = lam * (1.0 + 1e-9)
    if dE(lo) >= 0.0:
        return {"R_0": None}
    hi = lam * 2.0
    while dE(hi) <= 0.0:
        hi *= 2.0
        if hi > lam * 1e14:
            return {"R_0": None}

    R_0 = brentq(dE, lo, hi, rtol=rtol, xtol=rtol * lam)
    return {"R_0": float(R_0)}
