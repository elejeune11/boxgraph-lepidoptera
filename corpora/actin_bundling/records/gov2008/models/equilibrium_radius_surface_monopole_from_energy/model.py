import math
from scipy.optimize import brentq


def run(inputs, settings):
    a = inputs["a"]
    eps_b = inputs["eps_b"]
    k_s_prime = inputs["k_s_prime"]
    rho = inputs["rho"]
    rtol = float(settings.get("rtol", 1e-10))
    h_rel = float(settings.get("h_rel", 1e-6))

    def energy(R):
        x = R / a
        e_bind = -math.pi * x * eps_b * (x - 2.0) / (x * x)
        e_def = 2.0 * math.pi * k_s_prime * (1.0 / (x * x)) * (R * R) * rho * math.log(x)
        return e_bind + e_def

    def dE(R):
        h = R * h_rel
        return (energy(R + h) - energy(R - h)) / (2.0 * h)

    lo = a * 1e-4
    hi = a * 2.0
    while dE(lo) >= 0.0:
        lo *= 0.1
        if lo < a * 1e-14:
            return {"R_0_over_a": None, "R_0": None}
    while dE(hi) <= 0.0:
        hi *= 2.0
        if hi > a * 1e14:
            return {"R_0_over_a": None, "R_0": None}

    R_0 = float(brentq(dE, lo, hi, rtol=rtol, xtol=rtol * a))
    return {"R_0_over_a": R_0 / a, "R_0": R_0}
