import math
from scipy.optimize import brentq


def run(inputs, settings):
    a = inputs["a"]
    eps_b = inputs["eps_b"]
    k_b_prime = inputs["k_b_prime"]
    rho = inputs["rho"]
    rtol = float(settings.get("rtol", 1e-10))
    h_rel = float(settings.get("h_rel", 1e-6))

    def energy(R):
        x = R / a
        e_bind = -math.pi * x * eps_b * (x - 2.0) / (x * x)
        e_def = 2.0 * math.pi * k_b_prime * (1.0 / (x * x)) * (R ** 4) * (rho ** 2) * math.log(x)
        return e_bind + e_def

    def dE(R):
        h = R * h_rel
        return (energy(R + h) - energy(R - h)) / (2.0 * h)

    # Note: E_defects_b' changes sign at R/a = exp(-1/2) (~0.6065); for all R > exp(-1/2)*a,
    # dE increases monotonically from -inf to +inf.
    lo = a * math.exp(-0.5) * (1.0 + 1e-9)
    hi = a * 2.0
    while dE(lo) >= 0.0:
        lo = a * math.exp(-0.5) + 0.1 * (lo - a * math.exp(-0.5))
    while dE(hi) <= 0.0:
        hi *= 2.0
        if hi > a * 1e14:
            return {"R_0_over_a": None, "R_0": None}

    R_0 = float(brentq(dE, lo, hi, rtol=rtol, xtol=rtol * a))
    return {"R_0_over_a": R_0 / a, "R_0": R_0}
