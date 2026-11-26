import math
from scipy.optimize import brentq


def run(inputs, settings):
    k1 = float(inputs["k1"])
    k2 = float(inputs["k2"])
    xtol = float(settings.get("xtol", 1.0e-12))

    def p_adp_minus_half(t):
        if abs(k1 - k2) < 1e-12:
            p_adp = 1.0 - (1.0 + k1 * t) * math.exp(-k1 * t)
        else:
            p_adp = 1.0 - (k1 * math.exp(-k2 * t) - k2 * math.exp(-k1 * t)) / (k1 - k2)
        return p_adp - 0.5

    hi = max(10.0 / min(k1, k2), 1.0)
    while p_adp_minus_half(hi) <= 0.0:
        hi *= 2.0
    t_half = brentq(p_adp_minus_half, 0.0, hi, xtol=xtol)
    return {"t_half": float(t_half)}
