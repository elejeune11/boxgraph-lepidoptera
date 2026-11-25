import math
from scipy.optimize import brentq


def _W_unequal(V_ratio):
    if abs(V_ratio - 2.0) < 1e-10:
        return 8.0 * math.sqrt(2.0) / 5.0
    x2 = (V_ratio - 1.0) ** (1.0 / 3.0)
    num = 2.0 * (1.0 - 1.0 / x2)
    den = 1.0 / (2.0 * math.sqrt(2.0)) - 1.0 / (x2 * ((1.0 + x2 * x2) ** 1.5))
    return num / den


def run(inputs, settings):
    W = float(inputs["W"])
    tol = float(settings.get("tol", 1.0e-10))
    # W_unequal strictly increases from 8*sqrt(2)/5 (~2.2627) at V_ratio=2 to 4*sqrt(2) (~5.6569) as V_ratio->inf
    lo = 2.0 + 1e-9
    hi = 10.0
    while _W_unequal(hi) < W:
        hi *= 2.0
    V_ratio_unequal = float(brentq(lambda v: _W_unequal(v) - W, lo, hi, xtol=tol, rtol=tol))
    return {"V_ratio_unequal": V_ratio_unequal}
