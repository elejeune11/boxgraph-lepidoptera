from scipy.optimize import brentq


def _W_sym(V_ratio):
    x = (0.5 * V_ratio) ** (1.0 / 3.0)
    return 2.0 * ((1.0 + x * x) ** 2.5) / (1.0 + 4.0 * x * x)


def run(inputs, settings):
    W = float(inputs["W"])
    tol = float(settings.get("tol", 1.0e-10))
    # W_sym(V_ratio) is strictly increasing for V_ratio > 0.25 (x > 0.5)
    lo = 0.5
    hi = 10.0
    while _W_sym(hi) < W:
        hi *= 2.0
    V_ratio_sym = float(brentq(lambda v: _W_sym(v) - W, lo, hi, xtol=tol, rtol=tol))
    return {"V_ratio_sym": V_ratio_sym}
