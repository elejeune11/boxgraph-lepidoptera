import math
from scipy.optimize import brentq


def run(inputs, settings):
    t = inputs["t"]
    w1 = inputs["w1"]
    tau1 = inputs["tau1"]
    tau2 = inputs["tau2"]
    w2 = 1.0 - w1

    def survival(x):
        return w1 * math.exp(-x / tau1) + w2 * math.exp(-x / tau2)

    pct = 100.0 * (1.0 - survival(t))
    hi = max(tau1, tau2) * 20.0
    t_50 = float(brentq(lambda x: survival(x) - 0.5, 0.0, hi))
    t_80 = float(brentq(lambda x: survival(x) - 0.2, 0.0, hi))
    return {"percent": pct, "t_50": t_50, "t_80": t_80}


def symbolic():
    import sympy as sp

    t, w1, tau1, tau2 = sp.symbols("t w1 tau1 tau2", positive=True)
    return {
        "percent": 100 * (1 - w1 * sp.exp(-t / tau1) - (1 - w1) * sp.exp(-t / tau2)),
    }
