import math


def run(inputs, settings):
    x = inputs["x"]
    x0 = inputs["x0"]
    return {
        "survival": x0 * math.exp(-x / inputs["tau_1"])
        + (1.0 - x0) * math.exp(-x / inputs["tau_2"])
    }


def symbolic():
    import sympy as sp

    x, x0, tau_1, tau_2 = sp.symbols("x x0 tau_1 tau_2", positive=True)
    return {"survival": x0 * sp.exp(-x / tau_1) + (1 - x0) * sp.exp(-x / tau_2)}
