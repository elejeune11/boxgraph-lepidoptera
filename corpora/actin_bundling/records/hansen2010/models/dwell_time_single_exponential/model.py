import math


def run(inputs, settings):
    return {"survival": inputs["x0"] * math.exp(-inputs["x"] / inputs["tau_1"])}


def symbolic():
    import sympy as sp

    x, x0, tau_1 = sp.symbols("x x0 tau_1", positive=True)
    return {"survival": x0 * sp.exp(-x / tau_1)}
