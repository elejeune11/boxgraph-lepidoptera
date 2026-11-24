import math


def run(inputs, settings):
    L_f = inputs["L_f"]
    L_i = inputs["L_i"]
    return {"L": L_f + (L_i - L_f) * math.exp(-inputs["t"] / inputs["tau"])}


def symbolic():
    import sympy as sp

    t, L_i, L_f, tau = sp.symbols("t L_i L_f tau", positive=True)
    return {"L": L_f + (L_i - L_f) * sp.exp(-t / tau)}
