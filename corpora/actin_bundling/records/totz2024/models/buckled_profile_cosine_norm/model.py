import math


def run(inputs, settings):
    x_norm = inputs["x_norm"]
    a = inputs["a"]
    L = inputs["L"]
    L_0 = inputs["L_0"]
    return {"w": 0.5 * a * (1.0 + math.cos(2.0 * math.pi * x_norm * L_0 / L))}


def symbolic():
    import sympy as sp
    x_norm, a, L, L_0 = sp.symbols("x_norm a L L_0", positive=True)
    return {"w": (a / 2) * (1 + sp.cos(2 * sp.pi * x_norm * L_0 / L))}
