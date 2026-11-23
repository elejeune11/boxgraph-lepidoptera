import math


def run(inputs, settings):
    x = inputs["x"]
    a = inputs["a"]
    L = inputs["L"]
    return {"w": a * (1.0 + math.cos(2.0 * math.pi * x / L))}


def symbolic():
    import sympy as sp
    x, a, L = sp.symbols("x a L", positive=True)
    return {"w": a * (1 + sp.cos(2 * sp.pi * x / L))}
