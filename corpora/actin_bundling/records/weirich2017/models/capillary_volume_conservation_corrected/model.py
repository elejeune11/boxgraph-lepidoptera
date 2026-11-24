import math


def run(inputs, settings):
    R_0 = inputs["R_0"]
    u = inputs["u"]
    lam = inputs["lam"]
    V = math.pi * lam * (R_0**2 + (u**4) / (16.0 * (R_0**2)))
    return {"V": V}


def symbolic():
    import sympy as sp

    R_0, u, lam = sp.symbols("R_0 u lam", positive=True)
    return {"V": sp.pi * lam * (R_0**2 + u**4 / (16 * R_0**2))}
