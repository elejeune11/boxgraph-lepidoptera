import math


def run(inputs, settings):
    R_0 = inputs["R_0"]
    u = inputs["u"]
    lam = inputs["lam"]
    # Exact integral of 2*pi * int_0^lam (R_0 + u*sin(2*pi*z/lam) - u^2/(4*R_0))^2 dz
    V = 2.0 * math.pi * lam * (R_0**2 + (u**4) / (16.0 * (R_0**2)))
    return {"V": V}


def symbolic():
    import sympy as sp

    R_0, u, lam = sp.symbols("R_0 u lam", positive=True)
    return {"V": 2 * sp.pi * lam * (R_0**2 + u**4 / (16 * R_0**2))}
