import math


def run(inputs, settings):
    # EI in N*m^2, L in um -> F_b in pN (1e24 conversion factor)
    F_b = (math.pi**2) * inputs["EI"] * 1e24 / (inputs["L"] ** 2)
    return {"F_b": F_b}


def symbolic():
    import sympy as sp

    EI, L = sp.symbols("EI L", positive=True)
    return {"F_b": (sp.pi**2) * EI * sp.Integer(10) ** 24 / (L**2)}
