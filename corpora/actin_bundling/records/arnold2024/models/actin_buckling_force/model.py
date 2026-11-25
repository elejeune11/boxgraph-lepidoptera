import math


def run(inputs, settings):
    # EI in N*m^2, L in um -> L_m = L * 1e-6 m -> EI / L_m^2 in N -> * 1e12 for pN
    # Combined factor: 1e12 / (1e-6)^2 = 1e24
    F_b = math.pi * inputs["EI"] * 1e24 / (inputs["L"] ** 2)
    return {"F_b": F_b}


def symbolic():
    import sympy as sp

    EI, L = sp.symbols("EI L", positive=True)
    return {"F_b": sp.pi * EI * sp.Integer(10) ** 24 / (L**2)}
