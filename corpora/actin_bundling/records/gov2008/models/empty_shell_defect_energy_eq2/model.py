import math


def run(inputs, settings):
    R = inputs["R"]
    a = inputs["a"]
    rho = inputs["rho"]
    lam = inputs["lam"]
    E_c = inputs["E_c"]
    return {"E_defect": rho * math.pi * ((R - lam) ** 2) * E_c / ((R / a) ** 2)}


def symbolic():
    import sympy as sp

    R, a, rho, lam, E_c = sp.symbols("R a rho lam E_c", positive=True)
    return {"E_defect": rho * sp.pi * (R - lam) ** 2 * E_c / (R / a) ** 2}
