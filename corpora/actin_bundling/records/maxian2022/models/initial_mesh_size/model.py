import math


def run(inputs, settings):
    return {"ell_m": math.sqrt((inputs["Ld"] ** 3) / (inputs["F"] * inputs["L"]))}


def symbolic():
    import sympy as sp

    Ld, F, L = sp.symbols("Ld F L", positive=True)
    return {"ell_m": sp.sqrt(Ld**3 / (F * L))}
