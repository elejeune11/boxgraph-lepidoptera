import math


def run(inputs, settings):
    return {"delta_ell_B": 2.0 * math.sqrt(inputs["kBT"] / inputs["Kc"])}


def symbolic():
    import sympy as sp

    kBT, Kc = sp.symbols("kBT Kc", positive=True)
    return {"delta_ell_B": 2 * sp.sqrt(kBT / Kc)}
