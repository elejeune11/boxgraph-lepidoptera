import math


def run(inputs, settings):
    return {"delta_ell": math.sqrt(inputs["kBT"] / inputs["Kc"])}


def symbolic():
    import sympy as sp

    kBT, Kc = sp.symbols("kBT Kc", positive=True)
    return {"delta_ell": sp.sqrt(kBT / Kc)}
