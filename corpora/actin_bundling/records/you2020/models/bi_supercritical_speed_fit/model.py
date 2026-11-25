import math


def run(inputs, settings):
    k = float(inputs["kappa"])
    d = float(inputs["delta"])
    dc = float(inputs["delta_c_over_kappa"]) * k
    v = math.sqrt(d * d - dc * dc) if d > dc else 0.0
    return {"v": v}


def symbolic():
    import sympy as sp

    kappa, delta, delta_c_over_kappa = sp.symbols("kappa delta delta_c_over_kappa", positive=True)
    return {"v": sp.sqrt(delta**2 - (delta_c_over_kappa * kappa) ** 2)}
