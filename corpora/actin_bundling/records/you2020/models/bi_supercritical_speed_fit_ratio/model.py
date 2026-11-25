import math


def run(inputs, settings):
    k = float(inputs["kappa"])
    d = float(inputs["delta_over_kappa"]) * k
    dc = float(inputs["delta_c_over_kappa"]) * k
    v = math.sqrt(d * d - dc * dc) if d > dc else 0.0
    return {"v": v}


def symbolic():
    import sympy as sp

    kappa, delta_over_kappa, delta_c_over_kappa = sp.symbols(
        "kappa delta_over_kappa delta_c_over_kappa", positive=True
    )
    return {"v": kappa * sp.sqrt(delta_over_kappa**2 - delta_c_over_kappa**2)}
