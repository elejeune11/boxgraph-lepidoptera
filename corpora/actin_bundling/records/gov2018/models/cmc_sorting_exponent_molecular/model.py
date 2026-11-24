import math


def run(inputs, settings):
    chi = inputs["kappa_over_kBT"] * inputs["a_over_R"]
    return {"chi": chi, "enrichment": math.exp(chi)}


def symbolic():
    import sympy as sp

    kappa_over_kBT, a_over_R = sp.symbols("kappa_over_kBT a_over_R", positive=True)
    chi = kappa_over_kBT * a_over_R
    return {"chi": chi, "enrichment": sp.exp(chi)}
