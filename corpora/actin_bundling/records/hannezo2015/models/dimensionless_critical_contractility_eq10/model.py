import math


def run(inputs, settings):
    return {"chi_c": (1.0 + math.sqrt(inputs["phi"])) ** 2}


def symbolic():
    import sympy as sp

    phi = sp.Symbol("phi", positive=True)
    return {"chi_c": (1 + sp.sqrt(phi)) ** 2}
