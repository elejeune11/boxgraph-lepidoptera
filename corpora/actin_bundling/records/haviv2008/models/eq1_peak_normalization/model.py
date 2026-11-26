import math


def run(inputs, settings):
    f = float(inputs["f"])
    kon_D = float(inputs["kon_over_D"])
    koff_D = float(inputs["koff_over_D"])
    R_n = f * math.exp(-(kon_D ** 2) / koff_D)
    return {"R_n": R_n}


def symbolic():
    import sympy as sp

    f, kon_over_D, koff_over_D = sp.symbols("f kon_over_D koff_over_D", positive=True)
    return {"R_n": f * sp.exp(-(kon_over_D ** 2) / koff_over_D)}
