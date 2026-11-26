import math


def run(inputs, settings):
    beta = inputs["beta"]
    E_total = inputs["E_total"]
    E_0 = inputs["E_0"]
    D_R = inputs["D_R"]
    return {"P": math.exp(-2.0 * beta * (E_total - E_0) / D_R)}


def symbolic():
    import sympy as sp

    beta, E_total, E_0, D_R = sp.symbols("beta E_total E_0 D_R", positive=True)
    return {"P": sp.exp(-2 * beta * (E_total - E_0) / D_R)}
