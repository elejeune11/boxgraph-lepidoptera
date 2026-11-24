import math


def run(inputs, settings):
    w = inputs["omega"]
    if w >= 1.0:
        ar = 2.0 * math.sqrt(w)
    else:
        ar = 1.0 + w
    return {"aspect_ratio": ar}


def symbolic():
    import sympy as sp

    omega = sp.Symbol("omega", positive=True)
    return {"aspect_ratio": sp.Piecewise((2 * sp.sqrt(omega), omega >= 1), (1 + omega, True))}
