import math


def run(inputs, settings):
    L = inputs["L"]
    A = inputs["A"]
    lam = inputs["lam"]
    L0 = 2.0 * math.sqrt(math.pi * A)
    L_ex = L / L0 - 1.0
    dF = lam * (L - L0)
    return {"L0": L0, "L_ex": L_ex, "dF": dF}


def symbolic():
    import sympy as sp

    L, A, lam = sp.symbols("L A lam", positive=True)
    L0 = 2 * sp.sqrt(sp.pi * A)
    L_ex = L / L0 - 1
    dF = lam * (L - L0)
    return {"L0": L0, "L_ex": L_ex, "dF": dF}
