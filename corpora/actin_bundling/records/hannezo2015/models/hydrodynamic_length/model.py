import math


def run(inputs, settings):
    # eta [Pa*s*m] / xi [Pa*s/m] -> m^2; sqrt -> m; * 1e6 -> um
    return {"l_h": math.sqrt(inputs["eta"] / inputs["xi"]) * 1e6}


def symbolic():
    import sympy as sp

    eta, xi = sp.symbols("eta xi", positive=True)
    return {"l_h": sp.sqrt(eta / xi) * 10**6}
