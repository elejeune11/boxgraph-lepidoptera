import math


def run(inputs, settings):
    xi = float(inputs["xi"])
    D = float(inputs["D"])
    return {"chi_rho0": (1.0 + math.sqrt(D * xi)) ** 2}


def symbolic():
    import sympy as sp

    xi, D = sp.symbols("xi D", positive=True)
    return {"chi_rho0": (1 + sp.sqrt(D * xi)) ** 2}
