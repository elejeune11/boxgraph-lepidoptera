import math


def run(inputs, settings):
    xi = inputs["xi"]
    eta = inputs["eta"]
    y = inputs["y"]
    sxi = math.sin(xi)
    sc = sxi * math.cos(eta)
    term1 = sc + (1.0 + sc) * (y * y - 1.0) / 2.0
    term2 = (y * sxi * math.sin(eta)) ** 2
    return {"N": term1 * term1 + term2}


def symbolic():
    import sympy as sp

    xi, eta, y = sp.symbols("xi eta y", positive=True)
    sc = sp.sin(xi) * sp.cos(eta)
    N = (sc + (1 + sc) * (y**2 - 1) / 2) ** 2 + y**2 * sp.sin(xi) ** 2 * sp.sin(eta) ** 2
    return {"N": N}
