import math


def run(inputs, settings):
    k = 2.0 * math.pi * inputs["n"] / inputs["L"]
    wavelength = inputs["L"] / inputs["n"]
    return {"k": k, "wavelength": wavelength}


def symbolic():
    import sympy as sp

    n, L = sp.symbols("n L", positive=True)
    return {"k": 2 * sp.pi * n / L, "wavelength": L / n}
