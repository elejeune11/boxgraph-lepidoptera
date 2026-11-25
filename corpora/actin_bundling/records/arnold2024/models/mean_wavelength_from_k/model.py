import math


def run(inputs, settings):
    return {"wavelength": 2.0 * math.pi / inputs["k_mean"]}


def symbolic():
    import sympy as sp

    k_mean = sp.symbols("k_mean", positive=True)
    return {"wavelength": 2 * sp.pi / k_mean}
