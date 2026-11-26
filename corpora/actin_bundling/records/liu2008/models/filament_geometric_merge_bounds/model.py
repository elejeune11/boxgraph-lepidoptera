import math


def run(inputs, settings):
    L = float(inputs["L"])
    return {
        "x_max": 2.0 * L / math.pi,
        "D_max": 4.0 * L / math.pi,
    }


def symbolic():
    import sympy as sp

    L = sp.symbols("L", positive=True)
    return {
        "x_max": 2 * L / sp.pi,
        "D_max": 4 * L / sp.pi,
    }
