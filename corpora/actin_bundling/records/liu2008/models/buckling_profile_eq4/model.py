import math


def run(inputs, settings):
    a = float(inputs["a"])  # nm
    L = float(inputs["L"])  # um
    s = float(inputs["s"])  # um
    return {"r_x": a * math.cos(math.pi * s / (2.0 * L)) - a}


def symbolic():
    import sympy as sp

    a, L, s = sp.symbols("a L s", positive=True)
    return {"r_x": a * sp.cos(sp.pi * s / (2 * L)) - a}
