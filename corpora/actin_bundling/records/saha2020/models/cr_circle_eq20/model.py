import math


def run(inputs, settings):
    c1 = float(inputs["c1"])
    c2 = float(inputs["c2"])
    return {"radius": math.sqrt((c1 * c1 + c2 * c2) / 3.0)}


def symbolic():
    import sympy as sp

    c1, c2 = sp.symbols("c1 c2", real=True)
    return {"radius": sp.sqrt((c1**2 + c2**2) / 3)}
