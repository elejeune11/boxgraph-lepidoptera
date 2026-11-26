import math


def run(inputs, settings):
    a = inputs["a"]
    a_0 = inputs["a_0"]
    return {"alpha_b": math.pi * (a * a - a_0 * a_0)}


def symbolic():
    import sympy as sp

    a, a_0 = sp.symbols("a a_0", positive=True)
    return {"alpha_b": sp.pi * (a**2 - a_0**2)}
