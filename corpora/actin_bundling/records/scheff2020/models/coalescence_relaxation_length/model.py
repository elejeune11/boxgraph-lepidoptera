import math


def run(inputs, settings):
    return {"delta_l": math.sqrt(inputs["A_i"]) - math.sqrt(inputs["A_f"])}


def symbolic():
    import sympy as sp

    A_i, A_f = sp.symbols("A_i A_f", positive=True)
    return {"delta_l": sp.sqrt(A_i) - sp.sqrt(A_f)}
