import math


def run(inputs, settings):
    D_B = inputs["D_B"]
    k_f = inputs["k_f"]
    k_b = inputs["k_b"]
    return {"R_bar": math.sqrt(3.0 * D_B * k_f / ((k_f + k_b) * k_b))}


def symbolic():
    import sympy as sp

    D_B, k_f, k_b = sp.symbols("D_B k_f k_b", positive=True)
    return {"R_bar": sp.sqrt(3 * D_B * k_f / ((k_f + k_b) * k_b))}
