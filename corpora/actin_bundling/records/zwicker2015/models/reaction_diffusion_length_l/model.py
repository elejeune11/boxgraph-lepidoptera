import math


def run(inputs, settings):
    return {"l": math.sqrt(inputs["D_B"] / (inputs["k_f"] + inputs["k_b"]))}


def symbolic():
    import sympy as sp

    D_B, k_f, k_b = sp.symbols("D_B k_f k_b", positive=True)
    return {"l": sp.sqrt(D_B / (k_f + k_b))}
