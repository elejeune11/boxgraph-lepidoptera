import math


def run(inputs, settings):
    return {"k_l": inputs["beta"] * inputs["n"] * math.pi / inputs["L_over_l"]}


def symbolic():
    import sympy as sp

    beta, n, L_over_l = sp.symbols("beta n L_over_l", positive=True)
    return {"k_l": beta * n * sp.pi / L_over_l}
