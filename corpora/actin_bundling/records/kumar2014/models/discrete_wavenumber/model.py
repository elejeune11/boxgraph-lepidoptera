import math


def run(inputs, settings):
    return {"k_l": inputs["omega"] * inputs["n"] * math.pi / inputs["L_over_l"]}


def symbolic():
    import sympy as sp

    omega, n, L_over_l = sp.symbols("omega n L_over_l", positive=True)
    return {"k_l": omega * n * sp.pi / L_over_l}
