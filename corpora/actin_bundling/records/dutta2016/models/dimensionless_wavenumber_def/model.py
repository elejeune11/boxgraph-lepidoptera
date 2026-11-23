import math


def run(inputs, settings):
    return {"k_star": inputs["k"] * math.sqrt(inputs["L"] * inputs["kT"] / inputs["epsilon"])}


def symbolic():
    import sympy as sp
    k, L, kT, epsilon = sp.symbols("k L kT epsilon", positive=True)
    return {"k_star": k * sp.sqrt(L * kT / epsilon)}
