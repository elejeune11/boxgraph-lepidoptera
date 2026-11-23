import math


def run(inputs, settings):
    return {"R_g": math.sqrt(inputs["L"] * inputs["kT"] / inputs["epsilon"])}


def symbolic():
    import sympy as sp
    L, kT, epsilon = sp.symbols("L kT epsilon", positive=True)
    return {"R_g": sp.sqrt(L * kT / epsilon)}
