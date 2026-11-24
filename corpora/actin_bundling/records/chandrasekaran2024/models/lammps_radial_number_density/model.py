import math


def run(inputs, settings):
    return {"rho": inputs["N"] / (4.0 * math.pi * (inputs["r"] ** 2) * inputs["dr"])}


def symbolic():
    import sympy as sp

    N, r, dr = sp.symbols("N r dr", positive=True)
    return {"rho": N / (4 * sp.pi * r**2 * dr)}
