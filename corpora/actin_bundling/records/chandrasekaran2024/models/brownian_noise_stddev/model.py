import math


def run(inputs, settings):
    return {"sigma_B": math.sqrt(2.0 * inputs["D"] * inputs["dt"])}


def symbolic():
    import sympy as sp

    D, dt = sp.symbols("D dt", positive=True)
    return {"sigma_B": sp.sqrt(2 * D * dt)}
