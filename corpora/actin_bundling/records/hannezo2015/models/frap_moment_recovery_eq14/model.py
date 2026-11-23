import math


def run(inputs, settings):
    return {
        "m0_norm": math.exp(-inputs["t"] / inputs["tau"]),
        "variance": 2.0 * inputs["D"] * inputs["t"] + inputs["var0"],
    }


def symbolic():
    import sympy as sp

    t, tau, D, var0 = sp.symbols("t tau D var0", positive=True)
    return {
        "m0_norm": sp.exp(-t / tau),
        "variance": 2 * D * t + var0,
    }
