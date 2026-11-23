import math


def run(inputs, settings):
    return {"l": math.sqrt(inputs["eta"] / inputs["gamma"])}


def symbolic():
    import sympy as sp

    eta, gamma = sp.symbols("eta gamma", positive=True)
    return {"l": sp.sqrt(eta / gamma)}
