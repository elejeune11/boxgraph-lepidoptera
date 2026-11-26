import math


def run(inputs, settings):
    return {"ell": math.sqrt(inputs["eta"] / inputs["gamma"])}


def symbolic():
    import sympy as sp

    eta, gamma = sp.symbols("eta gamma", positive=True)
    return {"ell": sp.sqrt(eta / gamma)}
