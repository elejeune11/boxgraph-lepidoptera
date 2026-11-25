import math


def run(inputs, settings):
    eta = float(inputs["eta"])
    gamma = float(inputs["gamma"])
    return {"lam": math.sqrt(eta / gamma)}


def symbolic():
    import sympy as sp

    eta, gamma = sp.symbols("eta gamma", positive=True)
    return {"lam": sp.sqrt(eta / gamma)}
