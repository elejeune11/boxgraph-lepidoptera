import math


def run(inputs, settings):
    return {"U": inputs["zdm0"] / math.sqrt(inputs["eta"] * inputs["gamma"])}


def symbolic():
    import sympy as sp

    zdm0, eta, gamma = sp.symbols("zdm0 eta gamma", positive=True)
    return {"U": zdm0 / sp.sqrt(eta * gamma)}
