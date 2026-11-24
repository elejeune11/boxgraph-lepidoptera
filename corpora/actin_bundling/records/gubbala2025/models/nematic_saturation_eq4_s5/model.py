import math


def run(inputs, settings):
    beta = float(inputs["beta"])
    phi = float(inputs["phi"])
    A = 1.0 + phi
    S = math.sqrt(A / (2.0 * beta)) if A > 0.0 else 0.0
    return {"S": S}


def symbolic():
    import sympy as sp

    beta = sp.symbols("beta", positive=True)
    phi = sp.symbols("phi", real=True)
    return {"S": sp.sqrt((1 + phi) / (2 * beta))}
