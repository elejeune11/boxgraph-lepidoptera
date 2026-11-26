import math


def run(inputs, settings):
    c = inputs["c"]
    eps = inputs["eps"]
    ce = c * math.exp(eps)
    return {"eps_b": ce / (1.0 + ce)}


def symbolic():
    import sympy as sp

    c, eps = sp.symbols("c eps", positive=True)
    return {"eps_b": c * sp.exp(eps) / (1 + c * sp.exp(eps))}
