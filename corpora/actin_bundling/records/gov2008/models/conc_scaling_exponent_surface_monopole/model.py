import math


def run(inputs, settings):
    c = inputs["c"]
    eps = inputs["eps"]
    # R_{0,s} propto eps_b propto c*exp(eps)/(1 + c*exp(eps)),
    # so d(ln R_{0,s})/d(ln c) = 1 / (1 + c*exp(eps)), independent of unknown prefactor.
    return {"exp_c": 1.0 / (1.0 + c * math.exp(eps))}


def symbolic():
    import sympy as sp

    c, eps = sp.symbols("c eps", positive=True)
    eps_b = c * sp.exp(eps) / (1 + c * sp.exp(eps))
    return {"exp_c": sp.simplify(c * sp.diff(sp.log(eps_b), c))}
