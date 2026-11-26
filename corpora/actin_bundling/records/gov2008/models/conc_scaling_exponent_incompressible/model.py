import math


def run(inputs, settings):
    c = inputs["c"]
    eps = inputs["eps"]
    # R'_{0,b} propto eps_b^(1/3) propto [c*exp(eps)/(1 + c*exp(eps))]^(1/3),
    # so d(ln R'_{0,b})/d(ln c) = (1/3) / (1 + c*exp(eps)), independent of unknown prefactor.
    return {"exp_c": (1.0 / 3.0) / (1.0 + c * math.exp(eps))}


def symbolic():
    import sympy as sp

    c, eps = sp.symbols("c eps", positive=True)
    eps_b = c * sp.exp(eps) / (1 + c * sp.exp(eps))
    return {"exp_c": sp.simplify(c * sp.diff(sp.log(eps_b ** sp.Rational(1, 3)), c))}
