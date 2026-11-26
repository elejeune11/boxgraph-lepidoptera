import math


def run(inputs, settings):
    return {"tau_mu": (1.0 / inputs["a_mu"]) * math.log(1.0 / inputs["r_mu"])}


def symbolic():
    import sympy as sp
    a_mu, r_mu = sp.symbols("a_mu r_mu", positive=True)
    return {"tau_mu": (1 / a_mu) * sp.log(1 / r_mu)}
