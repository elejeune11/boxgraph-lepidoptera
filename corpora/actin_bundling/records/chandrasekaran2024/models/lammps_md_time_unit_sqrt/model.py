import math


def run(inputs, settings):
    sigma_m = inputs["sigma"] * 1e-9
    eps_J = inputs["eps_VV"] * 1e-21
    return {"tau_MD": sigma_m * math.sqrt(inputs["m"] / eps_J)}


def symbolic():
    import sympy as sp

    sigma, m, eps_VV = sp.symbols("sigma m eps_VV", positive=True)
    return {
        "tau_MD": (sigma * sp.Float("1e-9"))
        * sp.sqrt(m / (eps_VV * sp.Float("1e-21")))
    }
