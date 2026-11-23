import math


def run(inputs, settings):
    D_si = inputs["D"] * 1e-12
    val = (math.sqrt(inputs["eta"] / inputs["tau"]) + math.sqrt(D_si * inputs["xi"])) ** 2
    return {"chi_c": val}


def symbolic():
    import sympy as sp

    eta, tau, D, xi = sp.symbols("eta tau D xi", positive=True)
    D_si = D * sp.Rational(1, 10**12)
    return {"chi_c": (sp.sqrt(eta / tau) + sp.sqrt(D_si * xi)) ** 2}
