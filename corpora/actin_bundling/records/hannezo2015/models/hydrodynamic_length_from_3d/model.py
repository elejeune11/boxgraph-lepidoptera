import math


def run(inputs, settings):
    eta = inputs["eta_3D"] * (inputs["h"] * 1e-9)
    return {"l_h": math.sqrt(eta / inputs["xi"]) * 1e6}


def symbolic():
    import sympy as sp

    eta_3D, h, xi = sp.symbols("eta_3D h xi", positive=True)
    return {"l_h": sp.sqrt(eta_3D * h * sp.Rational(1, 10**9) / xi) * 10**6}
