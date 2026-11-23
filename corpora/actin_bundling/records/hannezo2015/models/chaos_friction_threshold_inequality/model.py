import math


def run(inputs, settings):
    # From 2*pi*(D*tau*eta/zeta)^(1/4) < (eta/zeta)^(1/2) => zeta < eta / ((2*pi)^4 * D * tau)
    D_si = inputs["D"] * 1e-12
    return {"zeta_c": inputs["eta"] / (((2.0 * math.pi) ** 4) * D_si * inputs["tau"])}


def symbolic():
    import sympy as sp

    eta, D, tau = sp.symbols("eta D tau", positive=True)
    return {"zeta_c": eta / ((2 * sp.pi) ** 4 * D * sp.Rational(1, 10**12) * tau)}
