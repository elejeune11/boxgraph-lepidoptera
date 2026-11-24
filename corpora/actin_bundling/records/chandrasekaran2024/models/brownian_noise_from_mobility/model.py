import math


def run(inputs, settings):
    # mu in um/(pN*s), kBT in pN*nm -> 1e-3 pN*um => D in um^2/s
    d_eff = inputs["mu"] * inputs["kBT"] * 1e-3
    return {"sigma_B": math.sqrt(2.0 * d_eff * inputs["dt"])}


def symbolic():
    import sympy as sp

    mu, kBT, dt = sp.symbols("mu kBT dt", positive=True)
    return {"sigma_B": sp.sqrt(2 * sp.Rational(1, 1000) * mu * kBT * dt)}
