import math


def run(inputs, settings):
    a = inputs["alpha"]
    phi = inputs["phi"]
    k = inputs["i"] * math.pi
    chi_i0 = phi + (k * k) / (a * a)
    chi_i2 = ((a * a * phi + k * k) ** 2) / (6.0 * a * a * k * k)
    return {
        "chi_i0": chi_i0,
        "chi_i0_over_phi": chi_i0 / phi,
        "chi_i2": chi_i2,
    }


def symbolic():
    import sympy as sp

    alpha, phi, i = sp.symbols("alpha phi i", positive=True)
    k = i * sp.pi
    chi_i0 = phi + k**2 / alpha**2
    chi_i2 = (alpha**2 * phi + k**2) ** 2 / (6 * alpha**2 * k**2)
    return {
        "chi_i0": chi_i0,
        "chi_i0_over_phi": chi_i0 / phi,
        "chi_i2": chi_i2,
    }
